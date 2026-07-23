#
# Copyright 2008-2017  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
"""
Tests for display_matplotlib
"""
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure  # noqa: E402
import matplotlib.text  # noqa: E402
from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402

import lsst.utils.tests  # noqa: E402
import lsst.geom  # noqa: E402
import lsst.afw.display as afwDisplay  # noqa: E402
import lsst.afw.geom as afwGeom  # noqa: E402
import lsst.afw.image as afwImage  # noqa: E402


def makeTestWcs():
    """Return a rotated TAN SkyWcs for testing."""
    return afwGeom.makeSkyWcs(
        crpix=lsst.geom.Point2D(100, 75),
        crval=lsst.geom.SpherePoint(30.0, -45.0, lsst.geom.degrees),
        cdMatrix=afwGeom.makeCdMatrix(scale=3.0*lsst.geom.arcseconds,
                                      orientation=30*lsst.geom.degrees),
    )


def makeFineTestWcs():
    """Return a fine-scale TAN SkyWcs whose decimal labels are wide.

    At this scale the numeric labels nearly fill the gap between ticks,
    so any mismatch between the axes geometry AST measures against and
    the geometry it is finally drawn at makes them overlap.
    """
    return afwGeom.makeSkyWcs(
        crpix=lsst.geom.Point2D(150, 150),
        crval=lsst.geom.SpherePoint(53.0876, -27.5207, lsst.geom.degrees),
        cdMatrix=afwGeom.makeCdMatrix(scale=0.2*lsst.geom.arcseconds,
                                      orientation=0*lsst.geom.degrees),
    )


class AstFrameSetTestCase(lsst.utils.tests.TestCase):
    """Tests for the SkyWcs to pyast FrameSet conversion."""

    def testRoundTrip(self):
        """The pyast FrameSet must reproduce pixelToSky exactly."""
        from lsst.display.matplotlib.wcsAxes import astFrameSetFromWcs

        wcs = makeTestWcs()
        frameSet = astFrameSetFromWcs(wcs)
        for x, y in [(0.0, 0.0), (100.0, 75.0), (199.5, 149.5), (-20.0, 300.0)]:
            sky = wcs.pixelToSky(x, y)
            result = frameSet.tran([[x], [y]])
            self.assertEqual(result[0][0], sky[0].asRadians())
            self.assertEqual(result[1][0], sky[1].asRadians())


def makeWcsAxes(**kwargs):
    """Create an Agg figure with pixel-coordinate axes and a manager.

    Returns the (axes, manager) pair.
    """
    from lsst.display.matplotlib.wcsAxes import WcsAxesManager, astFrameSetFromWcs

    fig = matplotlib.figure.Figure(figsize=(8, 6))
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    ax.set_xlim(-0.5, 199.5)
    ax.set_ylim(-0.5, 149.5)
    wcs = makeTestWcs()
    manager = WcsAxesManager(ax, astFrameSetFromWcs(wcs), **kwargs)
    return ax, manager


class WcsAxesManagerTestCase(lsst.utils.tests.TestCase):
    """Tests for WcsAxesManager drawing and artist lifecycle."""

    def testDrawCreatesArtists(self):
        ax, manager = makeWcsAxes()
        self.assertGreater(len(manager.artists), 0)
        # All tracked artists really are in the axes.
        axesArtists = set(ax.lines) | set(ax.texts)
        for artist in manager.artists:
            self.assertIn(artist, axesArtists)
        # Sky axis labels come from the AST SkyFrame.
        texts = [a.get_text() for a in manager.artists
                 if isinstance(a, matplotlib.text.Text)]
        self.assertIn("Right ascension", texts)
        self.assertIn("Declination", texts)

    def testFurnitureHidden(self):
        ax, manager = makeWcsAxes()
        self.assertFalse(ax.xaxis.get_visible())
        self.assertFalse(ax.yaxis.get_visible())
        for spine in ax.spines.values():
            self.assertFalse(spine.get_visible())

    def testTextNotClipped(self):
        ax, manager = makeWcsAxes()
        for artist in manager.artists:
            if isinstance(artist, matplotlib.text.Text):
                self.assertFalse(artist.get_clip_on())

    def testDecimalLabelsByDefault(self):
        ax, manager = makeWcsAxes()
        tickLabels = [a.get_text() for a in manager.artists
                      if isinstance(a, matplotlib.text.Text) and
                      a.get_text() not in ("Right ascension", "Declination")]
        self.assertGreater(len(tickLabels), 0)
        for label in tickLabels:
            self.assertNotIn(":", label)
        # AST chooses the precision; the labels must be decimal degrees
        # and adjacent labels must be distinct.
        self.assertTrue(any("." in label for label in tickLabels))
        self.assertEqual(len(tickLabels), len(set(tickLabels)))

    def testSexagesimalLabels(self):
        ax, manager = makeWcsAxes(useSexagesimal=True)
        tickLabels = [a.get_text() for a in manager.artists
                      if isinstance(a, matplotlib.text.Text)]
        self.assertTrue(any(":" in label for label in tickLabels))

    def testRedrawOnLimitChange(self):
        ax, manager = makeWcsAxes()
        oldArtists = list(manager.artists)
        ax.set_xlim(50, 150)
        ax.set_ylim(40, 110)
        self.assertGreater(len(manager.artists), 0)
        self.assertNotEqual(manager.artists, oldArtists)
        # The old artists must be gone from the axes.
        axesArtists = set(ax.lines) | set(ax.texts)
        for artist in oldArtists:
            self.assertNotIn(artist, axesArtists)

    def testDebouncedRedraw(self):
        """With a working timer, drag events coalesce into one redraw."""
        from matplotlib.backend_bases import TimerBase

        starts = []

        class FakeTimer(TimerBase):
            def _timer_start(self):
                starts.append(self)

            def _timer_stop(self):
                pass

        ax, manager = makeWcsAxes()
        canvas = ax.get_figure().canvas
        canvas.new_timer = lambda interval=None: FakeTimer(interval=interval)

        oldArtists = list(manager.artists)
        for i in range(5):
            ax.set_xlim(-0.5 + i, 199.5 + i)
            ax.set_ylim(-0.5 + i, 149.5 + i)
        # No redraw yet: the debounce timer is pending.
        self.assertEqual(manager.artists, oldArtists)
        self.assertGreater(len(starts), 0)
        # Fire the pending timer: exactly one rebuild for the new view.
        manager._redrawTimer._on_timer()
        self.assertNotEqual(manager.artists, oldArtists)
        axesArtists = set(ax.lines) | set(ax.texts)
        for artist in oldArtists:
            self.assertNotIn(artist, axesArtists)

    def testNoRedrawAfterRemove(self):
        ax, manager = makeWcsAxes()
        manager.remove()
        ax.set_xlim(50, 150)
        self.assertEqual(manager.artists, [])
        self.assertEqual(len(ax.lines), 0)

    def testSetSexagesimal(self):
        ax, manager = makeWcsAxes()

        def tickLabels():
            return [a.get_text() for a in manager.artists
                    if isinstance(a, matplotlib.text.Text)]

        self.assertFalse(any(":" in label for label in tickLabels()))
        manager.setSexagesimal(True)
        self.assertTrue(any(":" in label for label in tickLabels()))
        manager.setSexagesimal(False)
        self.assertFalse(any(":" in label for label in tickLabels()))

    def testConstructorFailureRestoresFurniture(self):
        from lsst.display.matplotlib.wcsAxes import WcsAxesManager, astFrameSetFromWcs

        fig = matplotlib.figure.Figure(figsize=(8, 6))
        FigureCanvasAgg(fig)
        ax = fig.add_subplot(111)
        ax.set_xlim(-0.5, 199.5)
        ax.set_ylim(-0.5, 149.5)
        with self.assertRaises(Exception):
            WcsAxesManager(ax, astFrameSetFromWcs(makeTestWcs()),
                           extraOptions="NoSuchAttr=1")
        self.assertTrue(ax.xaxis.get_visible())
        self.assertTrue(ax.yaxis.get_visible())
        for spine in ax.spines.values():
            self.assertTrue(spine.get_visible())
        self.assertEqual(len(ax.lines) + len(ax.texts), 0)

    def testRemove(self):
        ax, manager = makeWcsAxes()
        manager.remove()
        self.assertEqual(manager.artists, [])
        self.assertEqual(len(ax.lines), 0)
        self.assertEqual(len(ax.texts), 0)
        self.assertTrue(ax.xaxis.get_visible())
        self.assertTrue(ax.yaxis.get_visible())
        for spine in ax.spines.values():
            self.assertTrue(spine.get_visible())


class WcsAxesDisplayTestCase(lsst.utils.tests.TestCase):
    """Tests for WCS axes at the afwDisplay level."""

    frame = 100  # keep test figures separate from any other test

    def setUp(self):
        afwDisplay.setDefaultBackend("matplotlib")
        # A distinct frame per test so each test gets a fresh figure.
        WcsAxesDisplayTestCase.frame += 1
        self.frame = WcsAxesDisplayTestCase.frame

    def tearDown(self):
        afwDisplay.delAllDisplays()

    @staticmethod
    def makeExposure():
        exposure = afwImage.ExposureF(200, 150)
        exposure.image.array[:] = 1.0
        exposure.setWcs(makeTestWcs())
        return exposure

    def testWcsAxesByDefault(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        self.assertIn(ax, impl._wcsAxesManagers)
        self.assertGreater(len(impl._wcsAxesManagers[ax].artists), 0)
        self.assertFalse(ax.xaxis.get_visible())

    def testUseWcsAxesFalse(self):
        display = afwDisplay.Display(frame=self.frame, useWcsAxes=False)
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        self.assertEqual(impl._wcsAxesManagers, {})
        self.assertTrue(ax.xaxis.get_visible())

    def testNoWcsMeansPixelAxes(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(afwImage.ImageF(200, 150))
        impl = display._impl
        ax = impl._figure.gca()
        self.assertEqual(impl._wcsAxesManagers, {})
        self.assertTrue(ax.xaxis.get_visible())

    def testWcsGridOption(self):
        display = afwDisplay.Display(frame=self.frame, wcsGrid=True)
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        manager = impl._wcsAxesManagers[ax]
        self.assertIn("Grid=1", manager._plotOptions())

    def testLabelsDoNotOverlapColorbar(self):
        """Numeric labels must be spaced for the axes as finally drawn.

        The colorbar resizes the image axes, and AST must space the
        labels for that final size rather than the full-width axes it
        starts from, or wide decimal labels overlap.
        """
        exposure = afwImage.ExposureF(300, 300)
        exposure.image.array[:] = 1.0
        exposure.setWcs(makeFineTestWcs())

        display = afwDisplay.Display(frame=self.frame)
        display.mtv(exposure)
        fig = display._impl._figure
        ax = fig.gca()
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        # Numeric tick labels drawn horizontally (the bottom, RA axis).
        boxes = []
        for artist in ax.texts:
            text = artist.get_text()
            if text in ("Right ascension", "Declination"):
                continue
            if any(c.isdigit() for c in text) and artist.get_rotation() % 180 == 0:
                boxes.append(artist.get_window_extent(renderer))
        # Keep the row nearest the bottom of the axes.
        yBottom = min(box.y0 for box in boxes)
        bottom = sorted((box for box in boxes if box.y0 < yBottom + 20),
                        key=lambda box: box.x0)
        self.assertGreater(len(bottom), 1)
        for left, right in zip(bottom, bottom[1:]):
            self.assertLessEqual(left.x1, right.x0,
                                 "adjacent RA labels overlap")

    def testZoomRedraws(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        manager = impl._wcsAxesManagers[ax]
        oldArtists = list(manager.artists)
        display.zoom(4)
        display.pan(100, 75)
        self.assertGreater(len(manager.artists), 0)
        self.assertNotEqual(manager.artists, oldArtists)

    def testErasePreservesWcsAxes(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        manager = impl._wcsAxesManagers[ax]
        nAstLines = sum(1 for a in manager.artists if a in set(ax.lines))
        display.dot("+", 100, 75)
        self.assertGreater(len(ax.lines), nAstLines)
        display.erase()
        self.assertEqual(len(ax.lines), nAstLines)
        axesArtists = set(ax.lines) | set(ax.texts)
        for artist in manager.artists:
            self.assertIn(artist, axesArtists)

    def testFallbackRestoresPixelAxes(self):
        display = afwDisplay.Display(frame=self.frame, astPlotOptions="NoSuchAttr=1")
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        self.assertEqual(impl._wcsAxesManagers, {})
        self.assertTrue(ax.xaxis.get_visible())
        for spine in ax.spines.values():
            self.assertTrue(spine.get_visible())

    def testEraseRemovesPatches(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(self.makeExposure())
        ax = display._impl._figure.gca()
        display.dot("o", 100, 75, size=10)
        self.assertEqual(len(ax.patches), 1)
        display.erase()
        self.assertEqual(len(ax.patches), 0)

    def testRepeatedMtvDoesNotLeak(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        nTexts = len(ax.texts)
        display.mtv(self.makeExposure())
        ax = impl._figure.gca()
        self.assertEqual(len(impl._wcsAxesManagers), 1)
        self.assertEqual(len(ax.texts), nTexts)

    def testUseSexagesimalRedraws(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(self.makeExposure())
        impl = display._impl
        ax = impl._figure.gca()
        manager = impl._wcsAxesManagers[ax]

        def hasColon():
            return any(":" in a.get_text() for a in manager.artists
                       if isinstance(a, matplotlib.text.Text))

        self.assertFalse(hasColon())
        display.useSexagesimal(True)
        self.assertTrue(hasColon())

    def testCloseReleasesManagers(self):
        display = afwDisplay.Display(frame=self.frame)
        display.mtv(self.makeExposure())
        impl = display._impl
        impl._close()
        self.assertEqual(impl._wcsAxesManagers, {})


class DisplayMatplotlibTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testSetAfwDisplayBackend(self):
        """Set the default backend to this package"""
        afwDisplay.setDefaultBackend("matplotlib")

    def testSetImageColormap(self):
        """This is a stand-in for an eventual testcase for changing image
        colormap.

        The basic outline should look something like:

        afwDisplay.setDefaultBackend("matplotlib")
        display = afwDisplay.Display()
        display.setImageColormap('viridis')
        assert display._image_colormap == 'viridis'
        """
        pass


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
