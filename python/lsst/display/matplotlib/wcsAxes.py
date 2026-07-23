#
# LSST Data Management System
# Copyright 2026 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""Draw sky coordinate axes on a matplotlib Axes using AST.

The axes are drawn by `starlink.Ast.Plot` using the matplotlib grf
plugin, so arbitrary AST WCS transforms are supported with no FITS
approximation.
"""

__all__ = ["astFrameSetFromWcs", "WcsAxesManager"]

import astshim
import matplotlib.text
import starlink.Ast
from matplotlib.backend_bases import TimerBase
from starlink.Grf import grf_matplotlib


class _AstStringSource:
    """Present AST native serialization text to `starlink.Ast.Channel`.

    `starlink.Ast.Channel` reads via an object with an ``astsource``
    method returning one line per call, and `None` at end of input.

    Parameters
    ----------
    text : `str`
        AST native serialization of an AST object.
    """

    def __init__(self, text):
        self._lines = text.splitlines()
        self._pos = 0

    def astsource(self):
        if self._pos >= len(self._lines):
            return None
        line = self._lines[self._pos]
        self._pos += 1
        return line


def astFrameSetFromWcs(wcs):
    """Convert an afw SkyWcs to a starlink-pyast FrameSet.

    Parameters
    ----------
    wcs : `lsst.afw.geom.SkyWcs`
        The WCS to convert.

    Returns
    -------
    frameSet : `starlink.Ast.FrameSet`
        FrameSet whose base frame is LSST 0-based pixel coordinates
        (domain ``PIXELS``) and whose current frame is sky coordinates
        in radians (domain ``SKY``).

    Notes
    -----
    The conversion serializes the WCS's underlying `astshim.FrameDict`
    to AST native text and reads it back with pyast, so the result is
    exact; no FITS approximation is involved.
    """
    stream = astshim.StringStream()
    astshim.Channel(stream).write(wcs.getFrameDict())
    channel = starlink.Ast.Channel(_AstStringSource(stream.getSinkData()))
    return channel.read()


class WcsAxesManager:
    """Manage AST-drawn sky coordinate axes on a matplotlib Axes.

    All visible axis furniture (border, ticks, labels, optional grid)
    is drawn by AST; matplotlib's own axis furniture is hidden while
    this manager is active and restored by `remove`.

    The axes' data coordinates must be pixel coordinates in the LSST
    convention, matching the base frame of ``frameSet`` as produced by
    `astFrameSetFromWcs`: the center of the first pixel of the parent
    image is at (0, 0), and a subimage keeps its parent's coordinates
    (its pixel coordinates start at the subimage's XY0, not at zero).
    This differs from the FITS convention, in which the first pixel of
    the image is centered at (1, 1).

    Parameters
    ----------
    axes : `matplotlib.axes.Axes`
        The axes to draw on.
    frameSet : `starlink.Ast.FrameSet`
        Pixel-to-sky transform, base frame in axes data coordinates.
    grid : `bool`, optional
        Draw the full curvilinear coordinate grid across the image?
    useSexagesimal : `bool`, optional
        Format labels as e.g. HH:MM:SS.s rather than decimal degrees.
    extraOptions : `str`, optional
        Comma-separated AST Plot attribute settings appended after the
        defaults, so they take precedence (e.g.
        ``"Colour(grid)=2, Width(border)=2"``).

    Notes
    -----
    AST colour values are 1-based indices into the grf colour table
    (see `starlink.Grf.grf_matplotlib`), not colour names:
    1=default, 2=red, 3=green, 4=blue, 5=cyan, 6=magenta, 7=yellow,
    8=black, 9=dark grey, 10=grey, 11=light grey, 12=white.
    """

    def __init__(self, axes, frameSet, *,
                 grid=False, useSexagesimal=False, extraOptions=""):
        self._axes = axes
        self._frameSet = frameSet
        self._grid = grid
        self._useSexagesimal = useSexagesimal
        self._extraOptions = extraOptions
        self.artists = []

        self._savedVisibility = (
            axes.xaxis.get_visible(),
            axes.yaxis.get_visible(),
            {name: spine.get_visible() for name, spine in axes.spines.items()},
        )
        axes.xaxis.set_visible(False)
        axes.yaxis.set_visible(False)
        for spine in axes.spines.values():
            spine.set_visible(False)

        self._redrawing = False
        self._redrawTimer = None
        self._renderedSize = None
        self._callbackIds = [
            axes.callbacks.connect("xlim_changed", self._onLimitsChanged),
            axes.callbacks.connect("ylim_changed", self._onLimitsChanged),
        ]
        # A figure resize changes the axes size without changing the view
        # limits, so the labels must be re-spaced for the new geometry too.
        self._resizeCallbackId = axes.get_figure().canvas.mpl_connect(
            "resize_event", self._onResize)

        try:
            self.draw()
        except Exception:
            # Restore the matplotlib furniture so the axes remain
            # usable as ordinary pixel axes.
            self.remove()
            raise

    def _plotOptions(self):
        """Return the AST Plot options string for the current state."""
        options = ["DrawTitle=0"]
        options.append("Grid=1" if self._grid else "Grid=0")
        if not self._useSexagesimal:
            # Decimal degrees; ".*" lets the Plot choose the precision
            # needed to keep adjacent labels distinct.
            options.append("Format(1)=d.*")
            options.append("Format(2)=d.*")
        if self._extraOptions:
            options.append(self._extraOptions)
        return ", ".join(options)

    def draw(self):
        """Draw (or redraw) the AST axes for the current view limits."""
        self._removeArtists()

        xlo, xhi = self._axes.get_xlim()
        ylo, yhi = self._axes.get_ylim()
        box = (xlo, ylo, xhi, yhi)

        before = set(self._axes.lines) | set(self._axes.texts)
        try:
            plot = starlink.Ast.Plot(self._frameSet, box, box,
                                     grf_matplotlib(self._axes),
                                     self._plotOptions())
            plot.grid()
        except Exception:
            # grid() may have drawn some artists before failing.
            for artist in list(self._axes.lines) + list(self._axes.texts):
                if artist not in before:
                    artist.remove()
            raise
        self.artists = [a for a in list(self._axes.lines) + list(self._axes.texts)
                        if a not in before]

        # Axes.add_artist clips to the axes patch but the exterior
        # labels sit just outside it.
        for artist in self.artists:
            if isinstance(artist, matplotlib.text.Text):
                artist.set_clip_on(False)

        # Record the size the labels were spaced for, so a later resize
        # event that leaves the size unchanged can be ignored.
        self._renderedSize = self._canvasSize()

    def _canvasSize(self):
        """Return the current canvas size in device pixels."""
        return tuple(self._axes.get_figure().canvas.get_width_height())

    def _removeArtists(self):
        """Remove the AST artists from the axes."""
        for artist in self.artists:
            try:
                artist.remove()
            except (ValueError, NotImplementedError):
                # Already gone, e.g. the axes were cleared.
                pass
        self.artists = []

    # Idle time before redrawing after a view change.  Interactive
    # panning and zooming change the limits on every mouse-motion
    # event, and a full AST rebuild per event is far too slow.
    _redrawDelayMs = 200

    def _onLimitsChanged(self, axes):
        """Schedule a redraw for new view limits; matplotlib callback.

        The ``axes`` argument (the axes the callback fired for) is unused.
        """
        self._scheduleRedraw()

    def _onResize(self, event):
        """Schedule a redraw for a new figure size; matplotlib callback.

        The ``event`` argument (the resize event) is unused.  Interactive
        web backends (e.g. ipympl) emit ``resize_event`` repeatedly at the
        size the figure already has, and a redraw repaints the canvas,
        which prompts yet another such event; acting on every one produces
        an unbounded rebuild loop.  Rebuild only when the canvas has
        actually changed size since it was last drawn.
        """
        if self._canvasSize() == self._renderedSize:
            return
        self._scheduleRedraw()

    def _scheduleRedraw(self):
        """Debounce a rebuild of the AST axes after a view or size change.

        The redraw is debounced with a one-shot timer so that a stream of
        changes (interactive panning or zooming, the x/y pair from a single
        zoom, or a drag-resize) produces a single rebuild once the view
        settles.  Backends without a running event loop cannot fire timers,
        so there the redraw happens immediately.
        """
        if self._redrawing:
            return
        if self._redrawTimer is None:
            canvas = self._axes.get_figure().canvas
            timer = canvas.new_timer(interval=self._redrawDelayMs)
            if type(timer) is TimerBase:
                # The base class is a no-op that never fires.
                self._redrawNow()
                return
            timer.single_shot = True
            timer.add_callback(self._redrawNow)
            self._redrawTimer = timer
        self._redrawTimer.stop()
        self._redrawTimer.start()

    def _redrawNow(self):
        """Rebuild the AST axes for the current view and repaint."""
        if not self._callbackIds:
            # The manager was removed while a redraw was pending.
            return
        if self._redrawing:
            return
        self._redrawing = True
        try:
            self.draw()
        finally:
            self._redrawing = False
        self._axes.get_figure().canvas.draw_idle()

    def setSexagesimal(self, useSexagesimal):
        """Switch between sexagesimal and decimal degree labels.

        Parameters
        ----------
        useSexagesimal : `bool`
            Format labels as e.g. HH:MM:SS.s iff True.
        """
        useSexagesimal = bool(useSexagesimal)
        if useSexagesimal != self._useSexagesimal:
            self._useSexagesimal = useSexagesimal
            self.draw()

    def remove(self):
        """Remove the AST axes and restore matplotlib's axis furniture."""
        if self._redrawTimer is not None:
            self._redrawTimer.stop()
            self._redrawTimer = None
        for callbackId in self._callbackIds:
            self._axes.callbacks.disconnect(callbackId)
        self._callbackIds = []
        self._axes.get_figure().canvas.mpl_disconnect(self._resizeCallbackId)

        self._removeArtists()

        xVisible, yVisible, spines = self._savedVisibility
        self._axes.xaxis.set_visible(xVisible)
        self._axes.yaxis.set_visible(yVisible)
        for name, visible in spines.items():
            self._axes.spines[name].set_visible(visible)
