#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2015 LSST Corporation.
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

#
# \file
# \brief Definitions to talk to matplotlib from python using the "afwDisplay"
#        interface

import math
import sys
import unicodedata

import matplotlib
import matplotlib.cm
import matplotlib.figure
import matplotlib.pyplot as pyplot
import matplotlib.cbook
import matplotlib.colors as mpColors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import numpy.ma as ma

import lsst.afw.display as afwDisplay
import lsst.afw.math as afwMath
import lsst.afw.display.rgb as afwRgb
import lsst.afw.display.interface as interface
import lsst.afw.display.virtualDevice as virtualDevice
import lsst.afw.display.ds9Regions as ds9Regions
import lsst.afw.image as afwImage

import lsst.afw.geom as afwGeom
import lsst.geom as geom

#
# Set the list of backends which support _getEvent and thus interact()
#
try:
    interactiveBackends
except NameError:
    # List of backends that support `interact`
    interactiveBackends = [
        "Qt4Agg",
        "Qt5Agg",
    ]

try:
    matplotlibCtypes
except NameError:
    matplotlibCtypes = {
        afwDisplay.GREEN: "#00FF00",
    }

    def mapCtype(ctype):
        """Map the ctype to a potentially different ctype

        Specifically, if matplotlibCtypes[ctype] exists, use it instead

        This is used e.g. to map "green" to a brighter shade
        """
        return matplotlibCtypes[ctype] if ctype in matplotlibCtypes else ctype


class DisplayImpl(virtualDevice.DisplayImpl):
    """Provide a matplotlib backend for afwDisplay

    Recommended backends in notebooks are:
      %matplotlib notebook
    or
      %matplotlib ipympl
    or
      %matplotlib qt
      %gui qt
    or
      %matplotlib inline
    or
      %matplotlib osx

    Apparently only qt supports Display.interact(); the list of interactive
    backends is given by lsst.display.matplotlib.interactiveBackends
    """
    def __init__(self, display, verbose=False,
                 interpretMaskBits=True, mtvOrigin=afwImage.PARENT, fastMaskDisplay=False,
                 reopenPlot=False, useSexagesimal=False, dpi=None, *args, **kwargs):
        """
        Initialise a matplotlib display

        @param fastMaskDisplay      If True only show the first bitplane that's
                                    set in each pixel
                                    (e.g. if (SATURATED & DETECTED)
                                    ignore DETECTED)
                                    Not really what we want, but a bit faster
        @param interpretMaskBits    Interpret the mask value under the cursor
        @param mtvOrigin            Display pixel coordinates with LOCAL origin
                                    (bottom left == 0,0 not XY0)
        @param reopenPlot           If true, close the plot before opening it.
                                    (useful with e.g. %ipympl)
        @param useSexagesimal       If True, display coordinates in sexagesimal
                                    E.g. hh:mm:ss.ss (default:False)
                                    May be changed by calling
                                          display.useSexagesimal()
        @param dpi                  Number of dpi (passed to pyplot.figure)

        The `frame` argument to `Display` may be a matplotlib figure; this
        permits code such as
           fig, axes = plt.subplots(1, 2)

           disp = afwDisplay.Display(fig)
           disp.scale('asinh', 'zscale', Q=0.5)

           for axis, exp in zip(axes, exps):
              fig.sca(axis)    # make axis active
              disp.mtv(exp)
        """
        try:
            fig_class = matplotlib.figure.FigureBase
        except AttributeError:
            fig_class = matplotlib.figure.Figure

        if isinstance(display.frame, fig_class):
            figure = display.frame
        else:
            figure = None

        virtualDevice.DisplayImpl.__init__(self, display, verbose)

        if reopenPlot:
            pyplot.close(display.frame)

        if figure is not None:
            self._figure = figure
        else:
            self._figure = pyplot.figure(display.frame, dpi=dpi)
            self._figure.clf()

        self._display = display
        self._maskTransparency = {None: 0.7}
        self._interpretMaskBits = interpretMaskBits  # interpret mask bits in mtv
        self._fastMaskDisplay = fastMaskDisplay
        self._useSexagesimal = [useSexagesimal]  # use an array so we can modify the value in format_coord
        self._mtvOrigin = mtvOrigin
        self._mappable_ax = None
        self._colorbar_ax = None
        self._image_colormap = matplotlib.cm.gray
        #
        self.__alpha = unicodedata.lookup("GREEK SMALL LETTER alpha")  # used in cursor display string
        self.__delta = unicodedata.lookup("GREEK SMALL LETTER delta")  # used in cursor display string
        #
        # Support self._scale()
        #
        self._scaleArgs = dict()
        self._normalize = None
        #
        # Support self._erase(), reporting pixel/mask values, and
        # zscale/minmax; set in mtv
        #
        self._i_setImage(None)

    def _close(self):
        """!Close the display, cleaning up any allocated resources"""
        self._image = None
        self._mask = None
        self._wcs = None
        self._figure.gca().format_coord = lambda x, y: None  # keeps a copy of _wcs

    def _show(self):
        """Put the plot at the top of the window stacking order"""

        try:
            self._figure.canvas._tkcanvas._root().lift()  # tk
        except AttributeError:
            pass

        try:
            self._figure.canvas.manager.window.raise_()  # os/x
        except AttributeError:
            pass

        try:
            self._figure.canvas.raise_()  # qt[45]
        except AttributeError:
            pass

    #
    # Extensions to the API
    #
    def savefig(self, *args, **kwargs):
        """Defer to figure.savefig()

        Parameters
        ----------
        args : `list`
          Passed through to figure.savefig()
        kwargs : `dict`
          Passed through to figure.savefig()
        """
        self._figure.savefig(*args, **kwargs)

    def show_colorbar(self, show=True, where="right", axSize="5%", axPad=None, **kwargs):
        """Show (or hide) the colour bar

        Parameters
        ----------
        show : `bool`
          Should I show the colour bar?
        where : `str`
          Location of colour bar: "right" or "bottom"
        axSize : `float` or `str`
          Size of axes to hold the colour bar; fraction of current x-size
        axPad : `float` or `str`
          Padding between axes and colour bar; fraction of current x-size
        args : `list`
          Passed through to colorbar()
        kwargs : `dict`
          Passed through to colorbar()

        We set the default padding to put the colourbar in a reasonable
        place for roughly square plots, but you may need to fiddle for
        plots with extreme axis ratios.

        You can only configure the colorbar when it isn't yet visible, but
        as you can easily remove it this is not in practice a difficulty.
        """
        if show:
            if self._mappable_ax:
                if self._colorbar_ax is None:
                    orientationDict = dict(right="vertical", bottom="horizontal")

                    mappable, ax = self._mappable_ax

                    if where in orientationDict:
                        orientation = orientationDict[where]
                    else:
                        print(f"Unknown location {where}; "
                              f"please use one of {', '.join(orientationDict.keys())}")

                    if axPad is None:
                        axPad = 0.1 if orientation == "vertical" else 0.3

                    divider = make_axes_locatable(ax)
                    self._colorbar_ax = divider.append_axes(where, size=axSize, pad=axPad)

                    self._figure.colorbar(mappable, cax=self._colorbar_ax, orientation=orientation, **kwargs)
                    self._figure.sca(ax)

        else:
            if self._colorbar_ax is not None:
                self._colorbar_ax.remove()
                self._colorbar_ax = None

    def useSexagesimal(self, useSexagesimal):
        """Control the formatting coordinates as HH:MM:SS.ss

        Parameters
        ----------
        useSexagesimal : `bool`
           Print coordinates as e.g. HH:MM:SS.ss iff True

        N.b. can also be set in Display's ctor
        """

        """Are we formatting coordinates as HH:MM:SS.ss?"""
        self._useSexagesimal[0] = useSexagesimal

    def wait(self, prompt="[c(ontinue) p(db)] :", allowPdb=True):
        """Wait for keyboard input

        Parameters
        ----------
        prompt : `str`
           The prompt string.
        allowPdb : `bool`
           If true, entering a 'p' or 'pdb' puts you into pdb

        Returns the string you entered

        Useful when plotting from a programme that exits such as a processCcd
        Any key except 'p' continues; 'p' puts you into pdb (unless
        allowPdb is False)
        """
        while True:
            s = input(prompt)
            if allowPdb and s in ("p", "pdb"):
                import pdb
                pdb.set_trace()
                continue

            return s
    #
    # Defined API
    #

    def _setMaskTransparency(self, transparency, maskplane):
        """Specify mask transparency (percent)"""

        self._maskTransparency[maskplane] = 0.01*transparency

    def _getMaskTransparency(self, maskplane=None):
        """Return the current mask transparency"""
        return self._maskTransparency[maskplane if maskplane in self._maskTransparency else None]

    def _mtv(self, image, mask=None, wcs=None, title=""):
        """Display an Image and/or Mask on a matplotlib display
        """
        title = str(title) if title else ""

        #
        # Save a reference to the image as it makes erase() easy and permits
        # printing cursor values and minmax/zscale stretches.  We also save XY0
        #
        self._i_setImage(image, mask, wcs)

        #  We need to know the pixel values to support e.g. 'zscale' and
        # 'minmax', so do the scaling now
        if self._scaleArgs.get('algorithm'):  # someone called self.scale()
            self._i_scale(self._scaleArgs['algorithm'], self._scaleArgs['minval'], self._scaleArgs['maxval'],
                          self._scaleArgs['unit'], *self._scaleArgs['args'], **self._scaleArgs['kwargs'])

        ax = self._figure.gca()
        ax.cla()

        self._i_mtv(image, wcs, title, False)

        if mask:
            self._i_mtv(mask, wcs, title, True)

        self.show_colorbar()

        if title:
            ax.set_title(title)

        self._title = title

        def format_coord(x, y, wcs=self._wcs, x0=self._xy0[0], y0=self._xy0[1],
                         origin=afwImage.PARENT, bbox=self._image.getBBox(afwImage.PARENT),
                         _useSexagesimal=self._useSexagesimal):

            fmt = '(%1.2f, %1.2f)'
            if self._mtvOrigin == afwImage.PARENT:
                msg = fmt % (x, y)
            else:
                msg = (fmt + "L") % (x - x0, y - y0)

            col = int(x + 0.5)
            row = int(y + 0.5)
            if bbox.contains(geom.PointI(col, row)):
                if wcs is not None:
                    raDec = wcs.pixelToSky(x, y)
                    ra = raDec[0].asDegrees()
                    dec = raDec[1].asDegrees()

                    if _useSexagesimal[0]:
                        from astropy import units as u
                        from astropy.coordinates import Angle as apAngle

                        kwargs = dict(sep=':', pad=True, precision=2)
                        ra = apAngle(ra*u.deg).to_string(unit=u.hour, **kwargs)
                        dec = apAngle(dec*u.deg).to_string(unit=u.deg, **kwargs)
                    else:
                        ra = "%9.4f" % ra
                        dec = "%9.4f" % dec

                    msg += r" (%s, %s): (%s, %s)" % (self.__alpha, self.__delta, ra, dec)

                msg += ' %1.3f' % (self._image[col, row])
                if self._mask:
                    val = self._mask[col, row]
                    if self._interpretMaskBits:
                        msg += " [%s]" % self._mask.interpret(val)
                    else:
                        msg += " 0x%x" % val

            return msg

        ax.format_coord = format_coord
        # Stop images from reporting their value as we've already
        # printed it nicely
        for a in ax.get_images():
            a.get_cursor_data = lambda ev: None  # disabled

        # using tight_layout() is too tight and clips the axes
        self._figure.canvas.draw_idle()

    def _i_mtv(self, data, wcs, title, isMask):
        """Internal routine to display an Image or Mask on a DS9 display"""

        title = str(title) if title else ""
        dataArr = data.getArray()

        if isMask:
            maskPlanes = data.getMaskPlaneDict()
            nMaskPlanes = max(maskPlanes.values()) + 1

            planes = {}                      # build inverse dictionary
            for key in maskPlanes:
                planes[maskPlanes[key]] = key

            planeList = range(nMaskPlanes)

            maskArr = np.zeros_like(dataArr, dtype=np.int32)

            colorNames = ['black']
            colorGenerator = self.display.maskColorGenerator(omitBW=True)
            for p in planeList:
                color = self.display.getMaskPlaneColor(planes[p]) if p in planes else None

                if not color:            # none was specified
                    color = next(colorGenerator)
                elif color.lower() == afwDisplay.IGNORE:
                    color = 'black'     # we'll set alpha = 0 anyway

                colorNames.append(color)
            #
            # Convert those colours to RGBA so we can have per-mask-plane
            # transparency and build a colour map
            #
            # Pixels equal to 0 don't get set (as no bits are set), so leave
            # them transparent and start our colours at [1] --
            # hence "i + 1" below
            #
            colors = mpColors.to_rgba_array(colorNames)
            alphaChannel = 3            # the alpha channel; the A in RGBA
            colors[0][alphaChannel] = 0.0      # it's black anyway
            for i, p in enumerate(planeList):
                if colorNames[i + 1] == 'black':
                    alpha = 0.0
                else:
                    alpha = 1 - self._getMaskTransparency(planes[p] if p in planes else None)

                colors[i + 1][alphaChannel] = alpha

            cmap = mpColors.ListedColormap(colors)
            norm = mpColors.NoNorm()
        else:
            cmap = self._image_colormap
            norm = self._normalize

        ax = self._figure.gca()
        bbox = data.getBBox()
        extent = (bbox.getBeginX() - 0.5, bbox.getEndX() - 0.5,
                  bbox.getBeginY() - 0.5, bbox.getEndY() - 0.5)

        with matplotlib.rc_context(dict(interactive=False)):
            if isMask:
                for i, p in reversed(list(enumerate(planeList))):
                    if colors[i + 1][alphaChannel] == 0:  # colors[0] is reserved
                        continue

                    bitIsSet = (dataArr & (1 << p)) != 0
                    if bitIsSet.sum() == 0:
                        continue

                    maskArr[bitIsSet] = i + 1  # + 1 as we set colorNames[0] to black

                    if not self._fastMaskDisplay:  # we draw each bitplane separately
                        ax.imshow(maskArr, origin='lower', interpolation='nearest',
                                  extent=extent, cmap=cmap, norm=norm)
                        maskArr[:] = 0

                if self._fastMaskDisplay:  # we only draw the lowest bitplane
                    ax.imshow(maskArr, origin='lower', interpolation='nearest',
                              extent=extent, cmap=cmap, norm=norm)
            else:
                # If we're playing with subplots and have reset the axis
                # the cached colorbar axis belongs to the old one, so set
                # it to None
                if self._mappable_ax and self._mappable_ax[1] != self._figure.gca():
                    self._colorbar_ax = None

                mappable = ax.imshow(dataArr, origin='lower', interpolation='nearest',
                                     extent=extent, cmap=cmap, norm=norm)
                self._mappable_ax = (mappable, ax)

        self._figure.canvas.draw_idle()

    def _i_setImage(self, image, mask=None, wcs=None):
        """Save the current image, mask, wcs, and XY0"""
        self._image = image
        self._mask = mask
        self._wcs = wcs
        self._xy0 = self._image.getXY0() if self._image else (0, 0)

        self._zoomfac = None
        if self._image is None:
            self._width, self._height = 0, 0
        else:
            self._width, self._height = self._image.getDimensions()

        self._xcen = 0.5*self._width
        self._ycen = 0.5*self._height

    def _setImageColormap(self, cmap):
        """Set the colormap used for the image

        cmap should be either the name of an attribute of matplotlib.cm or an
        mpColors.Colormap (e.g. "gray" or matplotlib.cm.gray)

        """
        if not isinstance(cmap, mpColors.Colormap):
            try:
                cmap = matplotlib.colormaps[cmap]
            except AttributeError:
                cmap = getattr(matplotlib.cm, cmap)

        self._image_colormap = cmap

    #
    # Graphics commands
    #

    def _buffer(self, enable=True):
        if enable:
            pyplot.ioff()
        else:
            pyplot.ion()
            self._figure.show()

    def _flush(self):
        pass

    def _erase(self):
        """Erase the display"""

        for axis in self._figure.axes:
            axis.lines = []
            axis.texts = []

        self._figure.canvas.draw_idle()

    def _dot(self, symb, c, r, size, ctype,
             fontFamily="helvetica", textAngle=None):
        """Draw a symbol at (col,row) = (c,r) [0-based coordinates]
    Possible values are:
            +                        Draw a +
            x                        Draw an x
            *                        Draw a *
            o                        Draw a circle
            @:Mxx,Mxy,Myy            Draw an ellipse with moments
                                     (Mxx, Mxy, Myy) (argument size is ignored)
            An afwGeom.ellipses.Axes Draw the ellipse (argument size is
                                     ignored)

    Any other value is interpreted as a string to be drawn. Strings obey the
    fontFamily (which may be extended with other characteristics, e.g.
    "times bold italic".  Text will be drawn rotated by textAngle
    (textAngle is ignored otherwise).
    """
        if not ctype:
            ctype = afwDisplay.GREEN

        axis = self._figure.gca()
        x0, y0 = self._xy0

        if isinstance(symb, afwGeom.ellipses.Axes):
            from matplotlib.patches import Ellipse

            # Following matplotlib.patches.Ellipse documentation 'width' and
            # 'height' are diameters while 'angle' is rotation in degrees
            # (anti-clockwise)
            axis.add_artist(Ellipse((c + x0, r + y0), height=2*symb.getA(), width=2*symb.getB(),
                                    angle=90.0 + math.degrees(symb.getTheta()),
                                    edgecolor=mapCtype(ctype), facecolor='none'))
        elif symb == 'o':
            from matplotlib.patches import CirclePolygon as Circle

            axis.add_artist(Circle((c + x0, r + y0), radius=size, color=mapCtype(ctype), fill=False))
        else:
            from matplotlib.lines import Line2D

            for ds9Cmd in ds9Regions.dot(symb, c + x0, r + y0, size, fontFamily="helvetica", textAngle=None):
                tmp = ds9Cmd.split('#')
                cmd = tmp.pop(0).split()

                cmd, args = cmd[0], cmd[1:]

                if cmd == "line":
                    args = np.array(args).astype(float) - 1.0

                    x = np.empty(len(args)//2)
                    y = np.empty_like(x)
                    i = np.arange(len(args), dtype=int)
                    x = args[i%2 == 0]
                    y = args[i%2 == 1]

                    axis.add_line(Line2D(x, y, color=mapCtype(ctype)))
                elif cmd == "text":
                    x, y = np.array(args[0:2]).astype(float) - 1.0
                    axis.text(x, y, symb, color=mapCtype(ctype),
                              horizontalalignment='center', verticalalignment='center')
                else:
                    raise RuntimeError(ds9Cmd)

    def _drawLines(self, points, ctype):
        """Connect the points, a list of (col,row)
        Ctype is the name of a colour (e.g. 'red')"""

        from matplotlib.lines import Line2D

        if not ctype:
            ctype = afwDisplay.GREEN

        points = np.array(points)
        x = points[:, 0] + self._xy0[0]
        y = points[:, 1] + self._xy0[1]

        self._figure.gca().add_line(Line2D(x, y, color=mapCtype(ctype)))

    def _scale(self, algorithm, minval, maxval, unit, *args, **kwargs):
        """
        Set gray scale

        N.b.  Supports extra arguments:
        @param maskedPixels  List of names of mask bits to ignore
                             E.g. ["BAD", "INTERP"].
                             A single name is also supported
        """
        self._scaleArgs['algorithm'] = algorithm
        self._scaleArgs['minval'] = minval
        self._scaleArgs['maxval'] = maxval
        self._scaleArgs['unit'] = unit
        self._scaleArgs['args'] = args
        self._scaleArgs['kwargs'] = kwargs

        try:
            self._i_scale(algorithm, minval, maxval, unit, *args, **kwargs)
        except (AttributeError, RuntimeError):
            # Unable to access self._image; we'll try again when we run mtv
            pass

    def _i_scale(self, algorithm, minval, maxval, unit, *args, **kwargs):

        maskedPixels = kwargs.get("maskedPixels", [])
        if isinstance(maskedPixels, str):
            maskedPixels = [maskedPixels]
        bitmask = afwImage.Mask.getPlaneBitMask(maskedPixels)

        sctrl = afwMath.StatisticsControl()
        sctrl.setAndMask(bitmask)

        if minval == "minmax":
            if self._image is None:
                raise RuntimeError("You may only use minmax if an image is loaded into the display")

            mi = afwImage.makeMaskedImage(self._image, self._mask)
            stats = afwMath.makeStatistics(mi, afwMath.MIN | afwMath.MAX, sctrl)
            minval = stats.getValue(afwMath.MIN)
            maxval = stats.getValue(afwMath.MAX)
        elif minval == "zscale":
            if bitmask:
                print("scale(..., 'zscale', maskedPixels=...) is not yet implemented")

        if algorithm is None:
            self._normalize = None
        elif algorithm == "asinh":
            if minval == "zscale":
                if self._image is None:
                    raise RuntimeError("You may only use zscale if an image is loaded into the display")

                self._normalize = AsinhZScaleNormalize(image=self._image, Q=kwargs.get("Q", 8.0))
            else:
                self._normalize = AsinhNormalize(minimum=minval,
                                                 dataRange=maxval - minval, Q=kwargs.get("Q", 8.0))
        elif algorithm == "linear":
            if minval == "zscale":
                if self._image is None:
                    raise RuntimeError("You may only use zscale if an image is loaded into the display")

                self._normalize = ZScaleNormalize(image=self._image,
                                                  nSamples=kwargs.get("nSamples", 1000),
                                                  contrast=kwargs.get("contrast", 0.25))
            else:
                self._normalize = LinearNormalize(minimum=minval, maximum=maxval)
        else:
            raise RuntimeError("Unsupported stretch algorithm \"%s\"" % algorithm)
    #
    # Zoom and Pan
    #

    def _zoom(self, zoomfac):
        """Zoom by specified amount"""

        self._zoomfac = zoomfac

        if zoomfac is None:
            return

        x0, y0 = self._xy0

        size = min(self._width, self._height)
        if size < self._zoomfac:        # avoid min == max
            size = self._zoomfac
        xmin, xmax = self._xcen + x0 + size/self._zoomfac*np.array([-1, 1])
        ymin, ymax = self._ycen + y0 + size/self._zoomfac*np.array([-1, 1])

        ax = self._figure.gca()

        tb = self._figure.canvas.toolbar
        if tb is not None:              # It's None for e.g. %matplotlib inline in jupyter
            tb.push_current()           # save the current zoom in the view stack

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal', 'datalim')

        self._figure.canvas.draw_idle()

    def _pan(self, colc, rowc):
        """Pan to (colc, rowc)"""

        self._xcen = colc
        self._ycen = rowc

        self._zoom(self._zoomfac)

    def _getEvent(self, timeout=-1):
        """Listen for a key press, returning (key, x, y)"""

        if timeout < 0:
            timeout = 24*3600           # -1 generates complaints in QTimer::singleShot.  A day is a long time

        mpBackend = matplotlib.get_backend()
        if mpBackend not in interactiveBackends:
            print("The %s matplotlib backend doesn't support display._getEvent()" %
                  (matplotlib.get_backend(),), file=sys.stderr)
            return interface.Event('q')

        event = None

        # We set up a blocking event loop. On receipt of a keypress, the
        # callback records the event and unblocks the loop.

        def recordKeypress(keypress):
            """Matplotlib callback to record keypress and unblock"""
            nonlocal event
            event = interface.Event(keypress.key, keypress.xdata, keypress.ydata)
            self._figure.canvas.stop_event_loop()

        conn = self._figure.canvas.mpl_connect("key_press_event", recordKeypress)
        try:
            self._figure.canvas.start_event_loop(timeout=timeout)  # Blocks on keypress
        finally:
            self._figure.canvas.mpl_disconnect(conn)
        return event


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class Normalize(mpColors.Normalize):
    """Class to support stretches for mtv()"""

    def __call__(self, value, clip=None):
        """
        Return a MaskedArray with value mapped to [0, 255]

        @param value Input pixel value or array to be mapped
        """
        if isinstance(value, np.ndarray):
            data = value
        else:
            data = value.data

        data = data - self.mapping.minimum[0]
        return ma.array(data*self.mapping.mapIntensityToUint8(data)/255.0)


class AsinhNormalize(Normalize):
    """Provide an asinh stretch for mtv()"""
    def __init__(self, minimum=0, dataRange=1, Q=8):
        """Initialise an object able to carry out an asinh mapping

        @param minimum   Minimum pixel value (default: 0)
        @param dataRange Range of values for stretch if Q=0; roughly the
                         linear part (default: 1)
        @param Q Softening parameter (default: 8)

        See Lupton et al., PASP 116, 133
        """
        # The object used to perform the desired mapping
        self.mapping = afwRgb.AsinhMapping(minimum, dataRange, Q)

        vmin, vmax = self._getMinMaxQ()[0:2]
        if vmax*Q > vmin:
            vmax *= Q
        super().__init__(vmin, vmax)

    def _getMinMaxQ(self):
        """Return an asinh mapping's minimum and maximum value, and Q

        Regrettably this information is not preserved by AsinhMapping
        so we have to reverse engineer it
        """

        frac = 0.1                      # magic number in AsinhMapping
        Q = np.sinh((frac*self.mapping._uint8Max)/self.mapping._slope)/frac
        dataRange = Q/self.mapping._soften

        vmin = self.mapping.minimum[0]
        return vmin, vmin + dataRange, Q


class AsinhZScaleNormalize(AsinhNormalize):
    """Provide an asinh stretch using zscale to set limits for mtv()"""
    def __init__(self, image=None, Q=8):
        """Initialise an object able to carry out an asinh mapping

        @param image  image to use estimate minimum and dataRange using zscale
                      (see AsinhNormalize)
        @param Q Softening parameter (default: 8)

        See Lupton et al., PASP 116, 133
        """

        # The object used to perform the desired mapping
        self.mapping = afwRgb.AsinhZScaleMapping(image, Q)

        vmin, vmax = self._getMinMaxQ()[0:2]
        # n.b. super() would call AsinhNormalize,
        # and I want to pass min/max to the baseclass
        Normalize.__init__(self, vmin, vmax)


class ZScaleNormalize(Normalize):
    """Provide a zscale stretch for mtv()"""
    def __init__(self, image=None, nSamples=1000, contrast=0.25):
        """Initialise an object able to carry out a zscale mapping

        @param image to be used to estimate the stretch
        @param nSamples Number of data points to use (default: 1000)
        @param contrast Control the range of pixels to display around the
                        median (default: 0.25)
        """

        # The object used to perform the desired mapping
        self.mapping = afwRgb.ZScaleMapping(image, nSamples, contrast)

        super().__init__(self.mapping.minimum[0], self.mapping.maximum)


class LinearNormalize(Normalize):
    """Provide a linear stretch for mtv()"""
    def __init__(self, minimum=0, maximum=1):
        """Initialise an object able to carry out a linear mapping

        @param minimum  Minimum value to display
        @param maximum  Maximum value to display
        """
        # The object used to perform the desired mapping
        self.mapping = afwRgb.LinearMapping(minimum, maximum)

        super().__init__(self.mapping.minimum[0], self.mapping.maximum)
