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

##
## \file
## \brief Definitions to talk to matplotlib from python using the "afwDisplay" interface

from __future__ import absolute_import, division, print_function

import math
import os
import re
import sys
import time
import unicodedata
import warnings

import matplotlib.pyplot as pyplot
import matplotlib.cbook
import matplotlib.colors as mpColors
from matplotlib.blocking_input import BlockingInput

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

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Set the list of backends which support _getEvent and thus interact()
#
try:
    interactiveBackends
except NameError:
    ## List of backends that support `interact`
    interactiveBackends = [
        "Qt4Agg",
    ]

class DisplayImpl(virtualDevice.DisplayImpl):
    """Provide a matplotlib backend for afwDisplay

    Recommended backends in notebooks are:
      %matplotlib notebook
    or
      %matplotlib qt
      %gui qt
    or
      %matplotlib inline
    or 
      %matplotlib osx

    Apparently only qt supports Display.interact(); the list of interactive backends
    is given by lsst.display.matplotlib.interactiveBackends
    """
    def __init__(self, display, verbose=False, interpretMaskBits=True, mtvOrigin=afwImage.PARENT,
                 *args, **kwargs):
        """
        Initialise a matplotlib display

        @param interpretMaskBits    Interpret the mask value under the cursor
        @param mtvOrigin            Display pixel coordinates with LOCAL origin
                                    (bottom left == 0,0 not XY0)
        """
        virtualDevice.DisplayImpl.__init__(self, display, verbose)

        self._figure = pyplot.figure(display.frame)
        self._display = display
        self._maskTransparency = {None : 0.7}
        self._interpretMaskBits = interpretMaskBits # interpret mask bits in mtv
        self._mtvOrigin = mtvOrigin
        self._mappable = None
        #
        self.__alpha = unicodedata.lookup("GREEK SMALL LETTER alpha") # used in cursor display string
        self.__delta = unicodedata.lookup("GREEK SMALL LETTER delta") # used in cursor display string
        #
        # Support self._scale()
        #
        self._scaleArgs = dict()
        self._normalize = None
        #
        # Support self._erase(), reporting pixel/mask values, and zscale/minmax; set in mtv
        #
        self._i_setImage(None)
        #
        # Ignore warnings due to BlockingKeyInput
        #
        if not verbose:
            warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

    def _close(self):
        """!Close the display, cleaning up any allocated resources"""
        self._image = None
        self._mask = None
        self._wcs = None
        self._figure.gca().format_coord = None # keeps a copy of _wcs

    def _show(self):
        """Put the plot at the top of the window stacking order"""

        try:
            self._figure.canvas._tkcanvas._root().lift() # tk
        except AttributeError:
            pass

        try:
            self._figure.canvas.manager.window.raise_() # os/x
        except AttributeError:
            pass

        try:
            self._figure.canvas.raise_() # qt[45]
        except AttributeError:
            pass

    #
    # Extensions to the API
    #
    def show_colorbar(self, show=True):
        """Show (or hide) the colour bar"""
        if show:
            if self._mappable:
                self._figure.colorbar(self._mappable)
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
        # Save a reference to the image as it makes erase() easy and permits printing cursor values
        # and minmax/zscale stretches.  We also save XY0
        #
        self._i_setImage(image, mask, wcs)
        
        # We need to know the pixel values to support e.g. 'zscale' and 'minmax', so do the scaling now
        if self._scaleArgs.get('algorithm'): # someone called self.scale()
            self._i_scale(self._scaleArgs['algorithm'], self._scaleArgs['minval'], self._scaleArgs['maxval'],
                          self._scaleArgs['unit'], *self._scaleArgs['args'], **self._scaleArgs['kwargs'])

        self._figure.clf()              # calling erase() calls _mtv

        self._i_mtv(image, wcs, title, False)
        ax = self._figure.gca()

        if mask:
            self._i_mtv(mask, wcs, title, True)
            
        if title:
            ax.set_title(title)

        self._title = title
        #
        def format_coord(x, y, wcs=self._wcs, x0=self._xy0[0], y0=self._xy0[1],
                         origin=afwImage.PARENT, bbox=self._image.getBBox(afwImage.PARENT)):

            fmt = '(%1.2f, %1.2f)' 
            if self._mtvOrigin == afwImage.PARENT:
                msg = fmt % (x, y)
            else:
                msg = (fmt + "L") % (x - x0, y - y0)

            col = int(x + 0.5)
            row = int(y + 0.5)
            if bbox.contains(afwGeom.PointI(col, row)):
                if wcs is not None:
                    ra, dec = wcs.pixelToSky(x, y)
                    msg += r" (%s, %s): (%9.4f, %9.4f)" % (self.__alpha, self.__delta, ra, dec)

                msg += ' %1.3f' % (self._image[col, row])
                if self._mask:
                    val = self._mask[col, row]
                    if self._interpretMaskBits:
                        msg += " [%s]" % self._mask.interpret(val)
                    else:
                        msg += " 0x%x" % val

            return msg

        ax.format_coord = format_coord
        # Stop images from reporting their value as we've already printed it nicely
        from matplotlib.image import AxesImage
        for a in ax.mouseover_set:
            if isinstance(a, AxesImage):
                a.get_cursor_data = lambda ev: None # disabled

        self._figure.tight_layout()
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
            # Set the maskArr image to be an index into our colour map (cmap; see below)
            #
            for i, p in enumerate(planeList):
                color = colorNames[i]
                maskArr[(dataArr & (1 << p)) != 0] += i + 1 # + 1 as we set colorNames[0] to black

            #
            # Convert those colours to RGBA so we can have per-mask-plane transparency
            # and build a colour map
            #
            colors = mpColors.to_rgba_array(colorNames)
            colors[0][3] = 0.0          # it's black anyway
            for i, p in enumerate(planeList):
                if colorNames[i + 1] == 'black':
                    alpha = 0.0
                else:
                    alpha = 1 - self._getMaskTransparency(planes[p] if p in planes else None)

                colors[i + 1][3] = alpha

            dataArr = maskArr
            cmap = mpColors.ListedColormap(colors)
            norm = mpColors.NoNorm()
        else:
            cmap = pyplot.cm.gray
            norm = self._normalize

        ax = self._figure.gca()
        bbox = data.getBBox()
        mappable = ax.imshow(dataArr, origin='lower', interpolation='nearest',
                             extent=(bbox.getBeginX() - 0.5, bbox.getEndX() - 0.5,
                                     bbox.getBeginY() - 0.5, bbox.getEndY() - 0.5),
                             cmap=cmap, norm=norm)

        if not isMask:
            self._mappable = mappable

        self._figure.canvas.draw_idle()

    def _i_setImage(self, image, mask=None, wcs=None):
        """Save the current image, mask, wcs, and XY0"""
        self._image = image
        self._mask = mask
        self._wcs = wcs
        self._xy0 = self._image.getXY0() if self._image else (0, 0)

        self._zoomfac = 1.0
        if self._image is None:
            self._width, self._height = 0, 0
        else:
            self._width, self._height = self._image.getDimensions()

        self._xcen = 0.5*self._width
        self._ycen = 0.5*self._height

    #
    # Graphics commands
    #
    def _buffer(self, enable=True):
        pass

    def _flush(self):
        pass

    def _erase(self):
        """Erase the display"""
        #
        # Rather than erase only the glyphs we'll redraw the image.
        #
        # This isn't a great solution.
        #
        self._figure.clf()

        if self._image:
            zoomfac = self._zoomfac
            xcen = self._xcen
            ycen = self._ycen

            self._mtv(self._image, mask=self._mask, wcs=self._wcs, title=self._title)

            self._xcen = xcen
            self._ycen = ycen
            self._zoom(zoomfac)

        self._figure.canvas.draw_idle()

    def _dot(self, symb, c, r, size, ctype,
             fontFamily="helvetica", textAngle=None):
        """Draw a symbol at (col,row) = (c,r) [0-based coordinates]
    Possible values are:
            +                Draw a +
            x                Draw an x
            *                Draw a *
            o                Draw a circle
            @:Mxx,Mxy,Myy    Draw an ellipse with moments (Mxx, Mxy, Myy) (argument size is ignored)
            An object derived from afwGeom.ellipses.BaseCore Draw the ellipse (argument size is ignored)
    Any other value is interpreted as a string to be drawn. Strings obey the fontFamily (which may be extended
    with other characteristics, e.g. "times bold italic".  Text will be drawn rotated by textAngle (textAngle is
    ignored otherwise).

    N.b. objects derived from BaseCore include Axes and Quadrupole.
    """
        if not ctype:
            ctype = afwDisplay.GREEN

        axis = self._figure.gca()
        x0, y0 = self._xy0
        
        if isinstance(symb, afwGeom.ellipses.BaseCore):
            from matplotlib.patches import Ellipse

            # Following matplotlib.patches.Ellipse documentation 'width' and 'height' are diameters while 
            # 'angle' is rotation in degrees (anti-clockwise)
            axis.add_artist(Ellipse((c + x0, r + y0), height=2*symb.getA(), width=2*symb.getB(),
                                    angle=90.+math.degrees(symb.getTheta()), 
                                    edgecolor=ctype, facecolor='none'))
        elif symb == 'o':
            from matplotlib.patches import CirclePolygon as Circle

            axis.add_artist(Circle((c + x0, r + y0), radius=size, color=ctype, fill=False))
        else:
            from matplotlib.lines import Line2D

            for ds9Cmd in ds9Regions.dot(symb, c + x0, r + y0, size, fontFamily="helvetica", textAngle=None):
                tmp = ds9Cmd.split('#')
                cmd = tmp.pop(0).split()
                comment = tmp.pop(0) if tmp else ""

                cmd, args = cmd[0], cmd[1:]
                
                if cmd == "line":
                    args = np.array(args).astype(float) - 1.0

                    x = np.empty(len(args)//2)
                    y = np.empty_like(x)
                    i = np.arange(len(args), dtype=int)
                    x = args[i%2 == 0]
                    y = args[i%2 == 1]

                    axis.add_line(Line2D(x, y, color=ctype))
                elif cmd == "text":
                    x, y = np.array(args[0:2]).astype(float) - 1.0
                    axis.text(x, y, symb, color=ctype,
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

        self._figure.gca().add_line(Line2D(x, y, color=ctype))
    #
    # Set gray scale
    #
    def _scale(self, algorithm, minval, maxval, unit, *args, **kwargs):
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
        if minval == "minmax":
            if self._image is None:
                raise RuntimeError("You may only use minmax if an image is loaded into the display")

            stats = afwMath.makeStatistics(self._image, afwMath.MIN | afwMath.MAX)
            minval = stats.getValue(afwMath.MIN)
            maxval = stats.getValue(afwMath.MAX)

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

        x0, y0 = self._xy0
        x1, y1 = x0 + self._width, y0 + self._height

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
            
        self._figure.canvas.draw_idle()

    def _pan(self, colc, rowc):
        """Pan to (colc, rowc)"""

        self._xcen = colc
        self._ycen = rowc

        self._zoom(self._zoomfac)        

    def _getEvent(self, timeout=-1):
        """Listen for a key press, returning (key, x, y)"""

        mpBackend = matplotlib.get_backend()
        if mpBackend not in interactiveBackends:
            print("The %s matplotlib backend doesn't support display._getEvent()" % 
                  (matplotlib.get_backend(),), file=sys.stderr)
            return interface.Event('q')
        
        blocking_input = BlockingKeyInput(self._figure)
        return blocking_input(timeout=timeout)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class BlockingKeyInput(BlockingInput):
    """
    Callable class to retrieve a single keyboard click
    """
    def __init__(self, fig):
        """Create a BlockingKeyInput

        \param fig The figure to monitor for keyboard events
        """
        BlockingInput.__init__(self, fig=fig, eventslist=('key_press_event',))

    def post_event(self):
        """
        Return the event containing the key and (x, y)
        """
        try:
            event = self.events[-1]
        except IndexError:
            ## details of the event to pass back to the display
            self.ev = None
        else:
            self.ev = interface.Event(event.key, event.xdata, event.ydata)

    def __call__(self, timeout=-1):
        """
        Blocking call to retrieve a single key click
        Returns key or None if timeout
        """
        self.ev = None

        BlockingInput.__call__(self, n=1, timeout=timeout)

        return self.ev

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
        @param dataRange Range of values for stretch if Q=0; roughly the linear part (default: 1)
        @param Q Softening parameter (default: 8)

        See Lupton et al., PASP 116, 133
        """
        Normalize.__init__(self)

        ## The object used to perform the desired mapping
        self.mapping = afwRgb.AsinhMapping(minimum, dataRange, Q)

class AsinhZScaleNormalize(Normalize):
    """Provide an asinh stretch using zscale to set limits for mtv()"""
    def __init__(self, image=None, Q=8):
        """Initialise an object able to carry out an asinh mapping

        @param image  image to use estimate minimum and dataRange using zscale (see AsinhNormalize)
        @param Q Softening parameter (default: 8)

        See Lupton et al., PASP 116, 133
        """
        Normalize.__init__(self)

        ## The object used to perform the desired mapping
        self.mapping = afwRgb.AsinhZScaleMapping(image, Q)

class ZScaleNormalize(Normalize):
    """Provide a zscale stretch for mtv()"""
    def __init__(self, image=None, nSamples=1000, contrast=0.25):
        """Initialise an object able to carry out a zscale mapping

        @param image to be used to estimate the stretch
        @param nSamples Number of data points to use (default: 1000)
        @param contrast Control the range of pixels to display around the median (default: 0.25)
        """

        Normalize.__init__(self)

        ## The object used to perform the desired mapping
        self.mapping = afwRgb.ZScaleMapping(image, nSamples, contrast)

class LinearNormalize(Normalize):
    """Provide a linear stretch for mtv()"""
    def __init__(self, minimum=0, maximum=1):
        """Initialise an object able to carry out a linear mapping

        @param minimum  Minimum value to display
        @param maximum  Maximum value to display
        """

        Normalize.__init__(self)

        ## The object used to perform the desired mapping
        self.mapping = afwRgb.LinearMapping(minimum, maximum)
