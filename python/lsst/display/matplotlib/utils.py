import matplotlib.pyplot as pyplot
import numpy as np
import imageProc

try:
    _mpFigures
except NameError:
    _mpFigures = {0 : None}              # matplotlib (actually pyplot) figures
    eventHandlers = {}                  # event handlers for matplotlib figures

def getMpFigure(fig=None, clear=True):
    """Return a pyplot figure(); if fig is supplied save it and make it the default
    fig may also be a bool (make a new figure) or an int (return or make a figure (1-indexed;
    python-list style -n supported)
    """

    if not pyplot:
        raise RuntimeError("I am unable to plot as I failed to import matplotlib")

    if isinstance(fig, bool):       # we want a new one
        fig = len(_mpFigures) + 1    # matplotlib is 1-indexed

    if isinstance(fig, int):
        i = fig
        if i == 0:
            raise RuntimeError("I'm sorry, but matplotlib uses 1-indexed figures")
        if i < 0:
            try:
                i = sorted(_mpFigures.keys())[i] # simulate list's [-n] syntax
            except IndexError:
                if _mpFigures:
                    print >> sys.stderr, "Illegal index: %d" % i
                i = 1

        def lift(fig):
            fig.canvas._tkcanvas._root().lift() # == Tk's raise, but raise is a python reserved word

        if _mpFigures.has_key(i):
            try:
                lift(_mpFigures[i])
            except Exception, e:
                del _mpFigures[i]
                
        if not _mpFigures.has_key(i):
            for j in range(1, i):
                getMpFigure(j, clear=False)
                
            _mpFigures[i] = pyplot.figure()
            #
            # Modify pyplot.figure().show() to make it raise the plot too
            #
            def show(self, _show=_mpFigures[i].show):
                _show(self)
                try:
                    lift(self)
                except Exception, e:
                    pass
            # create a bound method
            import types
            _mpFigures[i].show = types.MethodType(show, _mpFigures[i], _mpFigures[i].__class__)

        fig = _mpFigures[i]

    if not fig:
        i = sorted(_mpFigures.keys())[0]
        if i > 0:
            fig = _mpFigures[i[-1]]
        else:
            fig = getMpFigure(1)

    if clear:
        fig.clf()

    pyplot.figure(fig.number)           # make it active

    return fig
            ims.append(np.load(os.path.join(dirName, "%s.npy" % f)))
    else:
        ims += [None, None]
        
    if eimage:
        ims.append(np.load(os.path.join(dirName, "eimage.npy")))
    else:
        ims.append(None)

    return ims

def mtv(im, I0=0, b=1, mask=None, isMask=False, alpha=None, fig=None, clear=True, evData=None):
    """Display an image, using an asinh stretch (softened by b)"""
    fig = getMpFigure(fig, clear=clear)

    try:
        mtv(im.image, I0=I0, b=b, fig=fig, evData=im)
        mtv(im.mask, isMask=True, alpha=alpha, fig=fig, clear=False)

        return
    except AttributeError:
        pass
    
    if isMask:
        if alpha is None:
            alpha = 0.7

        maskPlanes = imageProc.maskPlanes

        
        r = (im & (maskPlanes["BAD"] | maskPlanes["CR"])) != 0
        g = (im & (maskPlanes["INTRP"] | maskPlanes["SATUR"] | maskPlanes["EDGE"])) != 0
        b = (im & (maskPlanes["DETECTED"] | maskPlanes["EDGE"])) != 0

        alpha = alpha*np.ones_like(im)
        alpha[im == 0] = 0

        lim4 = np.dstack([r, g, b, alpha]).reshape([im.shape[0], im.shape[1], 4])
        pyplot.imshow(lim4, origin="lower", interpolation="nearest" , figure=fig)
    else:
        if b == 0:
            b = 1e-10
        ax = pyplot.imshow(np.arcsinh((im - I0)/b), origin='lower', interpolation='nearest',
                           cmap=pyplot.cm.gray, figure=fig)

        if mask is not None:
            mtv(mask, isMask=True, alpha=alpha, fig=fig)

    if evData:
        axes = fig.get_axes()[0]
        myText = axes.text(0.05, 1.05, 'Press "return" to show intensity here',
                           transform=axes.transAxes, va='top')

        global eventHandlers
        eventHandlers[fig] = EventHandler((evData, myText), fig)

class EventHandler(object):
    """A class to handle key strokes with matplotlib displays"""
    def __init__(self, data, fig):
        self.fig = fig

        im, text = data
        try:
            self.image = im.image
            self.mask = im.mask
        except AttributeError:
            self.image = im
            self.mask = None

        self.text = text

        self.cid = self.fig.canvas.mpl_connect('key_press_event', self)

    def __call__(self, ev):
        if ev.key != "\n":
            return

        if not (ev.xdata and ev.ydata):
            return
        x = np.clip(int(ev.xdata + 0.5), 0, self.image.shape[0])
        y = np.clip(int(ev.ydata + 0.5), 0, self.image.shape[1])
        str = "(%4d, %4d) %9.2f" % (x, y, self.image[y, x])
        if self.mask is not None:
            str += " 0x%02x" % (self.mask[y, x])

            mval = self.mask[y, x]
            for k, v in imageProc.maskPlanes.items():
                if mval & v:
                    str += " %s" % k


        self.text.set_text(str)
        self.fig.canvas.draw()
