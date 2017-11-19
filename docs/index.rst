.. _lsst-display_matplotlib-name:

#######################################################
display.matplotlib An matplotlib backend for afwDisplay
#######################################################

.. _lsst-display_matplotlib-intro:

Introduction
============

We have a simple interface to display images on a variety of devices;  the original example
was ds9.  ``display.matplotlib`` provides a matplotlib backend; not only does this not require
you to have ds9 running (which might have involved a couple of ssh tunnels), but it can be inlined
in a jupyter notebook.

For example, you can display an ``afwImage.Exposure`` and then overlay the positions of detected objects.

.. _lsst-package-getting-started:

Getting Started
===============

Let's assume that you want to work from a notebook.  Start with

.. code-block:: python
   :emphasize-lines: 3-4

   %matplotlib notebook

   import lsst.afw.display as afwDisplay
   afwDisplay.setDefaultBackend("matplotlib")

(you might want something like ``%config InlineBackend.figure_format = 'retina'`` too)   

Then create an image and display it:

.. code-block:: python
   :emphasize-lines: 5,16-19

   import lsst.afw.geom as afwGeom
   import lsst.afw.image as afwImage
   import lsst.afw.math as afwMath

   disp = afwDisplay.Display()

   im = afwImage.MaskedImageF(100, 50)
   afwMath.randomGaussianImage(im.image, afwMath.Random())
   im.setXY0(10, 5)
   im[0, 0] = 20
   im.mask[0, 0] = im.mask.getPlaneBitMask("EDGE")
   im[-1, -1] = 30
   im.image[50, 22:24] = 40
   im.mask[50, 22] = 0x5

   disp.scale("asinh", "zscale", Q=8)
   title="test"

   disp.mtv(im, title=title)

.. figure:: /_static/example.png
   :name: fig-example
   :target: ../../_static/example.png

   The resulting displayed image.  You set the magenta pixel's mask to ``0x5``.  If you put the mouse
   on that pixel you'll see that its value is 40.0 and the set bits are ``BAD``, ``INTRP``.

Let's go ahead and add some overlays and zoom/pan to see what we did:

.. code-block:: python
   :emphasize-lines: 1-4,6-7

   disp.dot('o', 52, 23, ctype='RED')
   disp.dot('hello world', 52, 23)
   disp.dot('+', 55, 23, ctype='blue')
   disp.line([(50, 20), (58, 20), (58, 26), (50, 26), (50, 20)], ctype='cyan')

   disp.pan(55, 23)
   disp.zoom(4)

.. figure:: /_static/example-small.png
   :name: fig-example-small
   :target: ../../_static/example-small.png

   The resulting image.

.. _lsst-package-getting-using:

Using display_matplotlib
========================

All the standard ``afwDisplay`` commands are supported; see :py:class:`lsst.afw.display.Display`

Additionally, you can use ``matplotlib.pyplot`` commands;
the ``Figure`` can be retrieved as ``display._figure``.
The most useful of these commands is probably ``display._figure.savefig(fileName)``.

Interactive Callbacks
---------------------

The ``afwDisplay`` interface supports callbacks, but only for a subset of possible backends;  in particular
the inline ```notebook` and ``osx`` drivers **don't** support ``interact`` (but ``qt`` does).

Using Qt
--------

From a jupyter notebook, specify

.. code-block:: python

   %matplotlib qt
   %gui qt

to use the ``qt`` backend, which supports ``interact()``.
