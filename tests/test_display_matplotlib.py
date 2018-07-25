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

import lsst.utils.tests
import lsst.afw.display as afwDisplay

class DisplayMatplotlibTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testSetAfwDisplayBackend(self):
        """Set the default backend to this package"""
        afwDisplay.setDefaultBackend("matplotlib")

    def testSetImageColormap(self):
        """This is a stand-in for an eventual testcase for changing image colormap
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
