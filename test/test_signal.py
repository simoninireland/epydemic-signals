# Tests of basic signals
#
# Copyright (C) 2021 Simon Dobson
#
# This file is part of epydemic-signals, an experiment in epidemics processes.
#
# epydemic-signals is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# epydemic-signals is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with epydemic-signals. If not, see <http://www.gnu.org/licenses/gpl.html>.

from epydemic_signals import *
import epyc
import unittest
import networkx


class SignalTests(unittest.TestCase):

    def setUp(self):
        self._g = networkx.Graph()
        self._g.add_nodes_from([1, 2, 3, 4, 5, 6])
        self._g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])
        self._signal = Signal(None)

    def testEmpty(self):
        '''Test the properties of an empty signal.'''
        self.assertCountEqual(self._signal.transitions(), [])
        self.assertCountEqual(self._signal[0].keys(), [])
        with self.assertRaises(ValueError):
            self._signal.getBounds()

    def testOnePointSignal(self):
        '''Test the behaviour of a signal with one entry.'''
        d = self._signal[0]
        d['a'] = 10
        s = self._signal[1]

        self.assertCountEqual(self._signal.transitions(), [0])
        self.assertEqual(self._signal.getBounds(), (10, 10))

    def testTwoPointSignal(self):
        '''Test the behaviour of a signal with two points.'''
        d = self._signal[0]
        d['a'] = 10
        d['b'] = 30
        d = self._signal[1]
        d['b'] = 20

        self.assertCountEqual(self._signal.transitions(), [0, 1])
        self.assertEqual(self._signal.getBounds(), (10, 30))


if __name__ == '__main__':
    unittest.main()
