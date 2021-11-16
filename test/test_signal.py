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
        self._signal = Signal(self._g)

    def test_NegativeTime(self):
        '''Test we can't set the time negative.'''
        with self.assertRaises(Exception):
            self._signal.setTime(-1)


if __name__ == '__main__':
    unittest.main()
