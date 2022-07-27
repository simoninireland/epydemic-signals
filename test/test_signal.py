# Tests of basic signals
#
# Copyright (C) 2021--2022 Simon Dobson
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
import unittest
import networkx


class SignalTests(unittest.TestCase):

    def setUp(self):
        self._g = networkx.Graph()
        self._g.add_nodes_from([1, 2, 3, 4, 5, 6])
        self._g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])
        self._signal = Signal()

    def testEmpty(self):
        '''Test the properties of an empty signal.'''
        self.assertCountEqual(self._signal.transitions(), [])
        self.assertCountEqual(self._signal[0].keys(), [])
        self.assertEqual(len(self._signal.values()), 0)

    def testOnePointSignal(self):
        '''Test the behaviour of a signal with one entry.'''
        d = self._signal[0]
        d['a'] = 10
        s = self._signal[1]

        self.assertCountEqual(self._signal.transitions(), [0])
        self.assertCountEqual(self._signal.values(), [10])

    def testTwoPointSignal(self):
        '''Test the behaviour of a signal with two points.'''
        d = self._signal[0]
        d['a'] = 10
        d['b'] = 30
        d = self._signal[1]
        d['b'] = 20

        self.assertCountEqual(self._signal.transitions(), [0, 1])
        self.assertCountEqual(self._signal.values(), [10, 20, 30])

    def testTimeSeries(self):
        '''Test we can convert a node signal to a set of time series.'''
        self._signal.setNetwork(self._g)
        ts = [0, 1, 2, 3]
        for t in ts:
            s_t = self._signal[t]
            for n in self._g.nodes():
                s_t[n] = t * n

        # extract the time series
        tss = self._signal.toTimeSeries()
        self.assertCountEqual(tss.keys(), self._g.nodes())
        for n in self._g.nodes():
            ts_n = tss[n]
            self.assertEqual(len(ts_n), len(self._signal.transitions()))
            self.assertCountEqual(ts_n, [t * n for t in ts])

    def testUpdates(self):
        '''Test we can convert a node signal to update lists.'''
        self._signal.setNetwork(self._g)
        ts = [0, 1, 2, 3]
        for t in ts:
            s_t = self._signal[t]
            for n in self._g.nodes():
                s_t[n] = t * n

        # extract the update lists
        (times, nodes, values) = self._signal.toUpdates()
        self.assertEqual(len(times), len(nodes))
        self.assertEqual(len(times), len(values))

        # compute what the updates should look like
        for (t, n, v) in zip(times, nodes, values):
            self.assertEqual(v, t * n)
            self.assertEqual(self._signal[t][n], v)

    def testNonUpdates(self):
        '''Test that values unchanged at a timestep don't get updates.'''
        self._signal.setNetwork(self._g)
        ts = [0, 1, 2, 3]
        for t in ts:
            s_t = self._signal[t]
            for n in self._g.nodes():
                if n == 1:
                    s_t[n] = 0       # node 1 has constant signal 0
                else:
                    s_t[n] = t * n

        # compute what the updates should look like
        (times, nodes, values) = self._signal.toUpdates()
        seenOne = False
        for (t, n, v) in zip(times, nodes, values):
            if n == 1:
                # should only appear in the first update
                self.assertEqual(t, 0)
                self.assertEqual(v, 0)
                seenOne = True              # ...but we need to see it once
            else:
                self.assertEqual(v, t * n)
                self.assertEqual(self._signal[t][n], v)
        self.assertTrue(seenOne)


if __name__ == '__main__':
    unittest.main()
