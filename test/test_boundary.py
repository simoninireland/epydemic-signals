# Tests of boundary signals
#
# Copyright (C) 2021--2022 Simon Dobson
#
# This file is part of epydemic-signals, an experiment in epidemic processes.
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

import unittest
from heapq import heappush, heappop
from epydemic_signals import *
from epydemic import SIR, StochasticDynamics, ProcessSequence, FixedNetwork
from epyc import Experiment
from networkx import Graph


class BoundarySignalTests(unittest.TestCase):

    def setUp(self):
        self._g = Graph()
        self._g.add_nodes_from([1, 2, 3, 4, 5, 6])
        self._g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])

        self._evs = [(1.0, SIR.INFECTED, (3, 1)),
                     (2.0, SIR.REMOVED, 1),
                     (3.0, SIR.INFECTED, (4, 3)),
                     (4.0, SIR.REMOVED, 3)]

        self._p = SIR()
        self._e = StochasticDynamics(self._p, FixedNetwork(self._g))
        self._e.setNetwork(self._g)
        self._params = dict({SIR.P_INFECTED: 0.0,
                             SIR.P_INFECT: 0.0,
                             SIR.P_REMOVE: 0.0})
        self._signal = Signal()
        self._generator = InfectionBoundarySignalGenerator(self._signal)
        self._generator.setExperiment(self._e)
        self._generator.setProcess(self._p)
        self._p.reset()
        self._p.build(self._params)
        self._p.setUp(self._params)
        self._p.changeCompartment(1, SIR.INFECTED)
        self._generator.setUp(self._g, self._params)

    def _playEventsTo(self, ft):
        '''Play all events up to and including time ft, against both the
        process and the signal generator.

        :param ft: the final event time
        :returns: the signal at time ft'''
        for (t, etype, e) in self._evs:
            if t <= ft:
                if etype == SIR.INFECTED:
                    self._p.infect(t, e)
                elif etype == SIR.REMOVED:
                    self._p.remove(t, e)
                self._generator.event(t, etype, e)
        return self._signal[ft]


    # ----------  Small tests ----------

    def testInitial(self):
        '''Test the original signal.'''
        s = self._signal[0.0]
        self.assertEqual(s[1], 2)
        self.assertEqual(s[2], 0)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 0)
        self.assertEqual(s[6], 0)

    def testAt1(self):
        '''Test the signal at t = 1.0.'''
        s = self._playEventsTo(1.0)
        self.assertEqual(s[1], 1)
        self.assertEqual(s[2], 0)
        self.assertEqual(s[3], 2)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 0)
        self.assertEqual(s[6], 0)

    def testAt2(self):
        '''Test the signal at t = 2.0.'''
        s = self._playEventsTo(2.0)
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 0)
        self.assertEqual(s[3], 2)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 0)
        self.assertEqual(s[6], 0)

    def testAt3(self):
        '''Test the signal at t = 3.0.'''
        s = self._playEventsTo(3.0)
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 0)
        self.assertEqual(s[3], 1)
        self.assertEqual(s[4], 3)
        self.assertEqual(s[5], 0)
        self.assertEqual(s[6], 0)

    def testAt4(self):
        '''Test the signal at t = 4.0.'''
        s = self._playEventsTo(4.0)
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 0)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 3)
        self.assertEqual(s[5], 0)
        self.assertEqual(s[6], 0)


if __name__ == '__main__':
    unittest.main()
