# Tests of progress signals
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

import unittest
from epydemic_signals import *
from epydemic import SIR, StochasticDynamics, ProcessSequence, FixedNetwork
from epyc import Experiment
from networkx import Graph, fast_gnp_random_graph, grid_graph, convert_node_labels_to_integers


class ProgressSignalTests(unittest.TestCase):

    def setUp(self):
        self._g = Graph()
        self._g.add_nodes_from([1, 2, 3, 4, 5, 6])
        self._g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])

        self._evs = [(0.0, SIR.INFECTED, 1),
                     (1.0, SIR.INFECTED, 3),
                     (2.0, SIR.REMOVED, 1),
                     (3.0, SIR.INFECTED, 4),
                     (4.0, SIR.REMOVED, 3)]

        self._signal = SIRProgressSignal(self._g, self._evs)


    # ----------  Small tests ----------

    def testBase(self):
        '''Test the base signal is correct.'''
        s = self._signal[0.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 1)
        self.assertEqual(s[4], 2)
        self.assertEqual(s[5], 3)
        self.assertEqual(s[6], 3)

    def testSlightlyBeyondBase(self):
        '''Test that times before the first transition stay like base.'''
        s = self._signal[0.2]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 1)
        self.assertEqual(s[4], 2)
        self.assertEqual(s[5], 3)
        self.assertEqual(s[6], 3)

    def testTransitionAt1(self):
        '''Test t=1.0.'''
        s = self._signal[1.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 1)
        self.assertEqual(s[5], 2)
        self.assertEqual(s[6], 2)

    def testTransitionAt2(self):
        '''Test t=2.0.'''
        s = self._signal[2.0]
        self.assertEqual(s[1], -1)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 1)
        self.assertEqual(s[5], 2)
        self.assertEqual(s[6], 2)

    def testTransitionAt3(self):
        '''Test t=3.0.'''
        s = self._signal[3.0]
        self.assertEqual(s[1], -1)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 1)
        self.assertEqual(s[6], 1)

    def testTransitionAt4(self):
        '''Test t=4.'''
        s = self._signal[4.0]
        self.assertEqual(s[1], -2)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], -1)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 1)
        self.assertEqual(s[6], 1)

    def testLate(self):
        '''Test the signal doesn't change after the last transition.'''
        s = self._signal[6.0]
        self.assertEqual(s[1], -2)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], -1)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 1)
        self.assertEqual(s[6], 1)

    def testBackAndForward(self):
        '''Test that the signal backs-up correctly.'''
        s = self._signal[2.0]
        s = self._signal[1.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 1)
        self.assertEqual(s[5], 2)
        self.assertEqual(s[6], 2)


    # ---------- Soak test ----------

    def testInvariantsAdd(self):
        '''Test invariants when adding infected nodes.'''
        g = Graph()
        g.add_nodes_from([1, 2, 3, 4, 5, 6])
        g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])
        evs = [(0.0, SIR.INFECTED, 1),
               (1.0, SIR.INFECTED, 3),
               (2.0, SIR.INFECTED, 5),
               (3.0, SIR.INFECTED, 4),
               (4.0, SIR.INFECTED, 2),
               (5.0, SIR.INFECTED, 6)]
        self.checkInvariants(g, evs)

    def testInvariantsRemove(self):
        '''Test invariants when adding and then removing infected nodes.'''
        g = Graph()
        g.add_nodes_from([1, 2, 3, 4, 5, 6])
        g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])
        evs = [(0.0, SIR.INFECTED, 1),
               (1.0, SIR.INFECTED, 3),
               (2.0, SIR.INFECTED, 5),
               (3.0, SIR.INFECTED, 4),
               (4.0, SIR.INFECTED, 2),
               (5.0, SIR.INFECTED, 6),
               (6.0, SIR.REMOVED, 4),
               (7.0, SIR.REMOVED, 5),
               (8.0, SIR.REMOVED, 1),
               (9.0, SIR.REMOVED, 3),
               (10.0, SIR.REMOVED, 2),
               (11.0, SIR.REMOVED, 6)]
        self.checkInvariants(g, evs)

    def testInvariantsLatticeManual(self):
        '''Test invariants on a larger network.'''
        g = convert_node_labels_to_integers(grid_graph(dim=(6, 6)), first_label=1)
        evs = [(0.0, SIR.INFECTED, 21),
               (1.0, SIR.INFECTED, 22)]
        self.checkInvariants(g, evs)

    def testInvariantsLattice(self):
        '''Test invariants for an epidemic on a lattice.'''
        g = convert_node_labels_to_integers(grid_graph(dim=(10, 10)), first_label=1)
        pInfect = 0.5
        pRemove = 0.01

        params = dict()
        params[SIR.P_INFECT] = pInfect
        params[SIR.P_REMOVE] = pRemove

        sir = OneInfectionSIR()
        p = ProcessSequence([sir, HittingHealingTimes(sir)])
        e = StochasticDynamics(p, FixedNetwork(g))
        rc = e.set(params).run(fatal=True)
        evs = HittingHealingTimes.timeline(rc)
        self.checkInvariants(g, evs)

    def testInvariantsLarge(self):
        '''Test invariants as we run a larger epidemic.'''
        N = 2000
        kmean = 50
        g = fast_gnp_random_graph(N, kmean / N)
        pInfect = 0.01
        pRemove = 0.001

        params = dict()
        params[SIR.P_INFECT] = pInfect
        params[SIR.P_REMOVE] = pRemove

        sir = OneInfectionSIR()
        p = ProcessSequence([sir, HittingHealingTimes(sir)])
        e = StochasticDynamics(p, FixedNetwork(g))
        rc = e.set(params).run()
        evs = HittingHealingTimes.timeline(rc)
        self.checkInvariants(g, evs)

    def checkInvariants(self, g, evs):
        signal = SIRProgressSignal(g, evs)
        ns = list(g.nodes())
        susceptibles = ns.copy()
        infecteds = set()
        removeds = set()
        for (t, e, s) in evs:
            #print(t, e, s)
            sig = signal[t]

            # advance the sets
            self.assertIn(e, [SIR.INFECTED, SIR.REMOVED])
            if e == SIR.INFECTED:
                susceptibles.remove(s)
                infecteds.add(s)
            elif e == SIR.REMOVED:
                infecteds.remove(s)
                removeds.add(s)

            # test all nodes have an entry in the signal
            self.assertCountEqual(ns, sig.keys())

            # test that all infecteds are zeros
            for n in infecteds:
                #print(f'inf check {n}')
                self.assertEqual(sig[n], 0)

            # check all susceptibles are the right distance from the boundary
            self.checkSusceptibles(g, sig, susceptibles, infecteds, removeds)

    def checkSusceptibles(self, g, sig, susceptibles, infecteds, removeds):
        ss = susceptibles.copy()
        while len(ss) > 0:
            n = ss.pop(0)
            #print(f'sus check {n}')
            d = sig[n]
            # all neighbours should have distances differing by at most one
            # from us (if they're susceptibles), or be infecteds (in which case
            # our distance should be 1), or be removeds
            for m in g.neighbors(n):
                if m in susceptibles:
                    #print(n, m, d, sig[m])
                    self.assertTrue(abs(sig[m] - d) <= 1)
                elif m in infecteds:
                    self.assertEqual(d, 1)


if __name__ == '__main__':
    unittest.main()
