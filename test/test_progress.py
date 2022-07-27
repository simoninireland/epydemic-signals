# Tests of progress signals
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
from networkx import Graph, fast_gnp_random_graph, grid_graph, convert_node_labels_to_integers


class ProgressSignalTests(unittest.TestCase):

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
        self._e.setNetwork((self._g))
        self._params = dict({SIR.P_INFECTED: 0.0,
                             SIR.P_INFECT: 0.0,
                             SIR.P_REMOVE: 0.0})
        self._signal = Signal()
        self._generator = SIRProgressSignalGenerator(self._signal)
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

    def testBase(self):
        '''Test the base signal is correct.'''
        self._playEventsTo(0.0)
        s = self._signal[0.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 1)
        self.assertEqual(s[4], 2)
        self.assertEqual(s[5], 3)
        self.assertEqual(s[6], 3)

    def testSlightlyBeyondBase(self):
        '''Test that times before the first transition stay like base.'''
        self._playEventsTo(0.0)
        s = self._signal[0.2]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 1)
        self.assertEqual(s[4], 2)
        self.assertEqual(s[5], 3)
        self.assertEqual(s[6], 3)

    def testTransitionAt1(self):
        '''Test t=1.0.'''
        self._playEventsTo(1.0)
        s = self._signal[1.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 1)
        self.assertEqual(s[5], 2)
        self.assertEqual(s[6], 2)

    def testTransitionAt2(self):
        '''Test t=2.0.'''
        self._playEventsTo(2.0)
        s = self._signal[2.0]
        self.assertEqual(s[1], -1)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 1)
        self.assertEqual(s[5], 2)
        self.assertEqual(s[6], 2)

    def testTransitionAt3(self):
        '''Test t=3.0.'''
        self._playEventsTo(3.0)
        s = self._signal[3.0]
        self.assertEqual(s[1], -1)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 1)
        self.assertEqual(s[6], 1)

    def testTransitionAt4(self):
        '''Test t=4.'''
        self._playEventsTo(4.0)
        s = self._signal[4.0]
        self.assertEqual(s[1], -2)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], -1)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 1)
        self.assertEqual(s[6], 1)

    def testLate(self):
        '''Test the signal doesn't change after the last transition.'''
        self._playEventsTo(6.0)
        s = self._signal[6.0]
        self.assertEqual(s[1], -2)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], -1)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 1)
        self.assertEqual(s[6], 1)

    def testBackAndForward(self):
        '''Test that the signal backs-up correctly.'''
        self._playEventsTo(4.0)
        s = self._signal[2.0]
        s = self._signal[1.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 1)
        self.assertEqual(s[5], 2)
        self.assertEqual(s[6], 2)

    def testBackToZero(self):
        '''Test we can return to zero.'''
        self._playEventsTo(4.0)
        s = self._signal[2.0]
        self.assertEqual(s[1], -1)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 0)
        self.assertEqual(s[4], 1)
        self.assertEqual(s[5], 2)
        self.assertEqual(s[6], 2)
        s = self._signal[4.0]
        self.assertEqual(s[1], -2)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], -1)
        self.assertEqual(s[4], 0)
        self.assertEqual(s[5], 1)
        self.assertEqual(s[6], 1)
        s = self._signal[0.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 1)
        self.assertEqual(s[4], 2)
        self.assertEqual(s[5], 3)
        self.assertEqual(s[6], 3)

    def testSetUp(self):
        '''Test we transfer initial compartments properly.'''
        self._p = SIR()
        self._e = StochasticDynamics(self._p, FixedNetwork(self._g))
        self._p.setNetwork(self._g)
        self._p.reset()
        self._p.build(self._params)
        self._p.setUp(self._params)
        self._p.setCompartment(1, SIR.INFECTED)
        signal = Signal()
        gen = SIRProgressSignalGenerator(signal)
        gen.setExperiment(self._e)
        gen.setUp(self._g, self._params)
        s = signal[0.0]
        self.assertEqual(s[1], 0)
        self.assertEqual(s[2], 1)
        self.assertEqual(s[3], 1)
        self.assertEqual(s[4], 2)
        self.assertEqual(s[5], 3)
        self.assertEqual(s[6], 3)

    # TODO Test that disconnnected sub-graphs generate +/- infinity


    # ---------- Soak tests ----------

    def testInvariantsAdd(self):
        '''Test invariants when adding infected nodes.'''
        g = Graph()
        g.add_nodes_from([1, 2, 3, 4, 5, 6])
        g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])
        evs = [(1.0, SIR.INFECTED, (3, 1)),
               (2.0, SIR.INFECTED, (4, 3)),
               (3.0, SIR.INFECTED, (5, 4)),
               (4.0, SIR.INFECTED, (2, 1)),
               (5.0, SIR.INFECTED, (6, 4))]
        self.checkInvariants(g, [1], evs)

    def testInvariantsRemove(self):
        '''Test invariants when adding and then removing infected nodes.'''
        g = Graph()
        g.add_nodes_from([1, 2, 3, 4, 5, 6])
        g.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])
        evs = [(6.0, SIR.REMOVED, 4),
               (7.0, SIR.REMOVED, 5),
               (8.0, SIR.REMOVED, 1),
               (9.0, SIR.REMOVED, 3),
               (10.0, SIR.REMOVED, 2),
               (11.0, SIR.REMOVED, 6)]
        self.checkInvariants(g, [1, 2, 3, 4, 5, 6], evs)

    def testInvariantsLatticeManual(self):
        '''Test invariants on a larger network.'''
        g = convert_node_labels_to_integers(grid_graph(dim=(6, 6)), first_label=1)
        evs = [(1.0, SIR.INFECTED, (22, 21))]
        self.checkInvariants(g, [21], evs)

    def checkInvariants(self, g, initialInfecteds, evs, endState=False):
        p = SIR()
        x = StochasticDynamics(p, FixedNetwork(g))
        x.setNetwork(g)
        params = dict({SIR.P_INFECTED: 0.0,
                       SIR.P_INFECT: 0.0,
                       SIR.P_REMOVE: 0.0})
        signal = Signal()
        generator = SIRProgressSignalGenerator(signal)
        generator.setExperiment(x)
        generator.setProcess(p)
        p.reset()
        p.build(params)
        p.setUp(params)
        for i in initialInfecteds:
            p.changeCompartment(i, SIR.INFECTED)
        generator.setUp(g, params)

        ns = list(g.nodes())
        susceptibles = set()
        infecteds = set()
        removeds = set()
        for n in ns:
            c = p.getCompartment(n)
            if c == SIR.SUSCEPTIBLE:
                susceptibles.add(n)
            elif c == SIR.INFECTED:
                infecteds.add(n)
            elif c == SIR.REMOVED:
                removeds.add(n)
            else:
                self.assertTrue(False, 'Invalid compartment {c}')

        for (t, etype, e) in evs:
            if etype == SIR.INFECTED:
                p.infect(t, e)
            elif etype == SIR.REMOVED:
                p.remove(t, e)
            generator.event(t, etype, e)
            sig = signal[t]

            # advance the sets
            self.assertIn(etype, [SIR.INFECTED, SIR.REMOVED])
            if etype == SIR.INFECTED:
                (s, _) = e
                susceptibles.remove(s)
                infecteds.add(s)
            elif etype == SIR.REMOVED:
                infecteds.remove(e)
                removeds.add(e)

            # test all nodes have an entry in the signal
            self.assertCountEqual(ns, sig.keys())

            # test that all infecteds are zeros
            for n in infecteds:
                #print(f'inf check {n}')
                self.assertEqual(sig[n], 0)

            # check all susceptibles are the right distance from the boundary
            self.checkSusceptibles(g, sig, susceptibles, infecteds, removeds)

            # check all removeds are the right distance from the boundary
            self.checkRemoveds(g, sig, susceptibles, infecteds, removeds)

        if endState:
            #  check the end state
            self.assertEqual(len(susceptibles), 0)
            self.assertEqual(len(susceptibles) + len(removeds), g.order())

    def checkSusceptibles(self, g, sig, susceptibles, infecteds, removeds):
        ss = susceptibles.copy()
        while len(ss) > 0:
            n = ss.pop()
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

            # check our distance to the infected boundary is correct
            self.assertEqual(d,
                             self.shortestPath(g, n, infecteds, susceptibles))

    def checkRemoveds(self, g, sig, susceptibles, infecteds, removeds):
        rr = removeds.copy()
        onpath = set(susceptibles).copy().union(set(removeds))
        while len(rr) > 0:
            n = rr.pop()
            #print(f'sus check {n}')
            d = sig[n]
            #print(f'd = {d}')
            # all neighbours should have distances differing by at most one
            # from us (if they're removeds), or be infecteds (in which case
            # our distance should be 1), or be susceptibles
            for m in g.neighbors(n):
                if m in removeds:
                    #print(n, m, d, sig[m])
                    self.assertTrue(abs(sig[m] - d) <= 1)
                elif m in infecteds:
                    self.assertEqual(d, -1)

            # check our distance to the infected boundary is correct
            dprime = self.shortestPath(g, n, infecteds, onpath)
            if dprime is None:
                # if we can't find a shortest path then there should be no infecteds left
                self.assertEqual(len(infecteds), 0)
            else:
                # signal is 0 - shortest path
                self.assertEqual(d, -dprime)

    def shortestPath(self, g, n, targets, onpath):
        distance = [(0, n)]
        visited = set()
        while len(distance) > 0:
            (d, n) = heappop(distance)
            visited.add(n)
            dprime = d + 1
            ms = g.neighbors(n)
            for m in ms:
                if m not in visited:
                    if m in targets:
                        # found a node in the target set
                        return dprime
                    else:
                        if m in onpath:
                            heappush(distance, (dprime, m))
                        else:
                            visited.add(m)
        return None


if __name__ == '__main__':
    unittest.main()
