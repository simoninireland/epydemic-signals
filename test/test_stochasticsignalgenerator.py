# Tests taps in stochastic signals
#
# Copyright (C) 2021 Simon Dobson
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
from epydemic import SIR, StochasticDynamics, ProcessSequence, FixedNetwork, ERNetwork
from epyc import Experiment, Lab, ParallelLab
from networkx import Graph, fast_gnp_random_graph, grid_graph, convert_node_labels_to_integers


class SIRProgressSignalInvariants(SignalGenerator):

    def __init__(self, s, ps):
        super().__init__(s)
        self._progressSignalGenerator = ps
        self.addEventTypeHandler(SIR.INFECTED, self.infect)
        self.addEventTypeHandler(SIR.REMOVED, self.remove)

    def setUp(self, g, params):
        super().setUp(g, params)
        g = self.network()
        p = self.experiment().process()

        self._ns = list(g.nodes()).copy()
        self._inf = g.order() + 1
        self._compartment = dict()
        self._compartment[SIR.SUSCEPTIBLE] = set()
        self._compartment[SIR.INFECTED] = set()
        self._compartment[SIR.REMOVED] = set()

        # extract the initial states
        for n in g.nodes():
            c = p.getCompartment(n)
            self._compartment[c].add(n)
        self.checkInvariants(0.0)

    def infect(self, t, e):
        (s, _) = e
        self._compartment[SIR.SUSCEPTIBLE].remove(s)
        self._compartment[SIR.INFECTED].add(s)

        p = self.experiment().process()
        # print('infect', s)
        # for n in self.network().nodes():
        #     print(n, p.getCompartment(n), self.signal()[t][n])
        self.checkInvariants(t)

    def remove(self, t, s):
        self._compartment[SIR.INFECTED].remove(s)
        self._compartment[SIR.REMOVED].add(s)

        p = self.experiment().process()
        # print('remove', s)
        # for n in self.network().nodes():
        #     print(n, p.getCompartment(n), self.signal()[t][n])
        self.checkInvariants(t)

    def checkInvariants(self, t):
        #print(f'check at {t}')
        g = self.network()
        signal = self.signal()
        sig = signal[t]

        # test all nodes have an entry in the signal
        for n in self._ns:
            if n not in sig:
                raise Exception(f'No key {n} in signal')

        # test that all infecteds are zeros
        for n in self._compartment[SIR.INFECTED]:
            s = sig[n]
            if s != 0:
                raise Exception(f'Infected node {n} signal should be 0 but is {s}')

        # check all susceptibles are the right distance from the boundary
        self.checkSusceptibles(g, sig)

        # check all removeds are the right distance from the boundary
        self.checkRemoveds(g, sig)

        # white-box testing of the algorithm
        self.checkBoundaries(t)

    def checkBoundaries(self, t):
        signal = self.signal()
        sig = signal[t]
        gen = self._progressSignalGenerator

        # check all susceptibles and removeds have a boundary
        for n in self._compartment[SIR.SUSCEPTIBLE]:
            if sig[n] == gen._inf:
                continue
            if n not in gen._boundary:
                raise Exception(f'No boundary for susceptible {n}')
        for n in self._compartment[SIR.REMOVED]:
            if sig[n] == -gen._inf:
                continue
            if n not in gen._boundary:
                raise Exception(f'No boundary for removed {n}')

        # check all infecteds have coboundaries
        for n in self._compartment[SIR.INFECTED]:
            if n not in gen._coboundary_S:
                raise Exception(f'No S coboundary for infected {n}')
            if n not in gen._coboundary_R:
                raise Exception(f'No R coboundary for infected {n}')

        # check all boundary nodes lie in the appropriate coboundary
        for n in self._compartment[SIR.SUSCEPTIBLE]:
            if sig[n] == gen._inf:
                continue
            if gen._boundary[n] not in gen._coboundary_S:
                raise Exception(f'No S coboundary for boundary of susceptible {n}', gen._boundary[n] )
            if n not in gen._coboundary_S[gen._boundary[n]]:
                raise Exception(f'S coboundary mismatch for susceptible {n}')
        for n in self._compartment[SIR.REMOVED]:
            if sig[n] == -gen._inf:
                continue
            if gen._boundary[n] not in gen._coboundary_R:
                raise Exception(f'No R coboundary for boundary of removed {n}', gen._boundary[n])
            if n not in gen._coboundary_R[gen._boundary[n]]:
                raise Exception(f'R coboundary mismatch for removed {n}')

    def checkSusceptibles(self, g, sig):
        onpath = set(self._compartment[SIR.SUSCEPTIBLE]).copy()
        for n in self._compartment[SIR.SUSCEPTIBLE]:
            #print(f'sus check {n}')
            d = sig[n]

            # susceptible signals should be > 0
            if d <= 0:
                raise Exception(f'Susceptible signal invalid {d}')

            # all neighbours should have distances differing by at most one
            # from us (if they're susceptibles), or be infecteds (in which case
            # our distance should be 1), or be removeds
            for m in g.neighbors(n):
                if m in self._compartment[SIR.SUSCEPTIBLE]:
                    #print(n, m, d, sig[m])
                    if not (abs(sig[m] - d) <= 1):
                        raise Exception(f'Susceptible {n} neighbour {m} signal diff too large', d, sig[m])
                elif m in self._compartment[SIR.INFECTED]:
                    if d != 1:
                        raise Exception(f'Susceptible {m} signal next to infected should be 1 but is {d}')

            # check our distance to the infected boundary is correct
            dprime = self.shortestPath(g, n, self._compartment[SIR.INFECTED], onpath)
            if dprime is not None:
                if d != dprime:
                    raise Exception(f'Susceptible {m} path should be {dprime} but is {d}')

    def checkRemoveds(self, g, sig):
        onpath = set(self._compartment[SIR.SUSCEPTIBLE]).copy().union(set(self._compartment[SIR.REMOVED]))
        for n in self._compartment[SIR.REMOVED]:
            #print(f'rem check {n}')
            d = sig[n]

            # removed signals should be < 0
            if d >= 0:
                raise Exception(f'Removed signal invalid {d}')

            #print(f'd = {d}')
            # all neighbours should have distances differing by at most one
            # from us (if they're removeds), or be infecteds (in which case
            # our distance should be 1), or be susceptibles
            for m in g.neighbors(n):
                if m in self._compartment[SIR.REMOVED]:
                    #print(n, m, d, sig[m])
                    if not (abs(sig[m] - d) <= 1):
                        raise Exception(f'Removed {n} neighbour {m} signal diff too large', d, sig[m])
                elif m in self._compartment[SIR.INFECTED]:
                    if d != -1:
                        raise Exception(f'Removed {n} signal should be -1 but is {d}')

            # check our distance to the infected boundary is correct
            dprime = self.shortestPath(g, n, self._compartment[SIR.INFECTED], onpath)
            if dprime is not None:
                if d != -dprime:
                    raise Exception(f'Removed {n} signal should be -{dprime} but is {d}')

    def shortestPath(self, g, s, targets, onpath):
        distance =[]
        heappush(distance, (0, s))
        seen = set([s])
        while len(distance) > 0:
            (d, n) = heappop(distance)
            if n in targets:
                # found a node in the target set
                return d

            # if we're potentially on the path, add all neighbours to be visited
            if n in onpath:
                dprime = d + 1
                ms = g.neighbors(n)
                for m in ms:
                    if m not in seen:
                        seen.add(m)
                        heappush(distance, (dprime, m))
        return None


class StochasticSignalDynamicsTests(unittest.TestCase):

    def testInvariantsLattice(self):
        '''Test invariants for an epidemic on a lattice.'''
        g = convert_node_labels_to_integers(grid_graph(dim=(10, 10)), first_label=1)
        pInfect = 0.8
        pRemove = 0.1

        params = dict()
        params[SIR.P_INFECT] = pInfect
        params[SIR.P_REMOVE] = pRemove

        sir = OneInfectionSIR()
        e = StochasticSignalDynamics(sir, FixedNetwork(g))

        sig = Signal()
        gen1 = SIRProgressSignalGenerator(sig)
        e.addSignalGenerator(gen1)
        gen2 = SIRProgressSignalInvariants(sig, gen1)    # checks the same signal
        e.addSignalGenerator(gen2)

        #rc = e.set(params).run(fatal=True)
        lab = Lab()
        lab[SIR.P_INFECT] = pInfect
        lab[SIR.P_REMOVE] = pRemove
        lab.runExperiment(e)

        self.assertTrue(len(sig) > 0)
        self.assertIsNotNone(sig.network())
        self.assertIsNotNone(gen1.network())
        self.assertIsNotNone(gen2.network())

    def testCreatWith(self):
        '''Test using createWith().'''

        def create(lab):
            g = convert_node_labels_to_integers(grid_graph(dim=(10, 10)), first_label=1)
            pInfect = 0.8
            pRemove = 0.1

            params = dict()
            params[SIR.P_INFECT] = pInfect
            params[SIR.P_REMOVE] = pRemove

            sir = OneInfectionSIR()
            e = StochasticSignalDynamics(sir, FixedNetwork(g))

            global sig, gen1, gen2
            sig = Signal()
            gen1 = SIRProgressSignalGenerator(sir, sig)
            e.addSignalGenerator(gen1)
            gen2 = SIRProgressSignalInvariants(sir, sig, gen1)    # checks the same signal
            e.addSignalGenerator(gen2)

            lab[SIR.P_INFECT] = pInfect
            lab[SIR.P_REMOVE] = pRemove
            lab.runExperiment(e)

        lab = Lab()
        lab.createWith('test', create)

        self.assertTrue(len(sig) > 0)
        self.assertIsNotNone(sig.network())
        self.assertIsNotNone(gen1.network())
        self.assertIsNotNone(gen2.network())

    def testInvariantsLarge(self):
        '''Test invariants as we run a larger epidemic.'''
        N = 1000
        kmean = 10
        g = fast_gnp_random_graph(N, kmean / N)
        pInfect = 0.01
        pRemove = 0.001

        params = dict()
        params[SIR.P_INFECT] = pInfect
        params[SIR.P_REMOVE] = pRemove

        signal = Signal()
        sir = OneInfectionSIR()
        e = StochasticSignalDynamics(sir, FixedNetwork(g))

        sig = Signal()
        gen1 = SIRProgressSignalGenerator(sir, sig)
        e.addSignalGenerator(gen1)
        gen2 = SIRProgressSignalInvariants(sir, sig, gen1)    # checks the same signal
        e.addSignalGenerator(gen2)

        rc = e.set(params).run(fatal=True)

        self.assertTrue(len(sig) > 0)
        self.assertIsNotNone(sig.network())
        self.assertIsNotNone(gen1.network())
        self.assertIsNotNone(gen2.network())

    def testInvariantER(self):
        '''Test invariants over an ER network, for less structure.'''
        lab = Lab()
        lab[ERNetwork.N] = 5000
        lab[ERNetwork.KMEAN] = 20
        lab[SIR.P_INFECTED] = 0.001
        lab[SIR.P_REMOVE] = 0.002
        lab[SIR.P_INFECT] = 0.00015

        sir = SIR()
        e = StochasticSignalDynamics(sir, ERNetwork())
        progress = Signal()
        gen = SIRProgressSignalGenerator(progress)
        e.addSignalGenerator(gen)
        gen2 = SIRProgressSignalInvariants(progress, gen)    # checks the same signal
        e.addSignalGenerator(gen2)

        lab.runExperiment(e)
        rc = lab.results()[0]
        if not rc[Experiment.METADATA][Experiment.STATUS]:
            print(rc[Experiment.METADATA][Experiment.TRACEBACK])
            raise Exception('Failed')


if __name__ == '__main__':
    unittest.main()
