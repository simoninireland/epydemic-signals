# Tests taps in synchronoussignals
#
# Copyright (C) 2021-2022 Simon Dobson
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
from test.progresssignalinvariants import SIRProgressSignalInvariants
from epydemic import SIR, StochasticDynamics, ProcessSequence, FixedNetwork, ERNetwork
from epyc import Experiment, Lab, ParallelLab
from networkx import Graph, fast_gnp_random_graph, grid_graph, convert_node_labels_to_integers

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
        e = SynchronousSignalDynamics(sir, FixedNetwork(g))

        sig = Signal()
        gen1 = SIRProgressSignalGenerator(sig)
        e.attachSignalGenerator(gen1, sir)
        gen2 = SIRProgressSignalInvariants(sig, gen1)    # checks the same signal
        e.attachSignalGenerator(gen2, sir)

        params = dict()
        params[SIR.P_INFECT] = pInfect
        params[SIR.P_REMOVE] = pRemove

        rc = e.set(params).run(fatal=True)

        self.assertTrue(len(sig) > 0)
        self.assertIsNotNone(sig.network())
        self.assertIsNotNone(gen1.network())
        self.assertIsNotNone(gen2.network())
