# An SIR model that starts with a single infection point
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

import numpy
from typing import Dict, Any
from epydemic_signals import HealingSIR


class OneInfectionSIR(HealingSIR):
    '''A :class:`HealingSIR` model with a single randomly-chosen point of
    infection.'''

    def build(self, params: Dict[str, Any]):
        '''Make sure there's no random infection.

        :param params: the experimental parameters'''
        params[self.P_INFECTED] = 0.0
        super().build(params)

    def initialCompartments(self):
        '''Select a single node to infect.'''
        g = self.network()
        ns = set(g.nodes())
        N = len(ns)

        # choose one node and infect it
        rng = numpy.random.default_rng()
        n = rng.integers(N)
        self.changeInitialCompartment(n, self.INFECTED)
        self.markHit(n, 0.0)

        # mark all other nodes as susceptible
        ns.remove(n)
        for n in ns:
            self.changeInitialCompartment(n, self.SUSCEPTIBLE)
