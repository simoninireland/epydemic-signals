# Extract an infection boundary signal from an SIR model
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

from heapq import heappush, heappop
from typing import Dict, Any, List, Tuple, cast
from networkx import Graph, single_source_shortest_path
from epydemic import Node, Edge, SIR, Process, CompartmentedModel
from epydemic_signals import Signal, SignalGenerator


class InfectionBoundarySignalGenerator(SignalGenerator):
    '''Create the infection boundary signal of an SIR
    epidemic. This signal is defined as the number of
    incident SI edges on infected nodes, and zero everywhere
    else. It therefore represents the "infection potential"
    of nodes: those with higher values are able to potentially
    infect mode nodes.

    :param s: the signal
    '''

    def __init__(self, s: Signal = None):
        super().__init__(s)

        # register the event handlers
        self.addEventTypeHandler(SIR.INFECTED, self.infect)
        self.addEventTypeHandler(SIR.REMOVED, self.remove)

    def process(self) -> Process:
        '''Return the process this generator is monitoring.

        :returns: the process'''
        return self.experiment().process()

    def setUp(self, g: Graph, params: Dict[str, Any]):
        '''Capture the initial signal.

        :param g: the network
        :param params: the experimental parameters'''
        super().setUp(g, params)

        signal = self.signal()
        g = self.network()
        p = self.process()
        cm = cast(CompartmentedModel, p)

        # initialise signal to 0 everywhere, to make sure
        # we have an entry for all nodes
        s = signal[0.0]
        for n in g.nodes():
            s[n] = 0

        # traverse all the infected nodes, counting incident edges
        for n in g.nodes():
            if cm.getCompartment(n) == SIR.INFECTED:
                # count the incident SI edges
                for m in g.neighbors(n):
                    if cm.getCompartment(m) == SIR.SUSCEPTIBLE:
                        s[n] += 1

    def infect(self, t: float, e: Edge):
        '''Change the signal on infection. This involves removing any
        edges to neighbouring infecteds, and adding those to neighbouring
        susceptibles.

        :param t: the simulation time
        :param e: the SI edge'''
        signal = self.signal()
        s = signal[t]
        g = self.network()
        cm = cast(CompartmentedModel, self.process())
        (n, _) = e

        # traverse all neighbours and count SI edges
        for m in g.neighbors(n):
            c = cm.getCompartment(m)
            if c == SIR.SUSCEPTIBLE:
                # new SI edge, increment us
                s[n] += 1
            elif c == SIR.INFECTED:
                # former SI edge, decrement the other end
                s[m] -= 1

    def remove(self, t: float, n: Node):
        '''Update signal on removal. This sets our value to 0.

        :param t: the simulation time
        :param n: the node'''
        signal = self.signal()
        s = signal[t]
        s[n] = 0
