# Extract an infection boundary signal from an SIR model
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

from heapq import heappush, heappop
from typing import Dict, Any, List, Tuple, cast
from networkx import Graph, single_source_shortest_path
from epydemic import Node, Edge, SIR, Process, CompartmentedModel
from epydemic_signals import Signal, SignalGenerator


class InfectionBoundarySignalGenerator(SignalGenerator):
    '''Create the infection boundary signal of an SIR
    epidemic. This signal is defined on nodes as the number of SI edges
    incident on that node, and so is 0 by definition on any other
    than S or I edges.

    :param s: the signal
    '''

    def __init__(self, p: Process, s: Signal):
        super().__init__(p, s)

        # register the event handlers
        self.addEventTypeHandler(SIR.INFECTED, self.infect)
        self.addEventTypeHandler(SIR.REMOVED, self.remove)

    def setUp(self, g: Graph):
        '''Capture the initial signal.'''
        super().setUp(g)

        signal = self.signal()
        g = self.network()
        p = self.process()
        cm = cast(CompartmentedModel, p)

        # initialise signal t 0 everywhere, to make sure
        # we have an entry for all nodes
        s = signal[0.0]
        for n in g.nodes():
            s[n] = 0

        # traverse all the infected nodes, counting incident edges
        for n in g.nodes():
            if cm.getCompartment(n) == SIR.INFECTED:
                # count the incident SI edges
                si = 0
                for (_, m) in g.edges(n):
                    if cm.getCompartment(m) == SIR.SUSCEPTIBLE:
                        si += 1
                        if m in s:
                            s[m] += 1
                        else:
                            s[m] = 1
                s[n] = si

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
        si = 0
        for m in g.neighbors(n):
            c = cm.getCompartment(m)
            if c == SIR.SUSCEPTIBLE:
                # new SI edge, incrment both us and the neighbour
                si += 1
                s[m] += 1
            elif c == SIR.INFECTED:
                # former SI edge, decrement the other end
                s[m] -= 1
        s[n] = si

    def remove(self, t: float, n: Node):
        '''Update signal on removal. This decrements any susceptibble
        neighbours' signal values and sets us to 0.


        :param t: the simulation time
        :param e: the SI edge'''
        signal = self.signal()
        s = signal[t]
        g = self.network()
        cm = cast(CompartmentedModel, self.process())

        for m in g.neighbors(n):
            c = cm.getCompartment(m)
            if c == SIR.SUSCEPTIBLE:
                # former SI edge, decrement the neighbour
                s[m] -= 1
        s[n] = 0
