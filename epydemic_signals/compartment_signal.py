# Extract the evolution of compartments as a signel
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

from typing import cast
from epydemic import Node, Edge, Element, Process, CompartmentedModel
from epydemic_signals import Signal, SignalGenerator


class CompartmentSignalGenerator(SignalGenerator):
    '''Create a signal from the way compartments change. Works for any
    compartmented model.

    :param s: the signal
    '''

    def __init__(self, s: Signal):
        super().__init__(s)

    def captureCompartments(self, t):
        '''Capture the state of the network in terms of compartments. The
        signal for a node at time t is simply its compartment at that time.

        :param t: the simulation time'''
        g = self.network()
        sig = self.signal()[t]
        cm = cast(CompartmentedModel, p)
        for n in g.nodes():
            sig[n] = cm.getCompartment(n)

    def setUp(self):
        '''Capture the initial state of the network.

        :param g: the network
        :param p: the process'''
        self.captureCompartments(0.0)

    def event(self, t: float, etype: str, e: Element):
        '''Respond to events by snapshotting the state of the compartments.

        :param t: the simulation time
        :param etype: the event type (not used)
        :parram e: the element (not used)'''
        self.captureCompartments(t)
