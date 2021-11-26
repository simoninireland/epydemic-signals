# Base class for signal generators
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

from typing import Callable
from networkx import Graph
from epydemic import Element, Node, Edge
from epydemic_signals import Signal


# Types for handling events
EventHandler =  Callable[[float, Element], None]   #: Type of event-handler function.


class SignalGenerator:
    '''The base class for signal generators. A signal generator is
    attached to a signal dynamics to repond to the events generated
    by a simulation, generating the corresponding signal.

    :param s: the signal being generated'''

    def __init__(self, s: Signal):
        self._signal = s
        self._network = s.network()
        self._typeHandler: Dict[str, EventHandler] = dict()

    def signal(self) -> Signal:
        '''The signal being generated.

        :returns: the signal'''
        return self._signal

    def network(self) -> Graph:
        '''Return the network the signal is being genearted for.

        :returns: the network'''
        return self._network

    def addEventTypeHandler(self, etype: str, eh: EventHandler):
        '''Register a handler for the given event. The handler is called
        whenever the event is fired. Several handlers can be registered
        to a given event type, in which case they are called in the order
        they are added.

        :param etype: the event type
        :param eh: the event handler'''
        if etype in self._typeHandler:
            # append the event to the list of handlers
            self._typeHandler[etype].append(eh)
        else:
            # add the handler
            self._typeHandler[etype] = [eh]

    def event(self, t: float, etype: str, e: Element):
        '''Respond to the given event. Events are dispatched to the
        handlers registered for them: any events for which there is no
        handler are silently ignored.

        :param t: the simulation time
        :param etype: the event type
        :param e: the element'''
        print(t, etype, e)
        ehs = self._typeHandler.get(etype, [])
        for eh in ehs:
            eh(t, e)
