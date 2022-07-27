# Base class for signal generators
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

from typing import Callable, Dict, Any
from networkx import Graph
from epydemic import Element, Node, Edge, Process, NetworkExperiment
from epydemic_signals import Signal


# Types for handling events
EventHandler = Callable[[float, Element], None]   #: Type of event-handler function.


class SignalGenerator:
    '''The base class for signal generators. A signal generator is
    attached to a signal dynamics to repond to the events generated
    by a simulation, generating the corresponding signal.

    By default a basic instance of :class:`Signal` is used for the signal.
    Specific instances and/or sub-classes can be provided if needed. It's
    often useful to set a default signal with a meaningful name in the
    constructor, so that different signal objects can be differentiated if
    written into a notebook.

    :param s: (optional) the signal being generated (creates one if missing)'''

    def __init__(self, s: Signal = None):
        if s is None:
            s = Signal()
        self._experiment = None
        self._process = None
        self._signal = s
        self._typeHandler: Dict[str, EventHandler] = dict()

    def setSignal(self, s: Signal):
        '''Set the signal being generated. This allows re-use of a signal generator
        to create more signals, or to change the name (type) of signal generated.

        :param s: the signal'''
        self._signal = s

    def signal(self) -> Signal:
        '''The signal being generated.

        :returns: the signal'''
        return self._signal

    def setExperiment(self, e: NetworkExperiment):
        '''Set the experiment the generator is attached to. This is called by
        :meth:`SignalExperiment.addSignalGenerator`.

        :param e: the experiment'''
        self._experiment = e

    def experiment(self) -> NetworkExperiment:
        '''Return the experiment the generator is attached to.

        :returns: the experiment'''
        return self._experiment

    def setProcess(self, p: Process):
        '''Set the process instance the generator is receiving events for.

        :param p: the process'''
        self._process = p

    def process(self) -> Process:
        '''Return the process instance that the generator is receiving events for.

        :returns: the process'''
        return self._process

    def network(self) -> Graph:
        '''Return the network the signal is being generated over. This
        is actually defined in the signal rather than the generator.

        :returns: the network'''
        return self._signal.network()


    # ---------- Event type registration ----------

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


    # ---------- Tap methods ----------

    def setUp(self, g: Graph, params: Dict[str, Any]):
        '''Notify the signal generator that the simulation has started.
        This can be overridden by sub-classes to get the generator ready
        to record the signal. The default sets the network for the signal,
        and so needs to be called by any overriding methods.

        :param g: the network
        :param params: the experimental parameters'''
        self._signal.setNetwork(g)

    def tearDown(self):
        '''Notify the signal generator that the simulation has ended.
        This can be overridden by sub-classes to destroy any working
        data structures, close or write to files. The default does nothing.
        '''
        pass

    def event(self, t: float, etype: str, e: Element):
        '''Respond to the given event. Events are dispatched to the
        handlers registered for them: any events for which there is no
        handler are silently ignored.

        Sub-classes can override this method to avoid the dispatch mechanism
        and simply handle all events in one place.

        :param t: the simulation time
        :param etype: the event type
        :param e: the element'''
        ehs = self._typeHandler.get(etype, [])
        for eh in ehs:
            eh(t, e)
