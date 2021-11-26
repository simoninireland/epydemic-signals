# A mixin to tap the event stream of a dynamics to a signal generator
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

from epydemic import Element
from epydemic_signals import SignalGenerator

class SignalDynamics:
    '''A mixin class used to add a "tap" for an event stream to a
    :class:`NetworkDynamics` class.'''

    def __init__(self):
        super().__init__()
        self._generator = None

    def setSignalGenerator(self, gen: SignalGenerator):
        '''Set the signal generator that will be passed events.

        :param gen: the signal generator'''
        self._generator = gen

    def eventFired(self, t: float, etype: str, e: Element):
        '''Pass all fired events to the signal generator.

        :param t: the simulation time
        :param etype: the event type
        :param e: the element'''
        self._generator.event(t, etype, e)
