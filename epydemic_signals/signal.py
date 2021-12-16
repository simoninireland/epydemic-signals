# Signals
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

from typing import Dict, TypeVar, Generic, List, Iterable
from networkx import Graph
from epydemic import Node, Process
from epydemic_signals import TimedDict


# Type variable for signal values
V = TypeVar('V')


class Signal(Generic[V]):
    '''Encode a signal on a network.

    A signal -- or strictly speaking a node signal -- associates a mapping
    from nodes to values for every point in time.'''

    def __init__(self):
        self._network: Graph = None                    # the domain of the signal
        self._dict: TimedDict[float, V] = TimedDict()  # the signal data structure


    # ---------- Accessing the signal ----------

    def setNetwork(self, g: Graph):
        '''Set the network over which this signal is defined.

        :param g: the network'''
        self._network = g

    def network(self) -> Graph:
        '''Return the network over which this signal is defined.

        :returns: the network'''
        return self._network

    def transitions(self) -> List[float]:
        '''Return a list of times at which the signal changes, in
        ascending order.

        :returns: a list of times'''
        return self._dict.updates()

    def values(self) -> Iterable[V]:
        '''Return all the values the signal takes, at any point and time.

        :returrns: the values'''
        return self._dict.valuesAtSomeTime()

    def __getitem__(self, t: float) -> Dict[Node, V]:
        '''Extract the mapping of the signal at the given time.

        :param t: the time
        :returns: a dict from nodes to value at the given times'''
        return self._dict[t]

    def __len__(self) -> int:
        '''Return the number of transition points in the signal, the
        times when it changed.

        :returns: the number of transitions'''
        return len(self._dict)
