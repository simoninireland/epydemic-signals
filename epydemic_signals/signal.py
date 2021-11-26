# Signals
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

from typing import Dict, Any, List
from networkx import Graph
from epydemic import Node
from epydemic_signals import TimedDict


class Signal:
    '''Encode a signal on a network.

    A signal -- or strictly speaking a node signal -- associates a mapping
    from nodes to values for every point in time.

    :param g: the network'''

    def __init__(self, g: Graph):
        self._network = g           # the network
        self._dict = TimedDict()    # the mapping


    # ---------- Accessing the signal ----------

    def network(self) -> Graph:
        '''Return the network over which this signal is defined.

        :returns: the network'''
        return self._network

    def transitions(self) -> List[float]:
        '''Return a list of times at which the signal changes, in
        ascending order.

        :returns: a list of times'''
        return self._dict.updates()

    def getBounds(self) -> List[float]:
        vs = list(self._dict.valuesAtSomeTime())
        if len(vs) > 2:
            vs.sort()
            return (vs[0], vs[-1])
        elif len(vs) == 1:
            return (vs[0], vs[0])
        else:
            raise ValueError('Empty signal')

    def __getitem__(self, t: float) -> Dict[Node, float]:
        '''Extract the mapping of the signal at the given time.

        :param t: the time
        :returns: a dict from nodes to values'''
        return self._dict[t]
