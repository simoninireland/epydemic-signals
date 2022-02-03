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

from uuid import uuid4
from typing import Dict, TypeVar, Generic, List, Tuple, Iterable
from networkx import Graph
from epydemic import Node, Process
from epydemic_signals import TimedDict


# Type variable for signal values
V = TypeVar('V')


class Signal(Generic[V]):
    '''Encode a time-varying signal on a network.

    A signal -- or strictly speaking a node signal -- associates a mapping
    from nodes to values for every point in time.

    Signals are recorded as experimental results, using the name of
    the signal as part of the key for the representation in the lab notebook.
    If no name is given then a UUID is generated, which is (pretty) unique but
    entirely uninformative and so to be avoided.

    :param g: (optional) the network over which the signal is defined
    :param name: (optional) the name of the signal'''

    def __init__(self, g: Graph = None, name: str = None):
        # fill in the defaults
        if name is None:
            name = uuid4()
        self._name = name                              # the signal name
        self._network: Graph = g                       # the domain of the signal
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

    def name(self ) -> str:
        '''Return the signal name.

        :returns: the signal name'''
        return self._name

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


    # ---------- Converting ----------

    def toTimeSeries(self) -> Dict[Node, List[V]]:
        '''Convert a node signal to a collection of time series for each
        node. The time series are all sampled at the same points, corresponding
        to the series returned by :meth:`transitions`.

        :returns: a dict of time series'''
        g = self.network()
        ns = list(g.nodes())
        N = g.order()
        tss = dict()
        for n in ns:
            tss[n] = []

        # we run through the times, retrieving the value for each node
        # this is more efficient that traversing per-node due to the way
        # TimedDict is implemented
        for t in self.transitions():
            s_t = self[t]
            for n in ns:
                tss[n].append(s_t[n])

        return tss

    def toUpdates(self) -> Tuple[List[float], List[Node], List[V]]:
        '''Convert a node signal to three lists encoding the updates
        made to the value at each node.

        At present this doesn't handle addition of deletion of nodes.

        :returns: a triple of update times, nodes, and values'''
        times = []
        nodes = []
        values = []

        ns = list(self.network().nodes())
        ts = self.transitions()

        # initial values
        s_t = self[ts[0]]
        for n in ns:
            times.append(ts[0])
            nodes.append(n)
            values.append(s_t[n])

        # updates
        s_t_old = s_t
        for t in ts[1:]:
            s_t = self[t]
            for n in ns:
                if s_t[n] != s_t_old[n]:
                    # value has changed, record as an update
                    times.append(t)
                    nodes.append(n)
                    values.append(s_t[n])
            s_t_old = s_t

        return (times, nodes, values)
