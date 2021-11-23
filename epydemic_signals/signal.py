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


class Signal:
    '''Encode a signal on a network. A signal is simply a value associated with
    a time and a node. (Strictly speaking these should be called "node signals".)

    :param g: the network'''

    def __init__(self, g: Graph):
        self._network = g      # the network

        # signal elements
        self._base = dict()    # the signal at t = 0
        self._diffs = dict()   # the diffs
        self._transitions = [] # list of transition times
        self._dirty = False    # True of the diffs need sorting
        self._bounds = None    # max, min of signal

        # signal access
        self._signal = dict()  # the signal at the current time
        self._index = None     # the current time index


    # ---------- Encoding the signal ----------

    def setBaseSignal(self, f: Dict[Node, float]):
        '''Set the signal at time t = 0.

        :parram f: a dict from node to value'''''
        self._base = f.copy()
        self._diffs = dict()
        self._diffs[0.0] = dict()
        self._transitions = [0.0]
        self._dirty = False

    def addDiff(self, t: float, vals: Dict[Node, float]):
        '''Add a set of differences to the signal.

        :param t: the time
        :param vals: a dict of nodes to changes in value'''
        self._diffs[t] = vals
        self._transitions.append(t)
        self._dirty = True

    def setBounds(self, mn: float, mx: float):
        self._bounds = [mn, mx]


    # ---------- Accessing the signal ----------

    def network(self) -> Graph:
        '''Return the network over which this signal is defined.

        :returns: the network'''
        return self._network

    def transitions(self) -> List[float]:
        '''Return a list of times at which the signal changes, in
        ascending order.

        :returns: a list of times'''
        if self._dirty:
            # sort the diffs by time
            self._transitions.sort()
            self._dirty = False
        return self._transitions

    def getBounds(self) -> List[float]:
        return self._bounds

    def _applyForwards(self, vals: Dict[Node, float]):
        '''Apply diffs moving forwards in time.

        :param vals: a dict of nodes to changes in value'''
        for n in vals.keys():
            self._signal[n] += vals[n]

    def _applyBackwards(self, vals: Dict[Node, float]):
        '''Apply diffs moving backwards in time.

        :param vals: a dict of nodes to changes in value'''
        for n in vals.keys():
            self._signal[n] -= vals[n]

    def getTime(self) -> float:
        '''Return the current signal time.

        :returns: the time'''''
        return self._transitions[self._index]

    def setTime(self, t: float):
        '''Set the current signal time, applying diffs as necessary. The
        time is set to the nearest transition boundary less than or
        equal to the time requested.

        :param t: the new time'''
        if t < 0:
            raise ValueError(f'Can\'t compute signal at negative time {t}')
        if self._index is None:
            # initialise to the base signal
            self._signal = self._base.copy()
            self._index = 0

        # apply the necessary diffs
        i = self._index
        current = self._transitions[i]
        if t < current:
            # before the current time
            while self._transitions[self._index] > t:
                self._applyBackwards(self._diffs[self._transitions[self._index]])
                self._index -= 1
        elif t > current:
            # after current time
            ndiffs = len(self._transitions)
            while self._index < ndiffs - 1 and self._transitions[self._index + 1] <= t:
                self._applyForwards(self._diffs[self._transitions[self._index + 1]])
                self._index += 1

    def __getitem__(self, t: float) -> Dict[Node, float]:
        '''Extract the mapping of the signal at the given time.

        :param t: the time
        :returns: a dict from nodes tto values'''
        self.setTime(t)
        return self._signal
