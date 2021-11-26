# Dicts with changes logged in time
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

from typing import Generic, TypeVar, Dict, Tuple, Any, List, Iterable
from networkx import Graph
from epydemic import Node


# Type variables for dict keys and values
K = TypeVar('K')
V = TypeVar('V')


class TimedDictView(Generic[K, V]):
    '''A view of a timed dict snapped at a particular time.

    :param d: the timed dict's diff structure
    :param t: the time'''

    def __init__(self, d: Dict[K, List[Tuple[float, bool, V]]], t : float):
        self._dict = d                     # diff from key to diff list
        self._time = t                     # projection time
        self._now: Dict[K, int] = dict()   # dict of key to index in diff list of last update in diff list
        self._project()


    # ---------- projection ----------

    def _project(self):
        '''Project-out the values in the dict ast the current time.'''
        self._now = dict()
        for k in self._dict:
            i = self._updateBefore(k)
            if i >= 0:
                (_, b, _) = self._dict[k][i]
                if b:
                    # update was a set, include this key in the projection
                    self._now[k] = i
                else:
                    # update was a delete
                    pass

    def _updateBefore(self, k: K) -> int:
        '''Return the index to the update that occurred on the given entry
        at or immediately preceeding the current time. Returns -1 if there
        is no such entry, meaning that the given k has never been added to
        the dict at any time.

        :param k: the key
        :returns: the index or -1'''
        if k not in self._dict:
            # no entry
            return -1
        else:
            vs = self._dict[k]

            # check times
            (lt, lb, _) = vs[-1]
            if lt <= self._time:
                # last update is before the current time
                return len(vs) - 1
            else:
                # last update is after current time, walk the diff list
                # from the start to find the appropriate update point
                if len(vs) == 1:
                    # only one entry that#'s after the current time,
                    # no point seareching for an appropriate one
                    return -1
                else:
                    for i in range(len(vs) - 1):
                        (ct, cb, _) = vs[i]
                        (nt, nb, _) = vs[i + 1]
                        if nt > self._time:
                            # we've found the last update
                            return i

            # if we get here, there's a problem with the data structures
            raise Exception(f'Corrupted diff list for {k}')


    # ---------- dict interface ----------

    def keys(self) -> Iterable[K]:
        '''Return the keys in the dict at the current time.

        :returns: a list of keys'''
        return self._now.keys()

    def __contains__(self, k: K) -> bool:
        '''Test whether the given k is defined at the current time.

        :param k: the keyed:returns: True if the key is in the dict at the current time'''
        return k in self._keys()

    def values(self) -> Iterable[V]:
        '''Return a list oif values in the dict at the current time.

        :returns: a list of values'''
        # sd: should be lazy?
        vs = set()
        for k in self._now:
            (_, _, v) = self._dict[k][self._now[k]]
            vs.add(v)
        return vs

    def __len__(self):
        '''Return the number of entries in the dict at the current time.

        :returns: the length of the dict'''
        return len(self._now)

    def __getitem__(self, k: K) -> V:
        '''Retrieve the value associated with the given key at the current time.
        Raises a KeyError if the key is not defined (even if it has been, or will be,
        at other times).

        :param k: the key
        :returns: the value'''
        if k in self._now:
            (_, _, v) = self._dict[k][self._now[k]]
            return v
        else:
            t = self._time
            raise KeyError(f'No key {k} at time {t}')

    def __setitem__(self, k: K, v: V):
        '''Set the value associated with the given key at the current time,
        overwriting any current value.

        :param k: the key
        :param v: the value'''
        if k in self._now:
            (ct, _, _) = self._dict[k][self._now[k]]
            if ct == self._time:
                # update at the current time
                self._dict[k][self._now[k]] = (self._time, True, v)
            else:
                # update at a time after the last update, insert a new entry
                self._dict[k].insert(self._now[k] + 1, (self._time, True, v))
                self._now[k] += 1
        else:
            # new element (at this time)
            i = self._updateBefore(k)
            if i < 0:
                # globally new entry, add to the main dict
                self._dict[k] = [(self._time, True, v)]
                self._now[k] = 0
            else:
                # add an entry for this time
                self._dict[k].insert(i + 1, (self._time, True, v))
                self._now[k] += 1

    def __delitem__(self, k: K):
        '''Delete the mapping for the given key at the current time. This
        does not affect values at earlier times, or assignments in the future.

        :param k: the key'''
        if k not in self._now:
            t = self._time
            raise KeyError(f'No key {k} at time {t}')
        else:
            i = self._updateBefore(k)
            self._dict[k].insert(i + 1, (self._time, False, None))
            del self._now[k]


    # ---------- conversion ----------

    def asdict(self) -> Dict[K, V]:
        '''Return a snapshot at the current time as a dict.

        :returns: a dict'''
        d = dict()
        for k in self.keys():
            d[k] = self[k]
        return d


class TimedDict(Generic[K, V]):
    '''A dict whose entries are all keyed by time, allowing the contents
    of the dict to be accessed at any timestamp.

    The implementation is optimised for access, on the assumption that
    we'll access a substantial number of nodes at a given time:
    access is "sparse in time" but "dense in space".'''

    def __init__(self):
        self._dict: Dict[K, List[Tuple[float, bool, V]]] = dict()
        self._time: float = 0.0


    # ---------- dict interface ----------

    def __getitem__(self, t: float) -> 'TimedDictView':
        '''Retrieve a view of the dict at the given time. The
        view reflect all changes made to the dict up to and
        including that time.

        :param t: the time
        :returns: a view of the dict at that time'''
        return TimedDictView(self._dict, t)
