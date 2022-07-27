# Dicts with changes logged in time
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

from numpy import array, zeros
from typing import Generic, TypeVar,Union, Dict, Tuple, List, Iterable, cast


# The top-level class for this code is the TimedDict. This is however a
# very thin wrapper onto the TimedDictView, which is the class that
# actually implements the sequence of diffs used to represent time-varying
# mappings.


# Type variables for dict keys and values
K = TypeVar('K')
V = TypeVar('V')

# Type variables for the zipFail iterator combinator
X = TypeVar('X')
Y = TypeVar('Y')


class TimedDictView(Generic[K, V]):
    '''A view of a timed dict snapped at a particular time.

    :param d: the timed dict's diff structure
    :param t: the time'''

    # TODO: This needs to be optimised to only look for the "now" state of nodes
    # that are actually accessed.

    def __init__(self, d: Dict[K, List[Tuple[float, bool, V]]], t : float):
        self._dict = d                     # dict from key to diff list
        self._time = t                     # projection time
        self._now: Dict[K, int] = dict()   # dict from key to index in diff list of last update in diff list
        self._project()


    # ---------- projection ----------

    def _project(self):
        '''Project-out the values in the dict at the current time.'''
        self._now = dict()
        for k in self._dict:
            i = self._updateBefore(k)
            if i >= 0:
                (_, b, _) = self._dict[k][i]
                if b:
                    # update was a set, include this key in the projection
                    self._now[k] = i
                else:
                    # update was a delete, don't include the key
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

    def _hasValueNow(self, k):
        '''Test whether a key currently has a value, meaning that it has
        been assigned to at some earrlier time and has not been subsequently deleted.

        :param k: the key
        :returns: True if the key has a value'''
        return k in self._now


    # ---------- dict interface ----------

    def keys(self) -> Iterable[K]:
        '''Return the keys in the dict at the current time.

        :returns: a list of keys'''
        return self._now.keys()

    def __contains__(self, k: K) -> bool:
        '''Test whether the given k is defined at the current time.

        :param k: the key
        :returns: True if the key is in the dict at the current time'''
        return k in self.keys()

    def values(self) -> Iterable[V]:
        '''Return a list of values in the dict at the current time.

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
        at other times). This is the dict form of :meth:`get`.

        :param k: the key
        :returns: the value'''
        return self.get(k)

    def __setitem__(self, k: K, v: V):
        '''Set the value associated with the given key at the current time,
        overwriting any current value.

        :param k: the key
        :param v: the value'''
        #t = self._time
        if self._hasValueNow(k):
            (ct, up, pv) = self._dict[k][self._now[k]]
            if ct == self._time:
                # update at the current time
                #print(f'overwritten {k}={v} at time {ct}')
                self._dict[k][self._now[k]] = (self._time, True, v)
            else:
                # only perform an update if the value differs from the last one
                if up and (pv != v):
                    # update at a time after the last update, insert a new entry
                    #print(f'changed {k}={v} at time {t}')
                    self._dict[k].insert(self._now[k] + 1, (self._time, True, v))
                    self._now[k] += 1
        else:
            # new element (at this time)
            i = self._updateBefore(k)
            if i < 0:
                # globally new entry, add to the main dict
                #print(f'initial {k}={v} at time {t}')
                self._dict[k] = [(self._time, True, v)]
                self._now[k] = 0
            else:
                # new element after a deletion, add an entry
                #print(f'new {k}={v} at time {t}')
                self._dict[k].insert(i + 1, (self._time, True, v))
                self._now[k] = i + 1

    @staticmethod
    def zipFail(v1s: Iterable[X], v2s: Iterable[Y]) -> Iterable[Tuple[X, Y]]:
        '''Return an iterator over a pair of child iterators that returns
        corresponding values from each, raising an exception if either of them
        runs out before the other. This contrasts with Python's `zip` function
        (which exhausts silently) and with `itertools.zip_longest` (which pads the
        shorter sequence).

        :param v1s the first iterable
        :param v2s: the second iterable
        :returns: an iterator of pairs'''
        i1 = iter(v1s)
        i2 = iter(v2s)
        i1_exhausted = False
        while True:
            # get a value from i1
            try:
                a = next(i1)
            except StopIteration:
                i1_exhausted = True

            # get a value from i2
            try:
                b = next(i2)
                if i1_exhausted:
                    # i1 finished before i2
                    raise ValueError('First iterator exhausted prematurely')
                else:
                    yield (a, b)
            except StopIteration:
                if i1_exhausted:
                    # both iterators exhausted together, so we're finished
                    return
                else:
                    # i2 finished before i1
                    raise ValueError('Second iterator exhausted prematurely')

    def setFrom(self, ks: Iterable[K], vs: Iterable[V]):
        '''Set the value at several keys. The keys and values are passed
        as two iterables, typically lists. These should be the same length: if
        not, the pairs that *can* be added, *will* be added, and an exception
        will then be raised.

        :param ks: the list of keys
        :param ss: the list of values'''
        for (k, v) in TimedDictView.zipFail(ks, vs):
            self[k] = v

    def __delitem__(self, k: K):
        '''Delete the mapping for the given key at the current time. This
        does not affect values at earlier times, or assignments in the future.
        It is silent if there is no entry for the given key at the given time.

        :param k: the key'''
        if self._hasValueNow(k):
            i = self._updateBefore(k)
            self._dict[k].insert(i + 1, (self._time, False, None))
            del self._now[k]

    def deleteFrom(self, ks: Iterable[K]):
        '''Delete the values associated with several keys.

        :param ks: the list of keys'''
        for k in ks:
            del self[k]

    def get(self, k: K, default: V = None) -> V:
        '''Get the value of the given key, returning the default value if
        the key doesn't have a value. Raise a KeyError if no default value
        is given and the key is missing.

        :param k: the key
        :param default: (optional) default value
        :returns: the key value of the default'''
        if self._hasValueNow(k):
            # key has a value, return it
            (_, _, v) = self._dict[k][self._now[k]]
            return v
        elif default is not None:
            # no value, return the default
            return default
        else:
            # no value and no default, raise an exception
            t = self._time
            raise KeyError(f'No key {k} at time {t}')


    # ---------- conversion ----------

    def asdict(self) -> Dict[K, V]:
        '''Return a snapshot at the current time as a dict.

        :returns: a dict'''
        d = dict()
        for k in self.keys():
            d[k] = self[k]
        return d

    def asarray(self, ks: Iterable[K] = None) -> array:
        '''Return a snapshot at the current time as a `numpy` array, with
        the order of the values being given by the list of keys. If no
        keys are given then all the keys with current values are used, in
        key order: this isn't likely to be too meaningful.

        :param ks: (optional) the keys (defaults to all)
        :returns: an array'''
        if ks is None:
            ks = list(self.keys())
        a = zeros(len(ks))
        i = 0
        for k in ks:
            a[i] = self[k]
            i += 1
        return a


class TimedDict(Generic[K, V]):
    '''A dict whose entries are all keyed by time, allowing the contents
    of the dict to be accessed at any timestamp.

    The implementation is optimised for access, on the assumption that
    we'll access a substantial number of nodes at a given time:
    access is "sparse in time" but "dense in space".'''

    def __init__(self):
        self._dict: Dict[K, List[Tuple[float, bool, V]]] = dict()
        self._time: float = 0.0


    # ---------- access ----------

    def updates(self) -> Iterable[float]:
        '''Return a list of update times, in ascending order.
        These are the points at which the signal changes in some way.

        The signal can be queried at any time, not just these, so they're
        not "keys" in any sense. However, the signal is guaranteed to have
        changed in some way between adjacent updates, so they do represent
        the "meaningful changes" to keys.

        :returns: a list of times'''
        ts = set()
        for k in self._dict:
            us = self._dict[k]
            for (t, _, _) in us:
                ts.add(t)
        sts = list(ts)
        sts.sort()
        return sts

    def keysAtSomeTime(self) -> Iterable[K]:
        '''Return the set of keys that appear at some time in the dict.

        :returns: a set of keys'''
        return self._dict.keys()

    def valuesAtSomeTime(self) -> Iterable[V]:
        '''Return a set of the values that have been assigned to some key
        at some time.

        The values are all the values that can be retrieved at *some* time
        from the dict. This *excludes* values that have been updated during
        a given time. For example, the code:

        .. code-block:: python

           td = TimedDict()
           d = td[0.0]
           d['a'] = 1
           d['b'] = 2
           d['c'] = 5

           d = td[1.0]
           d['b'] = 3
           d['b'] = 4

           print(td.valuesAtSomeTime())

        will print "1 2 5 4" (in some undefined order), as 2 and 4 could both
        be retrieved at some time but 3 couldn't be found in the dict at any
        time (since it was overwritten by 4).

        :returns: a set of values

        '''
        vs = set()
        for k in self._dict:
            us = self._dict[k]
            for (_, u, v) in us:
                if u:
                    # we're only concerned with updates, not deletions
                    vs.add(v)
        return vs


    # ---------- dict interface ----------

    def __getitem__(self, t: float) -> 'TimedDictView':
        '''Retrieve a view of the dict at the given time. The
        view reflect all changes made to the dict up to and
        including that time.

        :param t: the time
        :returns: a view of the dict at that time'''
        return TimedDictView(self._dict, t)

    def __len__(self):
        '''Return the number of transition points in the dict.
        This is not the number of keys per se, since the dict
        can be queried for *any* time; rather it is the length
        of the set returned by :meth:`updates`.

        :returns: the number of update times'''
        ts = set()
        for k in self._dict:
            us = self._dict[k]
            for (t, _, _) in us:
                ts.add(t)
        return len(ts)
