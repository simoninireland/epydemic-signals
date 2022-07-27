# Tests of timed dicts
#
# Copyright (C) 2021--2022 Simon Dobson
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

from epydemic_signals import *
import epyc
import unittest
import networkx


class TimedDictTests(unittest.TestCase):

    def setUp(self):
        self._dict = TimedDict()

    def testEmpty(self):
        '''Test the empty timed dict.'''
        d = self._dict[0]
        self.assertCountEqual(d.keys(), [])
        self.assertCountEqual(d.values(), [])
        self.assertEqual(len(self._dict), 0)
        self.assertEqual(len(d), 0)
        self.assertCountEqual(self._dict.updates(), [])

    def testAdd(self):
        '''Test adding at the same time.'''
        d = self._dict[0]
        d['a'] = 10
        d['b'] = 20
        self.assertCountEqual(d.keys(), ['a', 'b'])
        self.assertCountEqual(d.values(), [10, 20])
        self.assertEqual(len(d), 2)

    def testAddSequential(self):
        '''Test we can add at different times.'''
        d0 = self._dict[0]
        d0['a'] = 10
        self.assertCountEqual(d0.keys(), ['a'])
        d1 = self._dict[1]
        self.assertCountEqual(d1.keys(), ['a'])
        d1['a'] = 20
        d1['b'] = 30
        self.assertCountEqual(d1.keys(), ['a', 'b'])

        d = self._dict[0]
        self.assertCountEqual(d.keys(), ['a'])
        self.assertEqual(d['a'], 10)

        d = self._dict[1]
        self.assertCountEqual(d.keys(), ['a', 'b'])
        self.assertEqual(d['a'], 20)
        self.assertEqual(d['b'], 30)

    def testMidTimes(self):
        '''Test we can access mid-times.'''
        d0 = self._dict[0]
        d0['a'] = 10
        d1 = self._dict[1]
        d1['a'] = 20
        d1['b'] = 30

        d = self._dict[0.5]
        self.assertCountEqual(d.keys(), ['a'])

    def testDeleteNever(self):
        '''Test we can silently delete an element that's never been added.'''
        d0 = self._dict[0]
        d0['a'] = 10
        del d0['c']

        d1 = self._dict[1]
        d1['b'] = 20
        del d1['c']

    def deleteNoLonger(self):
        '''Test we can re-delete an element we already deleted.'''
        d0 = self._dict[0]
        d0['a'] = 10
        d0['b'] = 20

        d1 = self._dict[1]
        del d1['b']
        del d1['b']

        d = self._dict[0]
        self.assertCountEqual(d.keys(), ['a', 'b'])

    def testAddDeleteSameTimeGlobal(self):
        '''Test we can add and delete in the same time.'''
        d = self._dict[0]
        d['a'] = 10
        del d['a']
        self.assertCountEqual(d.keys(), [])

    def testAddDeleteSameTimeLater(self):
        '''Test we can add and delete in the same time.'''
        d0 = self._dict[0]
        d0['a'] = 10

        d1 = self._dict[1]
        d1['a'] = 20
        del d1['a']
        self.assertCountEqual(d1.keys(), [])

        d = self._dict[0]
        self.assertCountEqual(d.keys(), ['a'])

    def testDeleteInTheMiddle(self):
        '''Test we can delete at an intermediate time.'''
        d0 = self._dict[0]
        d0['a'] = 10
        self.assertEqual(len(d0), 1)

        d1 = self._dict[1]
        d1['a'] = 20
        self.assertEqual(len(d1), 1)

        d2 = self._dict[2]
        d2['a'] = 30
        self.assertEqual(len(d2), 1)

        d = self._dict[1]
        del d['a']
        self.assertCountEqual(d.keys(), [])
        self.assertEqual(len(d), 0)

        d = self._dict[0]
        self.assertCountEqual(d.keys(), ['a'])
        self.assertEqual(d['a'], 10)

        d = self._dict[1]
        self.assertCountEqual(d.keys(), [])

        d = self._dict[2]
        self.assertCountEqual(d.keys(), ['a'])
        self.assertEqual(d['a'], 30)

    def testDeletionDeletes(self):
        '''Test we can't retrieve a value we've deleted.'''
        d0 = self._dict[0]
        d0['a'] = 10

        d1 = self._dict[1]
        del d1['a']
        with self.assertRaises(Exception):
            d1['a']

        d2 = self._dict[2]
        with self.assertRaises(Exception):
            d2['a']

    def testGetGlobal(self):
        '''Test we cam't get a value that's never been added.'''
        d0 = self._dict[0]
        d0['a'] = 10
        with self.assertRaises(Exception):
            d0['b']

        d1 = self._dict[1]
        with self.assertRaises(Exception):
            d1['b']

    def testGetDefault(self):
        '''Test getting default values.'''
        d = self._dict[0]
        self.assertEqual(d.get('a', 10), 10)

    def testUpdateSame(self):
        '''Test that new updates to the same value don't add a transition.'''
        d0 = self._dict[0]
        d0['a'] = 10

        d1 = self._dict[1]
        self.assertEqual(d1['a'], 10)
        d1['a'] = 10
        self.assertEqual(d1['a'], 10)
        self.assertCountEqual(self._dict.updates(), [0])

        d1['a'] = 20
        self.assertEqual(d1['a'], 20)
        self.assertCountEqual(self._dict.updates(), [0, 1])

    def testNewValueAfterDelete(self):
        '''Test updating back to a value after an intermediate deletion.'''
        d0 = self._dict[0]
        d0['a'] = 10
        self.assertCountEqual(self._dict.updates(), [0])

        d1 = self._dict[1]
        del d1['a']
        self.assertCountEqual(self._dict.updates(), [0, 1])

        # should add an update, since we're after a delete
        d2 = self._dict[2]
        d2['a'] = 10
        self.assertCountEqual(self._dict.updates(), [0, 1, 2])

        # shouldn't add an update
        d3 = self._dict[3]
        d3['a'] = 10
        self.assertCountEqual(self._dict.updates(), [0, 1, 2])

    def testSnapshot(self):
        '''Test we can snap the timed dict.'''
        d0 = self._dict[0]
        self.assertEqual(len(d0.asdict()), 0)
        d0['a'] = 10
        self.assertEqual(len(d0.asdict()), 1)
        self.assertEqual(d0.asdict()['a'], 10)

        d1 = self._dict[1]
        self.assertEqual(len(d1.asdict()), 1)
        d1['b'] = 20
        self.assertEqual(len(d1.asdict()), 2)
        self.assertEqual(d1.asdict()['a'], 10)
        self.assertEqual(d1.asdict()['b'], 20)

    def testArrayEmpty(self):
        '''Test we can snapshot as an array, with no keys requested.'''
        d = self._dict[0]
        a = d.asarray([])
        self.assertEqual(a.shape[0], 0)

    def testArray(self):
        '''Test we can snapshot as an array.'''
        d = self._dict[0]
        d['a'] = 1
        d['b'] = 2
        d['c'] = 3

        a = d.asarray(['a', 'b', 'c'])
        self.assertCountEqual(a, [1, 2, 3])

        # different order
        a = d.asarray(['b', 'c', 'a'])
        self.assertCountEqual(a, [2, 3, 1])

        # duplicates
        a = d.asarray(['b', 'c', 'b'])
        self.assertCountEqual(a, [2, 3, 2])

    def testArrayWithNodes(self):
        '''Test we can get the values associated with the nodes of a graph.'''
        d = self._dict[0]
        g = networkx.Graph()
        ns = [0, 1, 2, 3, 4, 5]
        g.add_nodes_from(ns)
        keys = [5, 4, 2, 0, 3, 1]
        vals = [n * 8 for n in keys]
        d.setFrom(keys, vals)

        a = d.asarray(ns)
        for i in range(len(ns)):
            self.assertEqual(a[i], i * 8)

    def testUpdates(self):
        '''Test retrieving update times.'''
        d = self._dict[0]
        d['a'] = 0
        d = self._dict[1]
        d['b'] = 20
        del d['a']
        d = self._dict[2]
        del d['b']
        d = self._dict[3]
        # do nothing at this latest time

        self.assertCountEqual(self._dict.updates(), [0, 1, 2])

    def testAllKeys(self):
        '''Test we can retrieve all keys over all times.'''
        d = self._dict[0]
        d['a'] = 0
        d = self._dict[1]
        d['b'] = 20
        del d['a']
        d = self._dict[2]
        del d['b']
        d = self._dict[3]

        self.assertCountEqual(self._dict.keysAtSomeTime(), ['a', 'b'])

    def testAllValues(self):
        '''Test we can retrieve all values over all times.'''
        d = self._dict[0]
        d['a'] = 0
        d = self._dict[1]
        d['b'] = 20
        d['b'] = 30      # this overwrites the previous value at this time
        del d['a']
        d = self._dict[2]
        del d['b']
        d = self._dict[3]

        self.assertCountEqual(self._dict.valuesAtSomeTime(), [0, 30])

    def testSetFrom(self):
        '''Test we can set a selection of values in one go.'''
        d = self._dict[0]
        keys = ['a', 'b', 'c']
        vals = [25, 35, 45]
        d.setFrom(keys, vals)
        self.assertCountEqual(d.keys(), keys)
        self.assertCountEqual(d.values(), vals)
        for i in range(len(keys)):
            self.assertEqual(d[keys[i]], vals[i])

    def testSetFromFewerKeys(self):
        '''Test we can't have more values than keys.'''
        d = self._dict[0]
        keys = ['a', 'b']
        vals = [25, 35, 45]
        with self.assertRaises(ValueError):
            d.setFrom(keys, vals)

    def testSetFromFewerValues(self):
        '''Test we can't have more keys than values.'''
        d = self._dict[0]
        keys = ['a', 'b', 'c']
        vals = [25, 35]
        with self.assertRaises(ValueError):
            d.setFrom(keys, vals)

    def testDeleteFromNow(self):
        '''Test we can delete a set of signal values at the current moment.'''
        d = self._dict[0]
        keys = ['a', 'b', 'c']
        vals = [25, 35, 45]
        d.setFrom(keys, vals)
        d.deleteFrom(['a', 'c'])
        self.assertCountEqual(d.keys(), ['b'])

    def testDeleteFromLater(self):
        '''Test we can delete a set of signal values at a later moment.'''
        d = self._dict[0]
        keys = ['a', 'b', 'c']
        vals = [25, 35, 45]
        d.setFrom(keys, vals)

        d1 = self._dict[1]
        d1.deleteFrom(['a', 'c'])
        self.assertCountEqual(d1.keys(), ['b'])


if __name__ == '__main__':
    unittest.main()
