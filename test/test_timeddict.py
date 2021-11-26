# Tests of timed dicts
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
        '''Test we can't delete an element that's never been added.'''
        d0 = self._dict[0]
        d0['a'] = 10
        with self.assertRaises(Exception):
            del d0['c']

        d1 = self._dict[1]
        d1['b'] = 20
        with self.assertRaises(Exception):
            del d1['c']

    def deleteNoLonger(self):
        '''Test we can't delete an element we already deleted.'''
        d0 = self._dict[0]
        d0['a'] = 10
        d0['b'] = 20

        d1 = self._dict[1]
        del d1['b']
        with self.assertRaises(Exception):
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
        d['b'] = 30
        del d['a']
        d = self._dict[2]
        del d['b']
        d = self._dict[3]

        self.assertCountEqual(self._dict.valuesAtSomeTime(), [0, 30])


if __name__ == '__main__':
    unittest.main()
