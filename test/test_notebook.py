# Tests integration with epyc notebooks
#
# Copyright (C) 2021-2022 Simon Dobson
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
from epydemic import *
from epyc import *
import unittest
import networkx


class BondPercolationSignal(SignalExperiment, BondPercolation):
    '''A bond percolation experiment with a signal generator tap.'''

    def __init__(self, g, samples = None):
        super().__init__(g, samples)


class StaticGenerator(SignalGenerator):
    '''A signal generator that generates a worthwhile but simple signal.'''

    def __init__(self):
        super().__init__(Signal(name='test'))

    def event(self, t, etype, e):
        if etype == BondPercolation.SAMPLE:
            s_t = self.signal()[t]
            for n in self.network().nodes():
                s_t[n] = n * t


class NotebookTests(unittest.TestCase):

    def setUp(self):
        self._g = networkx.Graph()
        self._g.add_nodes_from([0, 1, 2, 3, 4, 5, 6])
        self._g.add_edges_from([(0, 1), (1, 2), (3, 0), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (4, 6)])
        self._signal = Signal()

    def testLab(self):
        '''Test the use of signals with a lab notebook, the "normal" use case.'''
        lab = Lab()
        nb = lab.notebook()
        lab['experiment'] = 1

        e = BondPercolationSignal(FixedNetwork(self._g))
        gen = StaticGenerator()
        e.addSignalGenerator(gen)

        rc = lab.runExperiment(e)
        df = nb.dataframe()
        self.assertEqual(len(df), 1)

        # check we have all the series we expect for the signal
        signal = gen.signal()
        (tn, en, vn) = BondPercolationSignal.signalSeries(signal)
        for k in [tn, en, vn]:
            self.assertIn(k, df.keys())

        # check the series agree with the signal
        for (t, n, v) in zip(df.iloc[0][tn], df.iloc[0][en], df.iloc[0][vn]):
            self.assertEqual(signal[t][n], v)


if __name__ == '__main__':
    unittest.main()
