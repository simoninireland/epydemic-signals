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

import os
from tempfile import NamedTemporaryFile
from epydemic_signals import *
from epydemic import *
from epyc import *
import unittest
import networkx


class StaticGenerator(SignalGenerator):
    '''A signal generator that generates a worthwhile but simple signal.'''

    SIGNALNAME = 'test'

    def __init__(self, name = None, multiplier = 1):
        if name is None:
            name = self.SIGNALNAME
        super().__init__(Signal(name=name))
        self._m = multiplier

    def event(self, t, etype, e):
        s_t = self.signal()[t]
        for n in self.network().nodes():
            s_t[n] = n * t * self._m


class NotebookTests(unittest.TestCase):

    def setUp(self):
        tf = NamedTemporaryFile()
        tf.close()
        self._fn = tf.name
        #self._fn = 'test.h5'

        self._g = networkx.complete_graph(50)

        self._nb = LabNotebook()
        self._lab = Lab(notebook=self._nb)
        self._lab[SIR.P_INFECTED] = 0.2
        self._lab[SIR.P_REMOVE] = 1.0
        self._lab[SIR.P_INFECT] = 0.8
        self._model = SIR()

    def tearDown( self ):
        try:
            os.remove(self._fn)
            #pass
        except OSError:
            pass

    def testLab(self):
        '''Test the use of signals with a lab notebook, the "normal" use case.'''
        gen = StaticGenerator()
        e = StochasticSignalDynamics(self._model, self._g)
        e.attachSignalGenerator(gen, self._model)

        rc = self._lab.runExperiment(e)
        df = self._nb.dataframe()
        self.assertEqual(len(df), 1)

        # check we have all the series we expect for the signal
        signal = gen.signal()
        (tn, en, vn) = StochasticSignalDynamics.signalSeries(signal)
        for k in [tn, en, vn]:
            self.assertIn(k, df.keys())

    def testNoSignal(self):
        '''Test that a null signal isn't added.'''
        gen = StaticGenerator()
        e = StochasticSignalDynamics(self._model, self._g)
        e.attachSignalGenerator(gen, self._model)

        self._lab[SIR.P_INFECTED] = 0.0  # no epidemic occurs

        rc = self._lab.runExperiment(e)
        df = self._nb.dataframe()
        self.assertEqual(len(df), 1)

        self.assertNotIn(SignalExperiment.SIGNALS, df)

    def testAllSignalsNone(self):
        '''Test extracting a signal from a series that has none.'''
        gen = StaticGenerator()
        e = StochasticSignalDynamics(self._model, self._g)  # no signal generator attached

        rc = self._lab.runExperiment(e)
        df = self._nb.dataframe()
        self.assertEqual(len(df), 1)
        self.assertCountEqual(SignalExperiment.signals(df.iloc[0]), [])

    def testAllSignalsOne(self):
        '''Test extracting a signal from a series that has one.'''
        gen = StaticGenerator()
        e = StochasticSignalDynamics(self._model, self._g)
        e.attachSignalGenerator(gen, self._model)

        rc = self._lab.runExperiment(e)
        df = self._nb.dataframe()
        self.assertEqual(len(df), 1)
        self.assertEqual(len(SignalExperiment.signals(df.iloc[0])), 1)

    def testAllSignalsTwo(self):
        '''Test extracting a signal from a series that has two.'''
        gen1 = StaticGenerator()
        gen2 = StaticGenerator(name='ttt', multiplier=10)
        e = StochasticSignalDynamics(self._model, self._g)
        e.attachSignalGenerator(gen1, self._model)
        e.attachSignalGenerator(gen2, self._model)

        rc = self._lab.runExperiment(e)
        df = self._nb.dataframe()
        self.assertEqual(len(df), 1)
        ss = SignalExperiment.signals(df.iloc[0])
        self.assertEqual(len(ss), 2)

        s1 = ss[0]
        s2 = ss[1]
        if s1.name() == 'ttt':
            s1, s2 = s2, s1      # normalise signal order
        for t in s1.transitions():
            s1_t = s1[t]
            for n in s1_t.keys():
                self.assertEqual(s1_t[n], n * t)
        for t in s2.transitions():
            s2_t = s2[t]
            for n in s2_t.keys():
                self.assertEqual(s2_t[n], n * t * 10)


    # ---------- Different notebook storage formats ----------

    def signalFromResults(self):
        '''Basic code for notebook tests.'''
        sir = SIR()
        gen = StaticGenerator()
        e = StochasticSignalDynamics(sir, self._g)
        e.attachSignalGenerator(gen, sir)

        self._lab[SIR.P_INFECTED] = 0.2
        self._lab[SIR.P_REMOVE] = 1.0
        self._lab[SIR.P_INFECT] = 1.0

        rc = self._lab.runExperiment(e)
        s1 = gen.signal()
        df = self._nb.dataframe()

        # check we can reconstruct the signal
        s2 = Signal(name=StaticGenerator.SIGNALNAME)
        s2.fromSeries(df.iloc[0])

        # check the original and reconstructed signal agree
        for t in s1.transitions():
            s1_t = s1[t]
            s2_t = s2[t]
            for n in s1_t.keys():
                self.assertEqual(s1_t[n], s2_t[n])

    def testSignalFromResults(self):
        '''Test we can re-construct a signal from a notebook.'''
        self._nb = LabNotebook()
        self._lab = Lab(notebook=self._nb)
        self.signalFromResults()

    def testJSON(self):
        '''Test we can save and restore from a JSON notebook.'''
        self._nb = JSONLabNotebook(self._fn)
        self._lab = Lab(notebook=self._nb)
        self.signalFromResults()

    def testHDF5(self):
        '''Test we can save and restore from an HDF5 notebook.'''
        self._nb = HDF5LabNotebook(self._fn)
        self._lab = Lab(notebook=self._nb)
        self.signalFromResults()


if __name__ == '__main__':
    unittest.main()
