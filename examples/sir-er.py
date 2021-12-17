# Example signals over a moderately-sized ER network
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

import pickle
from networkx import Graph
from epyc import Lab, Experiment
from epydemic import ERNetwork, SIR
from epydemic_signals import OneInfectionSIR, StochasticSignalDynamics, Signal
from epydemic_signals import CompartmentSignalGenerator, SIRProgressSignalGenerator

lab = Lab()

# network and disease parameters
lab[ERNetwork.N] = 10000
lab[ERNetwork.KMEAN] = 40
lab[SIR.P_INFECTED] = 0.001
lab[SIR.P_REMOVE] = 0.002
lab[SIR.P_INFECT] = 0.00015

# model and experiment
sir = SIR()
e = StochasticSignalDynamics(sir, ERNetwork())

# add signal generators
compartments = Signal()
gen = CompartmentSignalGenerator(sir, compartments)
e.addSignalGenerator(gen)
progress = Signal()
gen = SIRProgressSignalGenerator(sir, progress)
e.addSignalGenerator(gen)

# run the experiment (don't care about the actual results, just the signals)
lab.runExperiment(e)
rc = lab.results()[0]
if not rc[Experiment.METADATA][Experiment.STATUS]:
    print(rc[Experiment.METADATA][Experiment.TRACEBACK])

# save the signals
with open('sir-er-10000.pickle', 'wb') as fh:
    pickle.dump(progress, fh)
    pickle.dump(compartments, fh)
