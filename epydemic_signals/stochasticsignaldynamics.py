# A stochasitc dynamics with signal generation
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

from typing import Union
from networkx import Graph
from epydemic import StochasitcDynamics, Process, NetworkGenerator
from epydemic_signals import SignalDynamics


class StochasticSignalDynamics(StochasticDynamics, SignalDynamics):
    '''A stochastic (Gillespie) dynamics that passes all events to
    a signal generator.

    :param p: the process to run
    :param g: network or network generator (optional, can be provided later)'''

    def __init__(self, p: Process, g: Union[Graph, NetworkGenerator] = None):
        super().__init__(p, g)