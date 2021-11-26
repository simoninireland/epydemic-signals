# Initialistion for the epydemic_signals package
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

# Utilities
from .timeddict import TimedDict

# Epidemic process variants
from .sir_healing import HealingSIR
from .sir_healing_one import OneInfectionSIR

# Analysis mixins
from .hitting_healing import HittingHealingTimes

# Signals
from .signal import Signal
from .progress_signal import SIRProgressSignal
