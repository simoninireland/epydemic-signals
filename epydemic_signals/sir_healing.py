# An SIR model that also stores the healing times at which nodes are removed
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

import sys
if sys.version_info >= (3, 8):
    from typing import Final
else:
    from typing_extensions import Final
from epydemic import SIR, Node


class HealingSIR(SIR):
    '''An :class:`SIR` model that also captures healing times of nodes, when they
    are removed.'''

    T_HEALING: Final[str] = None   #: State variable for healing time.

    def __init__(self):
        super().__init__()
        self.T_HEALING = self.stateVariable('tHealing')

    def markHealed(self, n: Node, t: float):
        '''Mark a node as healed (removed), storing the healing time.

        :param n: the node:param t: the time'''
        g = self.network()
        g.nodes[n][self.T_HEALING] = t

    def remove(self, t: float, n: Node):
        '''Save the healing time for the removed node.

        :param t: the simulation time
        :param n: the node'''
        super().remove(t, n)
        self.markHealed(n, t)
