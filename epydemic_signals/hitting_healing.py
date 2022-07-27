# Process to collate hittign and healing times
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

import sys
from typing import Dict, List, Tuple, Any, Union, cast
if sys.version_info >= (3, 8):
    from typing import Final
else:
    from typing_extensions import Final
from networkx import get_node_attributes
from pandas import DataFrame
from epyc import ResultsDict, Experiment
from epydemic import SIR, Node, Process


class HittingHealingTimes(Process):
    '''Extract the sequence of hitting and healing times for an
    :class:`SIR` epidemic and store them in the results. This is an
    analysis mixin intended to be added to the end of a
    :class:`ProcessSequence` to be run after an epidemic process has run
    on the network.

    The process assumes that the hitting and healing times are stored
    as state on network nodes, in the :attr:`SIR.T_HITTING` and
    :attr:`HealingSIR.T_HEALING` attributes respectively. These are
    assembled into four sequences:

    - :attr:`HITTING_TIMES`: the hitting times
    - :attr:`HITTING_TIMES_NODES`: the nodes in infection order
    - :attr:`HEALING_TIMES`: the healing times
    - :attr:`HEALING_TIMES_NODES`: the nodes in removal order

    These can then be firtherr processed by :meth:`ttimeline` to construct
    a single event timeline.'''

    # Experimental results
    HITTING_TIMES: Final[str] = 'epydemic.hitting-healing.hitting-times'
    HITTING_TIMES_NODES: Final[str] = 'epydemic.hitting-healing.hitting-times-nodes'
    HEALING_TIMES: Final[str] = 'epydemic.hitting-healing.healing-times'
    HEALING_TIMES_NODES: Final[str] = 'epydemic.hitting-healing.healing-times-nodes'

    @classmethod
    def timeline(cls, ts: Union[DataFrame, ResultsDict]) -> List[Tuple[float, str, Node]]:
        '''Interleave hitting and healing to form a single timeline, in ascending
        order of time (essentially the merge step of a merge sort). Each element
        of the timeline is a triple consisting of the event time, a string describing
        the event type, and the affected node. The event types are represented by
        :attr:`SIR.INFECTED` (for infections) and :attr:`SIR.REMOVED` for healing (removal).

        :param ts: a DataFrame or results dict
        :returns: a list of (time, event, node) triples'''

        # extract the raw data
        if isinstance(ts, DataFrame):
            # dataframe, extract columns
            df = cast(DataFrame, ts)
            hitting_ts = list(df[cls.HITTING_TIMES].iloc[0])
            hitting_ns = list(df[cls.HITTING_TIMES_NODES].iloc[0])
            healing_ts = list(df[cls.HEALING_TIMES].iloc[0])
            healing_ns = list(df[cls.HEALING_TIMES_NODES].iloc[0])
        else:
            # results dict
            rc = cast(ResultsDict, ts)
            res = rc[Experiment.RESULTS]
            hitting_ts = res.get(cls.HITTING_TIMES, [])
            hitting_ns = res.get(cls.HITTING_TIMES_NODES, [])
            healing_ts = res.get(cls.HEALING_TIMES, [])
            healing_ns = res.get(cls.HEALING_TIMES_NODES, [])

        # merge events
        evs = []
        while len(hitting_ts) > 0 and len(healing_ts) > 0:
            if hitting_ts[0] < healing_ts[0]:
                evs.append((hitting_ts.pop(0), SIR.INFECTED, (hitting_ns.pop(0), None)))
            else:
                evs.append((healing_ts.pop(0), SIR.REMOVED, (healing_ns.pop(0), None)))
        while len(hitting_ts) > 0:
            evs.append((hitting_ts.pop(0), SIR.INFECTED, hitting_ns.pop(0)))
        while len(healing_ts) > 0:
            evs.append((healing_ts.pop(0), SIR.REMOVED, healing_ns.pop(0)))

        return evs

    def __init__(self, p: Process):
        super().__init__()
        self._disease = p

    def results(self) -> Dict[str, Any]:
        '''Add the sequence of healing and hitting times to the results.

        :returns: the experimental results'''
        res = super().results()
        g = self.network()

        # extract mappings of node to hitting and healiing times
        hitting_ns = get_node_attributes(g, self._disease.T_HITTING)
        healing_ns = get_node_attributes(g, self._disease.T_HEALING)

        # invert the mappings (safe since we know they're one-to-one) and sort
        hitting_tns = [(t, n) for (n, t) in hitting_ns.items()]
        hitting_tns.sort(key=lambda tn: tn[0])
        healing_tns = [(t, n) for (n, t) in healing_ns.items()]
        healing_tns.sort(key=lambda tn: tn[0])

        # add to the results
        res[self.HITTING_TIMES] = [t for (t, _) in hitting_tns]
        res[self.HITTING_TIMES_NODES] = [n for (_, n) in hitting_tns]
        res[self.HEALING_TIMES] = [t for (t, _) in healing_tns]
        res[self.HEALING_TIMES_NODES] = [n for (_, n) in healing_tns]
        return res
