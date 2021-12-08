# Extract a progress signal from an SIR epidemic
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

from heapq import heappush, heappop
from typing import Dict, Any, List, Tuple, cast
from networkx import Graph, single_source_shortest_path
from epydemic import Node, Edge, SIR, Process, CompartmentedModel
from epydemic_signals import Signal, SignalGenerator


class SIRProgressSignalGenerator(SignalGenerator):
    '''Create the progress signal for an SIR epidemic.

    The progress signal captures all the dynamics of an SIR epidemic into
    a single function. At a time :math:`t` the signnal is:

    - on infected nodes, 0
    - on susceptible nodes, the length of the shortest path from the node
      to an infected node where the path only traverses susceptible nodes
    - on removed nodes, zero minus the length of the shortest path to an infected node
      traversing susceptible or removed nodes

    Put another way, the signal on a susceptible node indicates the number
    of hops needed at a minimum to infect that node, which the signal on
    an infected node indicates the number of hops away from the nearest source
    of infection.

    :param p: the process
    :param s: the signal
    '''

    def __init__(self, p: Process, s: Signal = None):
        super().__init__(p, s)
        self._inf = None
        self._compartment = dict()
        self._boundary = dict()
        self._coboundary_S = dict()
        self._coboundary_R = dict()

        # register the event handlers
        self.addEventTypeHandler(SIR.INFECTED, self.infect)
        self.addEventTypeHandler(SIR.REMOVED, self.remove)

    def setUp(self, g: Graph):
        '''Capture the initial state of the network as susceptible, infected,
        and removed sets.

        :param g: the network'''
        super().setUp(g)

        s = self.signal()
        g = self.network()
        p = self.process()
        signal = s[0.0]

        # extract the initial state and signal
        self._inf = g.order() + 1           # a distance longer than the longest possible path
        self._compartment = dict()
        self._compartment[SIR.SUSCEPTIBLE] = set()
        self._compartment[SIR.INFECTED] = set()
        self._compartment[SIR.REMOVED] = set()
        cm = cast(CompartmentedModel, p)
        for n in g.nodes():
            # grab initial compartment
            self._compartment[cm.getCompartment(n)].add(n)

            # signal is initially infinite everywhere
            signal[n] = self._inf
        if len(self._compartment[SIR.REMOVED]) > 0:
            # don't handle initial removeds in the population for now
            raise ValueError('Initial network contains removed nodes')

        # compute the initial signal at t=0
        #print('initial infecteds')
        for s in self._compartment[SIR.INFECTED]:
            signal[s] = 0
            distance = [(0, s)]
            visited = set()
            self._coboundary_S[s] = set()
            self._coboundary_R[s] = set()
            while len(distance) > 0:
                (_, n) = heappop(distance)
                if n not in visited:
                    #print(f'visit {n}')
                    visited.add(n)
                    d = 1 + signal[n]
                    for m in g.neighbors(n):
                        if m not in visited:
                            if m in self._compartment[SIR.SUSCEPTIBLE]:
                                if d < signal[m]:
                                    # update the signal
                                    signal[m] = d
                                    heappush(distance, (d, m))
                                    #print(f'propose {m} distance {d}')

                                    # update the boundary
                                    if m in self._boundary:
                                        self._coboundary_S[self._boundary[m]].remove(m)
                                    self._boundary[m] = s
                                    self._coboundary_S[s].add(m)
                                    #print(f'Sus boundary of {m} now {s}')
                                else:
                                    # prune the tree
                                    #print(f'prune {m}')
                                    #visited.add(m)
                                    pass
                            else:
                                # prune the tree
                                #print(f'prune {m}')
                                visited.add(m)
        #print('initial signal')
        #for n in g.nodes():
        #    print(n, signal[n])

    def _shortestPath(self, n, target, onpath):
        '''Return the length of the shortest path from the node to a
        node in the target set, traversing only nodes in the path set.

        :param n: the node
        :param target: the compartment of the target set
        :param onpath: the compartments of nodes included in the path
        :returns: the node and the shortest path, or None if there is no path'''
        g = self.network()
        distance = [(0, n)]
        visited = set()
        while len(distance) > 0:
            (d, n) = heappop(distance)
            visited.add(n)
            dprime = d + 1
            ms = g.neighbors(n)
            for m in ms:
                if m not in visited:
                    if m in self._compartment[target]:
                        # found a node in the target set, return the distance
                        return (m, dprime)
                    else:
                        for c in onpath:
                            if m in self._compartment[c]:
                                # m can be on the path, search it
                                heappush(distance, (dprime, m))
                                break

                        # if we get here, prune the tree of this node
                        visited.add(m)
        return None

    def infect(self, t: float, e: Edge):
        '''Adjust the signal for an infection event.

        :param t: the event time
        :param e: the SI edge the infection passed over'''
        (s, _) = e
        g = self.network()
        signal = self.signal()[t]
        #print('infect', s)

        # update state
        #print('Phase I-1')
        self._compartment[SIR.SUSCEPTIBLE].remove(s)
        if s in self._boundary:
            # s has a boundary, remove it from that node's co-boundary
            self._coboundary_S[self._boundary[s]].remove(s)
            del self._boundary[s]

            # (The only way s will *not* have a boundary is if the initial
            # state of the network was all susceptibles with no infecteds.
            # It might be worth handling this as a special case?)
        self._compartment[SIR.INFECTED].add(s)

        # set signal at s
        signal[s] = 0

        # iterate all susceptible and removed nodes updating signal as the
        # shortest path length to an infected node
        # This is a modified Dijkstra's algorithm that prunes the
        # search tree if it encounters a node whose proposed distance
        # is greater than the distance it has in the signal already
        #print('Phase I-2')
        distance = [(0, s)]
        visited = set()
        self._coboundary_S[s] = set()
        self._coboundary_R[s] = set()
        while len(distance) > 0:
            (_, n) = heappop(distance)
            if n not in visited:
                #print(f'visit {n}')
                visited.add(n)
                d = 1 + signal[n]
                for m in g.neighbors(n):
                    if m not in visited:
                        if m in self._compartment[SIR.SUSCEPTIBLE]:
                            if d < signal[m]:
                                # update the signal
                                signal[m] = d
                                heappush(distance, (d, m))
                                #print(f'propose {m} distance {d}')

                                # update the boundary
                                if m in self._boundary:
                                    self._coboundary_S[self._boundary[m]].remove(m)
                                self._boundary[m] = s
                                self._coboundary_S[s].add(m)
                                #print(f'Sus boundary of {m} now {s}')
                            else:
                                # prune the tree
                                #print(f'prune {m}')
                                #visited.add(m)
                                pass
                        elif m in self._compartment[SIR.REMOVED]:
                            if -d > signal[m]:
                                # update the signal
                                #print('update', m, signal[m], -d)
                                signal[m] = -d
                                heappush(distance, (d, m))

                                # update the boundary
                                self._coboundary_R[self._boundary[m]].remove(m)
                                self._boundary[m] = s
                                self._coboundary_R[s].add(m)
                                #print(f'Rem boundary of {m} now {s}')
                            else:
                                #visited.add(m)
                                pass
                        else:
                            # prune the tree
                            #print(f'prune {m}')
                            visited.add(m)
        #print(f'Sus coboundary of {s} now', self._coboundary_S[s], 'signal', signal[s])
        #print(f'Rem coboundary of {s} now', self._coboundary_R[s])

    def remove(self, t: float, s: Node):
        '''Adjust the signal for a removal event.

        :param t: the event time
        :param n: the node'''
        #print(f'remove {s}')
        g = self.network()
        signal = self._signal[t]

        # removal event, update state
        self._compartment[SIR.INFECTED].remove(s)
        self._compartment[SIR.REMOVED].add(s)

        # re-compute all susceptible distances affected by our removal
        #print('Phase R-1')
        for q in self._coboundary_S[s]:
            sp = self._shortestPath(q, SIR.INFECTED, [SIR.SUSCEPTIBLE])
            if sp is None:
                # no infected nodes found, set to infinity
                #print('no infected left')
                signal[s] = self._inf
            else:
                (n, d) = sp

                # update signal at this node if needed
                if d != signal[q]:
                    if d < signal[q]:
                        raise ValueError('Signal at {q} got smaller {before} {after}???'.format(q=q, after=d, before=signal[q]))
                    signal[q] = d
                    #print(f'update sus distance {q} {d}')
                else:
                    #print(f'sus distance {q} unchanged {d}')
                    pass
                self._boundary[q] = n
                self._coboundary_S[n].add(q)
        del self._coboundary_S[s]

        # find distance from removed node to boundary
        #print('Phase R-2')
        sp = self._shortestPath(s, SIR.INFECTED, [SIR.SUSCEPTIBLE, SIR.REMOVED])
        if sp is None:
            # no infected nodes found, set to minus infinity
            #print('no infected left')
            signal[s] = -self._inf
        else:
            (n, d) = sp

            # update signal and boundary
            signal[s] = -d
            #print('update to boundary', s, signal[s])
            self._boundary[s] = n
            self._coboundary_R[n].add(s)

        # update the signal for all other removed nodes affected by our removal
        #print('Phase R-3')
        for q in self._coboundary_R[s]:
            sp = self._shortestPath(q, SIR.INFECTED, [SIR.SUSCEPTIBLE, SIR.REMOVED])
            if sp is None:
                # no infected nodes found, set to minus infinity
                #print('no infected left')
                signal[q] = -self._inf
            else:
                (n, d) = sp

                # update signal and boundary
                if d != signal[q]:
                    if d < signal[q]:
                        raise ValueError('Signal at {q} got smaller {before} {after}???'.format(q=q, after=d, before=signal[q]))
                    signal[q] = -d
                    #print('update in coboundary', q, signal[q], -d)
                self._boundary[q] = n
                self._coboundary_R[n].add(q)
        del self._coboundary_R[s]

        # for n in g.nodes():
        #     print(f'node {n}:')
        #     if n in susceptibles:
        #         print('   S')
        #     if n in infecteds:
        #         print('   I')
        #     if n in removeds:
        #         print('   R')
        #     print('   signal ', signal[n])
        #     print('   boundary ', boundary[n] if n in boundary else '')
        #     print('   coboundary_S', coboundary_S[n] if n in coboundary_S else '')
        #     print('   coboundary_R', coboundary_R[n] if n in coboundary_R else '')
