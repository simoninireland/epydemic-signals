# Extract a progress signal from an SIR epidemic
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

from heapq import heappush, heappop
from typing import Dict, Any, List, Tuple
from networkx import Graph, single_source_shortest_path
from epydemic import Node, SIR
from epydemic_signals import Signal


class SIRProgressSignal(Signal):
    '''Create the progress signal from the network and the epidemic data.

    :param g: the network
    :param df: the event stream of (time, event, node) triples
    '''

    def __init__(self, g: Graph, evs: List[Tuple[float, str, Node]]):
        super().__init__(g)
        self._createSignal(g, evs)

    def _createSignal(self, g: Graph, evs: List[Tuple[float, str, Node]]):
        inf = g.order() + 1           # a distance longer than the longest possible path
        susceptibles = set(g.nodes())
        removeds = set()
        infecteds = set()

        # create the initial signal, whiich should be an infection event
        (_, e, n) = evs[0]
        if e != SIR.INFECTED:
            raise ValueError(f'First event (on node {n}) is {e}, not an infection')
        boundary = dict()                               # the nearset infected node to a susceptible
        signal = single_source_shortest_path(g, n)
        for m in signal.keys():
            signal[m] = len(signal[m]) - 1
            boundary[m] = n
        susceptibles.remove(n)
        infecteds.add(n)
        coboundary = dict()                             # all the susceptibles we're closest to
        coboundary[n] = susceptibles.copy()
        minv = min(signal.values())
        maxv = max(signal.values())
        self.setBaseSignal(signal)

        # iterate through the rest of the events
        for i in range(1, len(evs)):
            (t, e, s) = evs[i]
            print(t, e, s)

            diff = dict()
            if e == SIR.INFECTED:
                # infection event at n, set signal to zero
                diff[s] = -signal[s]
                signal[s] = 0

                # iterate all susceptible nodes updating signal as the
                # shortest path length to an infected node
                # This is a modified Dijkstra's algorithm that prunes the
                # search tree if it encounters a node whose proposed distance
                # is greater than the distance it has in the signal already
                distance = []
                for m in susceptibles:
                    heappush(distance, (0 if m == s else inf, m))
                unvisited = susceptibles.copy()
                coboundary[s] = set()
                while len(distance) > 0:
                    (_, n) = heappop(distance)
                    if n in unvisited:
                        #print(f'visit {n}')
                        d = 1 + signal[n]
                        for m in g.neighbors(n):
                            if m in susceptibles and m in unvisited:
                                if d < signal[m]:
                                    # update the signal
                                    diff[m] = d - signal[m]
                                    signal[m] = d
                                    heappush(distance, (d, m))

                                    # update the boundary
                                    coboundary[boundary[m]].remove(m)
                                    boundary[m] = s
                                    if n not in coboundary.keys():
                                        coboundary[n] = set()
                                    coboundary[s].add(m)
                                    print(f'Boundary of {m} now {s}')
                                else:
                                    # prune the tree
                                    #print(f'prune {m}')
                                    unvisited.remove(m)
                        unvisited.remove(n)
                susceptibles.remove(s)
                coboundary[boundary[s]].remove(s)
                del boundary[s]
                infecteds.add(s)
                print(f'Coboundary of {s} now', coboundary[s])

            elif e == SIR.REMOVED:
                # removal event, re-compute all distances affected by our removal
                # This is a breadth-first traverse, through susceptible nodes only,
                # from the nodes in the coboundary of the removed node looking
                # for the closest infected node to establish as their boundary
                infecteds.remove(s)
                removeds.add(s)
                for q in coboundary[s]:
                    distance = [(0, q)]
                    visited = set()
                    while len(distance) > 0:
                        (d, n) = distance.pop(0)
                        if n in infecteds:
                            # found an infected, store new closest node
                            boundary[q] = n
                            coboundary[n].add(q)

                            # update signal at this node if needed
                            if d != signal[q]:
                                if d < signal[q]:
                                    raise ValueError(f'Signal at {d} got smaller???')
                                diff[q] = d - signal[q]
                                signal[q] = d

                            # terminate traverse
                            break
                        for m in g.neighbors(n):
                            if m not in removeds and m not in visited:
                                distance.append((d + 1, m))
                        visited.add(n)
                del coboundary[s]

                # find distance from removed node to boundary
                # This is a breadth-first traverse
                # that stops when it finds an infected node at distance d
                distance = [(0, s)]
                visited = set()
                while len(distance) > 0:
                    (d, n) = distance.pop(0)
                    if n in infecteds:
                        # found an infected,  update signal at newly removed node
                        diff[s] = -(signal[s] + d)
                        signal[s] = -d

                        # stop the traverse
                        break
                    for m in g.neighbors(n):
                        if m not in visited:
                            distance.append((d + 1, m))
                    visited.add(n)

                # update the signal for all other removed nodes
                # Again, this is Dijkstra's algorithm modified to prune branches
                # that can't change the signal. There's an additional difficulty
                # that the removeds sub-graph may be disconnected, in which case
                # we will check components that can't have their values updated
                distance = []
                for m in removeds:
                    heappush(distance, (-d if m == s else inf, m))
                unvisited = removeds.copy()
                while len(distance) > 0:
                    (_, n) = distance.pop(0)
                    if n in unvisited:
                        #print(f'visit {n}')
                        d = signal[n] - 1
                        for m in g.neighbors(n):
                            if m in unvisited:
                                if d > signal[m]:
                                    #print(f'{m} -> {d}')
                                    diff[m] = d - signal[m]
                                    signal[m] = d
                                    heappush(distance, (d, m))
                                else:
                                    # prune the tree
                                    #print(f'prune {m}')
                                    unvisited.remove(m)
                        unvisited.remove(n)
            else:
                raise ValueError(f'Unrecognised event {e}')

            # update diffs
            minv = min(minv, min(signal.values()))
            maxv = max(maxv, max(signal.values()))
            self.addDiff(abs(t), diff)
            #print(diff, signal)

        self.setBounds(minv, maxv)
