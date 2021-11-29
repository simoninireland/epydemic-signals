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
from typing import Dict, Any, List, Tuple, cast
from networkx import Graph, single_source_shortest_path
from epydemic import Node, Edge, SIR, Process, CompartmentedModel
from epydemic_signals import Signal, SignalGenerator


class SIRProgressSignalGenerator(SignalGenerator):
    '''Create the progress signal for an SIR epidemic.

    :param s: the signal
    '''

    def __init__(self, s: Signal):
        super().__init__(s)
        self._inf = None
        self._compartment = dict()
        self._boundary = dict()
        self._coboundary_S = dict()
        self._coboundary_R = dict()

        # register the event handlers
        self.addEventTypeHandler(SIR.INFECTED, self.infect)
        self.addEventTypeHandler(SIR.REMOVED, self.remove)

    def setUp(self):
        '''Capture the initial state of the network as susceptible, infected,
        and removed sets.

        :param g: the network
        :param p: the process'''
        s = self.signal()
        g = s.network()
        p = s.process()
        signal = s[0.0]

        self._inf = g.order() + 1           # a distance longer than the longest possible path
        self._compartment = dict()
        self._compartment[SIR.SUSCEPTIBLE] = set()
        self._compartment[SIR.INFECTED] = set()
        self._compartment[SIR.REMOVED] = set()

        # extract the initial states
        cm = cast(CompartmentedModel, p)
        for n in g.nodes():
            self._compartment[cm.getCompartment(n)].add(n)
            #print(n,  cm.getCompartment(n))

        # signal is initially infinite everywhere
        for n in g.nodes():
            signal[n] = self._inf

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

    def infect(self, t: float, e: Edge):
        '''Adjust the signal for an infection event.

        :param t: the event time
        :param e: the SI edge the infection passed over'''
        (s, _) = e
        g = self.network()
        signal = self.signal()[t]

        # update state
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
        #print('Phase I-1')
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
        #print(f'Sus coboundary of {s} now', coboundary_S[s], 'signal', signal[s], diff[s])
        #print(f'Rem coboundary of {s} now', coboundary_R[s])

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
        # This is a breadth-first traverse, through susceptible nodes only,
        # from the nodes in the coboundary of the removed node looking
        # for the closest infected node to re-establish their boundary
        #print('Phase R-1')
        for q in self._coboundary_S[s]:
            distance = [(0, q)]
            visited = set()
            while len(distance) > 0:
                (d, n) = distance.pop(0)
                if n not in visited:
                    #print(f'visit {n}')
                    visited.add(n)
                    if n in self._compartment[SIR.INFECTED]:
                        # found an infected, store new closest node
                        self._boundary[q] = n
                        self._coboundary_S[n].add(q)

                        # update signal at this node if needed
                        if d != signal[q]:
                            if d < signal[q]:
                                raise ValueError('Signal at {q} got smaller {before} {after}???'.format(q=q, after=d, before=signal[q]))
                            signal[q] = d
                            #print(f'update sus distance {q} {d}')
                        else:
                            #print(f'sus distance {q} unchanged {d}')
                            pass

                        # terminate traverse
                        break
                    for m in g.neighbors(n):
                        if m not in self._compartment[SIR.REMOVED] and m not in visited:
                            distance.append((d + 1, m))
        del self._coboundary_S[s]

        # find distance from removed node to boundary
        # This is a breadth-first traverse
        # that stops when it finds an infected node at distance d
        #print('Phase R-2')
        distance = [(0, s)]
        visited = set()
        while len(distance) > 0:
            (d, n) = distance.pop(0)
            if n not in visited:
                visited.add(n)
                if n in self._compartment[SIR.INFECTED]:
                    # found an infected,  update signal at newly removed node
                    signal[s] = -d

                    # store this as our boundary
                    self._boundary[s] = n
                    self._coboundary_R[n].add(s)
                    #print(f'boundary of {s} {n} signal', signal[s])

                    # stop the traverse
                    break
                for m in g.neighbors(n):
                    if m not in visited:
                        distance.append((d + 1, m))

        # update the signal for all other removed nodes
        # Again, this is a breadth-first traverse modified to prune branches
        # that can't change the signal
        #print('Phase R-3')
        for q in self._coboundary_R[s]:
            distance = [(0, q)]
            visited = set()
            while len(distance) > 0:
                (d, n) = distance.pop(0)
                if n not in visited:
                    #print(f'visit {n}')
                    visited.add(n)
                    if n in self._compartment[SIR.INFECTED]:
                        # found an infected, update signal
                        signal[q] = -d

                        # store this as our boundary
                        self._boundary[q] = n
                        self._coboundary_R[n].add(q)
                        #print(f'boundary of {q} {n} distance -{d}')

                        # stop the traverse
                        break
                    for m in g.neighbors(n):
                        if m not in visited:
                            distance.append((d + 1, m))
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
        #     print('   diff ', diff[n] if n in diff else '')
        #     print('   boundary ', boundary[n] if n in boundary else '')
        #     print('   coboundary_S', coboundary_S[n] if n in coboundary_S else '')
        #     print('   coboundary_R', coboundary_R[n] if n in coboundary_R else '')


    def _createSignal(self, g: Graph, evs: List[Tuple[float, str, Node]]):
        inf = g.order() + 1           # a distance longer than the longest possible path
        susceptibles = set(g.nodes())
        removeds = set()
        infecteds = set()

        # create the initial signal, which should be an infection event
        (_, e, n) = evs[0]
        if e != SIR.INFECTED:
            raise ValueError(f'First event (on node {n}) is {e}, not an infection')
        boundary = dict()                               # the nearest infected node to a non-infected
        signal = single_source_shortest_path(g, n)
        for m in signal.keys():
            signal[m] = len(signal[m]) - 1
            if m != n:
                boundary[m] = n                         # n is the only infected at this stage
        susceptibles.remove(n)
        infecteds.add(n)
        coboundary_S = dict()                           # all the susceptibles an infected is closest to
        coboundary_S[n] = susceptibles.copy()
        coboundary_R = dict()                           # all the removeds an infected is closest to
        coboundary_R[n] = set()
        minv = min(signal.values())
        maxv = max(signal.values())
        self.setBaseSignal(signal)

        # iterate through the rest of the events
        for i in range(1, len(evs)):
            (t, e, s) = evs[i]
            #print(t, e, s)

            diff = dict()
            if e == SIR.INFECTED:
                # infection event at n, update state
                susceptibles.remove(s)
                coboundary_S[boundary[s]].remove(s)
                del boundary[s]
                infecteds.add(s)

                # set signal at s
                diff[s] = -signal[s]
                signal[s] = 0

                # iterate all susceptible and removed nodes updating signal as the
                # shortest path length to an infected node
                # This is a modified Dijkstra's algorithm that prunes the
                # search tree if it encounters a node whose proposed distance
                # is greater than the distance it has in the signal already
                #print('Phase I-1')
                distance = [(0, s)]
                visited = set()
                coboundary_S[s] = set()
                coboundary_R[s] = set()
                while len(distance) > 0:
                    (_, n) = heappop(distance)
                    if n not in visited:
                        #print(f'visit {n}')
                        visited.add(n)
                        d = 1 + signal[n]
                        for m in g.neighbors(n):
                            if m not in visited:
                                if m in susceptibles:
                                    if d < signal[m]:
                                        # update the signal
                                        diff[m] = d - signal[m]
                                        signal[m] = d
                                        heappush(distance, (d, m))
                                        #print(f'propose {m} distance {d}')

                                        # update the boundary
                                        coboundary_S[boundary[m]].remove(m)
                                        boundary[m] = s
                                        coboundary_S[s].add(m)
                                        #print(f'Sus boundary of {m} now {s}')
                                    else:
                                        # prune the tree
                                        #print(f'prune {m}')
                                        #visited.add(m)
                                        pass
                                elif m in removeds:
                                    if -d > signal[m]:
                                        # update the signal
                                        diff[m] = -d - signal[m]
                                        signal[m] = -d
                                        heappush(distance, (d, m))

                                        # update the boundary
                                        coboundary_R[boundary[m]].remove(m)
                                        boundary[m] = s
                                        coboundary_R[s].add(m)
                                        #print(f'Rem boundary of {m} now {s}')
                                    else:
                                        #visited.add(m)
                                        pass
                                else:
                                    # prune the tree
                                    #print(f'prune {m}')
                                    visited.add(m)
                #print(f'Sus coboundary of {s} now', coboundary_S[s], 'signal', signal[s], diff[s])
                #print(f'Rem coboundary of {s} now', coboundary_R[s])

            elif e == SIR.REMOVED:
                # removal event, update state
                infecteds.remove(s)
                removeds.add(s)

                # re-compute all susceptible distances affected by our removal
                # This is a breadth-first traverse, through susceptible nodes only,
                # from the nodes in the coboundary of the removed node looking
                # for the closest infected node to re-establish their boundary
                #print('Phase R-1')
                for q in coboundary_S[s]:
                    distance = [(0, q)]
                    visited = set()
                    while len(distance) > 0:
                        (d, n) = distance.pop(0)
                        if n not in visited:
                            #print(f'visit {n}')
                            visited.add(n)
                            if n in infecteds:
                                # found an infected, store new closest node
                                boundary[q] = n
                                coboundary_S[n].add(q)

                                # update signal at this node if needed
                                if d != signal[q]:
                                    if d < signal[q]:
                                        raise ValueError('Signal at {q} got smaller {before} {after}???'.format(q=q, after=d, before=signal[q]))
                                    diff[q] = d - signal[q]
                                    signal[q] = d
                                    #print(f'update sus distance {q} {d}')
                                else:
                                    #print(f'sus distance {q} unchanged {d}')
                                    pass

                                # terminate traverse
                                break
                            for m in g.neighbors(n):
                                if m not in removeds and m not in visited:
                                    distance.append((d + 1, m))
                del coboundary_S[s]

                # find distance from removed node to boundary
                # This is a breadth-first traverse
                # that stops when it finds an infected node at distance d
                #print('Phase R-2')
                distance = [(0, s)]
                visited = set()
                while len(distance) > 0:
                    (d, n) = distance.pop(0)
                    if n not in visited:
                        visited.add(n)
                        if n in infecteds:
                            # found an infected,  update signal at newly removed node
                            diff[s] = -d - signal[s]
                            signal[s] = -d

                            # store this as our boundary
                            boundary[s] = n
                            coboundary_R[n].add(s)
                            #print(f'boundary of {s} {n} signal', signal[s], diff[s])

                            # stop the traverse
                            break
                        for m in g.neighbors(n):
                            if m not in visited:
                                distance.append((d + 1, m))

                # update the signal for all other removed nodes
                # Again, this is a breadth-first traverse modified to prune branches
                # that can't change the signal
                #print('Phase R-3')
                for q in coboundary_R[s]:
                    distance = [(0, q)]
                    visited = set()
                    while len(distance) > 0:
                        (d, n) = distance.pop(0)
                        if n not in visited:
                            #print(f'visit {n}')
                            visited.add(n)
                            if n in infecteds:
                                # found an infected, update signal
                                diff[q] = -d - signal[q]
                                signal[q] = -d

                                # store this as our boundary
                                boundary[q] = n
                                coboundary_R[n].add(q)
                                #print(f'boundary of {q} {n} distance -{d}')

                                # stop the traverse
                                break
                            for m in g.neighbors(n):
                                if m not in visited:
                                    distance.append((d + 1, m))
                del coboundary_R[s]

            else:
                raise ValueError(f'Unrecognised event {e}')

            # update diffs
            minv = min(minv, min(signal.values()))
            maxv = max(maxv, max(signal.values()))
            self.addDiff(t, diff)
            #print(diff, signal)

            # for n in g.nodes():
            #     print(f'node {n}:')
            #     if n in susceptibles:
            #         print('   S')
            #     if n in infecteds:
            #         print('   I')
            #     if n in removeds:
            #         print('   R')
            #     print('   signal ', signal[n])
            #     print('   diff ', diff[n] if n in diff else '')
            #     print('   boundary ', boundary[n] if n in boundary else '')
            #     print('   coboundary_S', coboundary_S[n] if n in coboundary_S else '')
            #     print('   coboundary_R', coboundary_R[n] if n in coboundary_R else '')

        self.setBounds(minv, maxv)
