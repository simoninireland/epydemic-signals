# Plotting routine to show compartments
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

from typing import Dict, Any, List, Tuple
from networkx import spring_layout
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.axes import Axes
from epydemic_signals import Signal


def plot_compartments(s: Signal, t: float, compartment: List[str],
                      ax: Axes = None,
                      compartment_cmap: Dict[str, Any] = None,
                      pos: Dict[Any, Tuple[float, float]] = None,
                      title: str = None,
                      fontsize: int = 5,
                      marker: str = '.',
                      markersize: float = 0.75):
    '''Draw a colour-coded diagram of the state of an epidemic.

    :param s: the signal
    :param t: the simulation time
    :param compartment: a list of compartments in node order
    :param ax: (optional) axes to draw into
    :param compartment_cmap: (optional) mapping from compartments to colours
    :param pos: (optional) a mapping of nodes to positions
    :param title: (optional) title for plot
    :param fontsize: (optional) size of label font
    :param marker: (optional) marker style for nodes
    :param markersize: (optional) marker size
    '''
    g = s.network()

    # fill in defaults
    if ax is None:
        # default to drraw into the global main axes
        ax = plt.gca()
    if pos is None:
        # default to spring layout
        pos = spring_layout(g)
    if compartment_cmap is None:
        # default to an arbitrary map if none is provided (not very useful...)
        compartment_cmap = dict()
        cmap = get_cmap('tab20')
        i = 1.0 / 40
        for m in compartment:
            if m not in compartment_cmap:
                compartment_cmap[n] = cmap(i)
                i += 1.0 / 20
                if i > 1.0:
                    # roll around for more than 20 compartments (unlikely...)
                    i = 1.0 / 40
    if title is None:
        title = f'Compartments ($t = {t:.2f}$)'

    # draw network coloured by compartment
    ax.set_title(title, fontsize=fontsize)
    ax.scatter(x=[p[0] for p in pos.values()], y=[p[1] for p in pos.values()],
               marker=marker, s=markersize,
               color=[compartment_cmap[m] for m in compartment])
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    # count nodes in each compartment
    nn = dict()
    for m in compartment:
        if m not in nn:
            nn[m] = 1
        else:
            nn[m] += 1

    # draw sidebar divided by fraction per compartment
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='2%', pad=0.05)
    cax.xaxis.set_ticks([])
    cax.yaxis.set_ticks([])
    by = 0.0
    for m in nn.keys():
        ty = by + 1.0 * (nn[m] / len(compartment))
        cax.add_patch(Rectangle((0.0, by), 1.0, ty,
                      fill=True, facecolor=compartment_cmap[m]))
        by=ty
