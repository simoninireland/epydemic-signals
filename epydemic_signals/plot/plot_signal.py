# Plotting routine to show signal
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

from typing import Dict, Any, List, Tuple
from networkx import spring_layout
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, Normalize, TwoSlopeNorm
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.axes import Axes
from epydemic_signals import Signal


def plot_signal(s: Signal, t: float,
                ax: Axes = None,
                cmap: Colormap = None,
                vmin: float = None, vmax: float = None,
                norm: Normalize = None,
                pos: Dict[Any, Tuple[float, float]] = None,
                title: str = None,
                fontsize: int = 5,
                tickfontsize: int = 3,
                marker: str = '.',
                markersize: float = 0.75):
    '''Draw a colour-coded diagram of a signal.

    :param s: the signal
    :param t: the simulation time
    :param ax: (optional) axes to draw into
    :param cmap: (optional) mapping from signal values to colours
    :param vmin: (optional) minimum signal value
    :param vmax (optional) maximum signal value
    :param norm: (optional) normaliser of signal values
    :param pos: (optional) a mapping of nodes to positions
    :param title: (optional) title for plot
    :param fontsize: (optional) size of label font
    :param tickfontsize: (optional) size of tick font
    :param marker: (optional) marker style for nodes
    :param markersize: (optional) marker size
    '''
    g = s.network()
    s_t = s[t]
    colours = [s_t[n] for n in g.nodes()]
    if vmin is None and vmax is None:
        vs = list(s.values())
        vs.sort()
        (vmin, vmax) = (vs[1], vs[-2])  # avoid endpoints in case of infinities
    elif vmin is None or vmax is None:
        raise ValueError('Need to provide both minimum and maximum signal value, or neither')

    # fil in defaults
    if ax is None:
        # default to global main axes
        ax = plt.gca()
    if pos is None:
        # default to spring layout
        pos = spring_layout(g)
    if cmap is None:
        cmap = 'viridis'
    if isinstance(cmap, str):
        cmap = get_cmap(cmap)
    if norm is None:
        # normalise the colourmap of the signal to be centred on 0
        norm = TwoSlopeNorm(0, vmin, vmax)
    if title is None:
        title = f'Signal ($t = {t:.2f}$)'

    # draw the network coloured by signal value
    ax.set_title(title, fontsize=fontsize)
    ax.scatter(x=[p[0] for p in pos.values()], y=[p[1] for p in pos.values()],
               marker=marker, s=markersize,
               c=colours,
               cmap=cmap, norm=norm)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    # draw the colour bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='2%', pad=0.05)
    m = ScalarMappable(norm=norm, cmap=cmap)
    m.set_array(colours)
    cbar = plt.colorbar(m, cax=cax, ticks=[vmin, 0, vmax], format='%.0f')
    #cbar.set_ticks()
    cax.tick_params(labelsize=tickfontsize, width=0.5, length=3, pad=0.3)
