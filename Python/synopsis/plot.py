# HistoricDisks - Synopsis of pressure data, sampling algorithms and pressure estimators for the hard-disk model of
# statistical physics
# https://github.com/jellyfysh/HistoricDisks
# Copyright (C) 2022 Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, and Werner Krauth
#
# This file is part of HistoricDisks.
#
# HistoricDisks is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# HistoricDisks is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with HistoricDisks in the LICENSE file.
# If not, see <https://www.gnu.org/licenses/>.
#
# If you use HistoricDisks in published work, please cite the following reference (see [Li2022] in References.bib):
# Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, Werner Krauth;
# Hard-disk pressure computations—a historic perspective.
# J. Chem. Phys. 21 December 2022; 157 (23): 234111. https://doi.org/10.1063/5.0126437
#
"""Executable Python script that plots data groups"""
from enum import Enum
import math
from typing import List, Optional, Sequence, Union, Tuple
import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
from read import convert_x, convert_y, DataGroup, read, XUnit, YUnit


def to_string(data: Union[Sequence[float], Sequence[int]], precision: int = 4) -> List[str]:
    """
    Convert a sequence of float data to a sequence of strings with the specified precision.

    By default, the precision is four digits after the decimal point.

    Parameters
    ----------
    data : Union[Sequence[float], Sequence[int]]
        The sequence of data.
    precision : int, optional
        The precision of the floats in the returned strings.

    Returns
    -------
    List[str]
        The list of converted data.
    """
    return [f"{point:.{precision}f}" for point in data]


class Zoom(Enum):
    """
    Enumeration class for the possible zoom settings of the interval of the ticks in the plot.
    """
    NONE = "None"  # Default interval of ticks.
    IN = "In"  # Smaller interval of ticks.
    OUT = "Out"  # Larger interval of ticks.


def plot_init(xmin: float = 1.24, xmax: float = 1.36, ymin: float = 6.0, ymax: float = 9.0,
              zoomed: Zoom = Zoom.NONE) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Set up the canvas of the plot.

    There are four x-axes and two y-axes. Each axis corresponds to a unit.
    The plotted region is specified by the parameters xmin, xmax, ymin, and ymax.
    The parameter 'zoomed' controls the interval of the ticks on the axes. Possible values are described in the Zoom
    enumeration class.

    Parameters
    ----------
    xmin : float, optional
        The lower limit of x.
    xmax : float, optional
        The upper limit of x.
    ymin : float, optional
        The lower limit of y.
    ymax : float, optional
        The upper limit of y.
    zoomed : Zoom, optional
        The zoom setting for the interval of the ticks.

    Returns
    -------
    matplotlib.figure.Figure
        The matplotlib figure object.
    matplotlib.axes.Axes
        The matplotlib axes object.
    """

    fig, host = plt.subplots(figsize=(6.5, 5.5))

    # Creation of the axes.
    ay1 = host.twinx()
    ax1 = host.twiny()
    ax2 = host.twiny()
    ax3 = host.twiny()

    # Labels of the axes.
    host.set_xlabel(r"$V/V_0$", fontsize=14)
    host.set_ylabel(r"$\beta P V_0 / N$", fontsize=14, labelpad=-10)
    ay1.set_ylabel(r"$\beta P (2\sigma)^2$", fontsize=14, labelpad=-30)
    ax1.set_xlabel(r"$\frac{v}{(2\sigma)^2}$", labelpad=15, fontsize=18)
    ax2.set_xlabel(r"$\eta$", fontsize=14, labelpad=5)
    ax3.set_xlabel(r"$\rho$", fontsize=14, labelpad=0)

    # Put the eta axis away.
    ax3.xaxis.tick_bottom()
    ax3.xaxis.set_label_position('bottom')
    ax1.spines['top'].set_position(('axes', 0))
    ax3.spines['bottom'].set_position(('axes', 1))

    # Ticks.

    # Invisible ticks.
    ay1.xaxis.set_ticks([])
    ax1.yaxis.set_ticks([])
    ax2.yaxis.set_ticks([])
    ax3.yaxis.set_ticks([])

    # Visible ticks.
    xticks_eta = [0.77 - 0.02 * i for i in range(8)]
    if ymax - ymin > 4:
        delta_y = 2
    elif ymax - ymin > 2:
        delta_y = 1
    elif ymax - ymin > 1:
        delta_y = 0.5
    elif ymax - ymin > 0.5:
        delta_y = 0.2
    elif ymax - ymin > 0.2:
        delta_y = 0.1
    else:
        delta_y = 0.05
    yticks = [0 + i * delta_y for i in range(int((ymax + 2 * delta_y) / delta_y))]
    if zoomed is Zoom.NONE:
        xticks = [1.20 + 0.02 * i for i in range(20)]
    elif zoomed is Zoom.IN:
        xticks = [1.24 + 0.01 * i for i in range(20)]
        xticks_eta = [0.77 - 0.01 * i for i in range(16)]
    else:
        assert zoomed is Zoom.OUT
        xticks = [1.16 + 0.04 * i for i in range(20)]

    xticks_grid = [math.pi / (2 * math.sqrt(3)) / x for x in xticks_eta]
    ay1.yaxis.set_ticks(yticks)
    ax1.xaxis.set_ticks(xticks)

    ax2.xaxis.set_ticks(xticks_grid)
    ax3.xaxis.set_ticks(xticks_grid)

    # Minor ticks.
    yminor_ticks = []
    ystep = (yticks[1] - yticks[0]) / 10
    for tick in yticks[:-1]:
        for i in range(1, 10):
            yminor_ticks.append(i * ystep + tick)

    xminor_ticks = []
    xstep = (xticks_eta[0] - xticks_eta[1]) / 10
    for tick in xticks_eta[:-1]:
        for i in range(0, 10):
            xminor_ticks.append(math.pi / (2 * math.sqrt(3)) / (-i * xstep + tick))

    ay1.yaxis.set_ticks(yminor_ticks, minor=True)
    ax1.xaxis.set_ticks(xminor_ticks, minor=True)
    ax2.xaxis.set_ticks(xminor_ticks, minor=True)
    ax3.xaxis.set_ticks(xminor_ticks, minor=True)

    ay1.yaxis.set_ticklabels(to_string(convert_y(yticks, YUnit.LEGACY, YUnit.MODERN)))
    ax1.xaxis.set_ticklabels(to_string(convert_x(xticks, XUnit.SPECIFIC_VOLUME, XUnit.VOLUME_PER_PARTICLE)))
    ax2.xaxis.set_ticklabels([str(tick)[:4] for tick in xticks_eta])
    ax3.xaxis.set_ticklabels([str(4 * tick / math.pi)[:6] for tick in xticks_eta])
    ax1.tick_params(direction='inout')
    ax1.tick_params(which='minor', length=0)
    ax2.tick_params(direction='inout', labelsize=15)
    ax3.tick_params(direction='inout')
    ay1.tick_params(direction='inout')

    # Host ticks.
    host.xaxis.set_ticks(xticks)
    host.yaxis.set_ticks(yticks)
    host.xaxis.set_ticks(xminor_ticks, minor=True)
    host.yaxis.set_ticks(yminor_ticks, minor=True)
    host.tick_params(axis='x', which='minor', length=0)
    host.tick_params(direction='inout', labelsize=15)

    # General aspect of the graph.
    host.grid(which="major")
    host.grid(which="minor", lw=0.5, ls="--")
    fig.tight_layout()

    # Boundaries.
    host.set_xlim(xmin, xmax)
    host.set_ylim(ymin, ymax)
    ay1.set_ylim(ymin, ymax)
    ax1.set_xlim(xmin, xmax)
    ax2.set_xlim(xmin, xmax)
    ax3.set_xlim(xmin, xmax)
    return fig, host


def plot(data_group: DataGroup, host: matplotlib.axes.Axes, marker: str, indices: Optional[List[int]] = None,
         color: Optional[List[str]] = None, markersize: int = 5, zorder: int = 2, alpha: float = 1.0,
         label: Optional[List[str]] = None) -> None:
    """
    Plot the data sets from a data group on an initialized canvas.

    ----------
    data_group : DataGroup
        The data group being plotted in the figure.
    host : matplotlib.axes.Axes
        The object of the main axes of the figure whose axes are in SPECIFIC_VOLUME and HISTORIC.
    marker : str
        The marker type that should be used in the plot for the given data group.
    indices : List[int] or None, optional
        The list of indices of the data to be plotted in the figure. By default (None), all the data is plotted.
    color : List[str] or None, optional
        The list of colors for the data to be plotted. By default (None), the default Python color scheme is used.
    markersize : int, optional
        The value of the marker size, by default 5.
    zorder : int, optional
        The zorder of the plotted data group, by default 2.
    alpha : float, optional
        The alpha of the plotted data group, by default 1.0.
    label : List[str] or None, optional
        The labels for the data to be plotted. By default (None), the format is
        "(data group title), ($\alpha = aspect ratio$)".
    """
    dataset = data_group.data
    title = data_group.title
    if indices is None:
        indices = np.arange(0, len(dataset), 1, dtype=int)
    for i in indices:
        data = dataset[i]
        tag = data.tag
        plot_x = convert_x(data.x_values, tag.x_unit, XUnit.SPECIFIC_VOLUME)
        error_bar = np.array(convert_y(data.y_errors, tag.y_unit, YUnit.LEGACY))
        plot_y = np.array(convert_y(data.y_values, tag.y_unit, YUnit.LEGACY))
        plot_label = ''
        if label is None:
            plot_label += title
            plot_label += ', $\\alpha = ' + tag.aspect_ratio + '$'
        else:
            plot_label = label[i]
        if color is not None:
            plot_color = color[i]
        else:
            plot_color = None
        if tag.type == 'error bar':
            host.errorbar(plot_x, plot_y, yerr=error_bar, label=plot_label, linestyle="None",
                          capsize=7.5, color=plot_color, marker=marker, markersize=markersize, zorder=zorder,
                          alpha=alpha)
        else:
            host.errorbar(plot_x, plot_y, label=plot_label, linestyle="None",
                          color=plot_color, marker=marker, markersize=markersize, zorder=zorder, alpha=alpha)


def main() -> None:
    """
    Plot all provided data in the repository.
    """
    directory = '../../DigitizedData/'
    file_list = ['Metropolis_1953/Metropolis_1953.csv', 'Alder_1962/Alder_1962.csv', 'Zollweg_1992/Zollweg_1992.csv', 
                 'Jaster_1999/Jaster_1999.csv', 'Jaster_2004/Jaster_2004.csv', 'Mak_2006/Mak_2006.csv', 
                 'Bernard_2011/Bernard_2011.csv', 'Engel_2013/Engel_2013.csv']
    file_list_this_work = ['ThisWork.csv']
    fig, host = plot_init()
    for file_name in file_list:
        data = read(directory + file_name)
        plot(data, host, 'x')
    fig.legend()
    plt.show()
    plt.close()
    fig, host = plot_init()
    for file_name in file_list_this_work:
        data = read(directory + file_name)
        plot(data, host, 'x')
    fig.legend()
    plt.show()


if __name__ == "__main__":
    main()
