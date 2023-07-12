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
"""
Python functions for psi6 computations.
"""
import sys
import numpy as np
from scipy.spatial import Delaunay
from typing import List, Tuple


def extend(config: np.ndarray, box: List[float], n: int) -> np.ndarray:
    """
    Add periodic images for Voronoi construction, trimmed to distance "edge".
    The origin of the coordinates is chosen as the lower-left corner of the simulation box.

    Parameters
    ----------
    config : np.ndarray
        Particle configuration.
    box : List[float]
        Linear size of the box.
    n : int
        Number of particles.

    Returns
    -------
    np.ndarray
        Particles with a duplicated shell of width = edge, used for local-orientational-order calculation.
    """
    edge = 2 * max(box) / np.sqrt(n)
    per = np.concatenate((config, config + [box[0], 0], config + [0, box[1]], config + box, config + [-box[0], 0],
                          config + [0, -box[1]], config - box, config + [box[0], -box[1]], config + [-box[0], box[1]]))
    inx = (per[:, 1] > -edge) & (per[:, 0] > -edge) & (per[:, 1]
                                                       < box[0] + edge) & (per[:, 0] < box[1] + edge)
    per = per[inx, :]
    return per


def hexatic(per: np.ndarray) -> Tuple[np.ndarray, np.ndarray, Delaunay]:
    """
    Find local psi6 using Delaunay tessellation.

    Parameters
    ----------
    per : np.ndarray
        Extended set of particles. (to avoid edge corrections)
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, Delaunay]
        A tuple containing the local psi6 and colormap index of each particle, and the Delaunay tessellation.
    """
    npart = per.shape[0]
    cx = np.zeros(npart, dtype='complex128')  # local orientational order
    counter = np.zeros(npart)  # number of links to site
    tr = Delaunay(per)

    p0 = tr.simplices[:, 0]  # 2N triangles, 3 links per triangle
    p1 = tr.simplices[:, 1]
    p2 = tr.simplices[:, 2]

    d01 = per[p0, :] - per[p1, :]  # 2 N vectors * 3, each link is twice here
    d20 = per[p2, :] - per[p0, :]
    d12 = per[p1, :] - per[p2, :]

    theta01 = 6 * np.arctan2(d01[:, 1], d01[:, 0])  # 2 N angles
    theta20 = 6 * np.arctan2(d20[:, 1], d20[:, 0])  #
    theta12 = 6 * np.arctan2(d12[:, 1], d12[:, 0])  # making 6N contributions

    # phase on the three sides of each triangle
    expt01 = np.exp(complex(0, 1) * theta01)
    expt20 = np.exp(complex(0, 1) * theta20)
    expt12 = np.exp(complex(0, 1) * theta12)

    np.add.at(cx, p0, expt01)  # distribute phase from links to 2 nodes
    np.add.at(cx, p1, expt01)

    np.add.at(cx, p2, expt20)
    np.add.at(cx, p0, expt20)

    np.add.at(cx, p1, expt12)
    np.add.at(cx, p2, expt12)

    np.add.at(counter, p0, 2)  # counter number of links to node
    np.add.at(counter, p1, 2)
    np.add.at(counter, p2, 2)

    # Normalize local orientational order by the number of neighbors
    cx = np.divide(cx, counter)
    index_mapping = (np.angle(cx) + np.pi) / 2. / np.pi  # colormap index in [0 1]
    return cx, index_mapping, tr


def global_psi6(config: np.ndarray, box: List[float]) -> np.cdouble:
    """
    Calculate the orientational-order parameter of a hard-disk configuration.
    The origin of the coordinates is chosen as the lower-left corner of the simulation box.

    Parameters
    ----------
    config : np.ndarray
        The hard-disk configuration.
    box : List[float]
        The geometry of the simulation box.

    Returns
    -------
    np.cdouble
        The global psi6 of the provided configuration.
    """
    n = config.shape[0]
    extended = extend(config, box, n)
    cx, _, _ = hexatic(extended)
    return np.average(cx[:n])


if __name__ == "__main__":
    # Calculate the psi6 of a configuration in a square box.
    # The origin of coordinates is set to be the lower left corner of the simulation box.
    # Example usage: "python3 psi6.py configuration.dat"
    box = [1, 1]
    configuration = np.loadtxt(sys.argv[1])
    print(global_psi6(configuration, box))
