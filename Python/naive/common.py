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
"""Module for common functions to simulate hard disks in a periodic box."""
import math
import random
from typing import List, Sequence


def correct_periodic_position(position: Sequence[float], box: Sequence[float]) -> List[float]:
    """
    Return the given position corrected for periodic boundary conditions in the given simulation box.

    Parameters
    ----------
    position : Sequence[float]
        The position vector.
    box : Sequence[float]
        The geometry of the simulation box.

    Returns
    -------
    List[float]
        The position vector after considering periodic boundary condition.
    """
    return [p % b for p, b in zip(position, box)]


def separation_vector(position_one: Sequence[float], position_two: Sequence[float],
                      box: Sequence[float]) -> List[float]:
    """
    Return the shortest separation vector position_one - position_two between the two given positions under
    consideration of periodic boundary conditions of the given simulation box.

    Parameters
    ----------
    position_one : Sequence[float]
        The first position.
    position_two : Sequence[float]
        The second position.
    box : Sequence[float]
        The geometry of the simulation box.

    Returns
    -------
    List[float]
        The shortest separation vector.
    """
    delta = correct_periodic_position([one - two for one, two in zip(position_one, position_two)], box)
    for i in range(len(delta)):
        if delta[i] > box[i] / 2.0:
            delta[i] -= box[i]
    return delta


def create_crystal(n_x: int, n_y: int, sigma: float, box: Sequence[float]) -> List[List[float]]:
    """
    Create an initial crystalline hard-disk configuration in the given simulation box so that the disks are located on
    the triangular lattice of a fully packed configuration.

    Parameters
    ----------
    n_x : int
        The number of disks per row in the lattice.
    n_y : int
        The number of rows in the lattice.
    sigma : float
        The radius of the disks.
    box : Sequence[float]
        The geometry of the box.

    Returns
    -------
    List[List[float]]
        The list of the initial two-dimensional hard-disk positions.

    Raises
    ------
    RuntimeError
        If the n_x * n_y hard disks of radius sigma do not fit in the specified simulation box.
    """
    n = n_x * n_y
    pos = [[0.0, 0.0] for _ in range(n)]
    distance_x = box[0] / n_x
    if distance_x < 2 * sigma:
        raise RuntimeError("The specified number of hard disks do not fit into the given simulation box.")
    distance_y = box[1] / n_y
    for i in range(n_y):
        for j in range(n_x):
            pos[i * n_x + j] = correct_periodic_position(
                [distance_x * j + 0.5 * distance_x * (i % 2), i * distance_y], box)
    return pos
    

def create_packed(n: int, sigma: float, box: Sequence[float]) -> List[List[float]]:
    """
    Create and initial hard-disk configuration in the given simulation box so the disks are located on a triangular
    lattice with edge length 2.05 * sigma.

    Parameters
    ----------
    n : int
        The number of disks.
    sigma : float
        The radius of the disks.
    box : Sequence[float]
        The geometry of the box.

    Returns
    -------
    List[List[float]]
        The list of the initial two-dimensional hard-disk positions.

    Raises
    ------
    RuntimeError
        If the n hard disks of radius sigma do not fit in the specified simulation box.
    """
    pos = [[0.0, 0.0] for _ in range(n)]
    i = 1
    j = 0
    displacement = 2.05 * sigma
    filling_low = True
    sqrt3 = math.sqrt(3)
    while i + j < n:
        previous = pos[i + j - 1]
        if filling_low:
            if previous[0] + displacement + 2.0 * sigma >= box[0]:
                move = [displacement / 2.0 - previous[0], displacement * sqrt3 / 2.0]
                filling_low = False
                i += 1
            else:
                move = [displacement, 0.0]
                i += 1
        else:
            if previous[0] + displacement + 1.0 * sigma >= box[0]:
                move = [-previous[0], displacement * sqrt3 / 2.0]
                filling_low = True
                j += 1
            else:
                move = [displacement, 0.0]
                j += 1
        pos[j + i - 1] = [p + m for p, m in zip(previous, move)]
    if pos[-1][1] >= box[1]:
        raise RuntimeError("The specified number of hard disks do not fit into the given simulation box.")
    return pos


def sample_vel(n: int) -> List[List[float]]:
    """
    Sample n uniformly distributed two-dimensional unit vectors as initial velocities.

    Parameters
    ----------
    n : int
        The number of sampled velocities.

    Returns
    -------
    List[List[float]]
        The list of sampled two-dimensional unit velocities.
    """
    vel = []
    for _ in range(n):
        theta = random.uniform(0.0, 2.0 * math.pi)
        vel.append([math.cos(theta), math.sin(theta)])
    return vel
