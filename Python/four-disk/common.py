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
"""Module for common functions to simulate four hard disks in a non-periodic square box of side length 1.0."""
import math
import random
from typing import List, Sequence


def wall_time(pos_comp: float, vel_comp: float, sigma: float) -> float:
    """
    Compute the time when the hard disk with radius sigma with the given position and velocity components hits the wall
    of the square box of side length 1.0.

    Parameters
    ----------
    pos_comp : float
        The current position component of the hard disk.
    vel_comp : float
        The current velocity component of the hard disk.
    sigma : float
        The radius of the hard disk.

    Returns
    -------
    float
        The time of the collision of the hard disk with the wall.
    """
    if vel_comp > 0.0:
        del_t = (1.0 - sigma - pos_comp) / vel_comp
    elif vel_comp < 0.0:
        del_t = (pos_comp - sigma) / abs(vel_comp)
    else:
        del_t = math.inf
    return del_t


def pair_time(pos_a: Sequence[float], vel_a: Sequence[float], pos_b: Sequence[float], vel_b: Sequence[float],
              sigma: float) -> float:
    """
    Compute the time when the two hard disks of radius sigma at the given positions with the given velocities collide.

    Parameters
    ----------
    pos_a : Sequence[float]
        The position of the first hard disk.
    vel_a : Sequence[float]
        The velocity of the first hard disk.
    pos_b : Sequence[float]
        The position of the second hard disk.
    vel_b : Sequence[float]
        The velocity of the second hard disk.
    sigma : float
        The radius of the hard disks.

    Returns
    -------
    float
        The time of the collision of the two disks.
    """
    del_x = [pos_b[0] - pos_a[0], pos_b[1] - pos_a[1]]
    del_x_sq = del_x[0] ** 2 + del_x[1] ** 2
    del_v = [vel_b[0] - vel_a[0], vel_b[1] - vel_a[1]]
    del_v_sq = del_v[0] ** 2 + del_v[1] ** 2
    scal = del_v[0] * del_x[0] + del_v[1] * del_x[1]
    upsilon = scal ** 2 - del_v_sq * (del_x_sq - 4.0 * sigma ** 2)
    if upsilon > 0.0 > scal:
        del_t = - (scal + math.sqrt(upsilon)) / del_v_sq
    else:
        del_t = float('inf')
    return del_t


def sample_vel(temperature: float = 1.0) -> List[List[float]]:
    """
    Sample the velocity of four hard disks on an 8-dimensional sphere, whose radius is determined by the temperature.

    Parameters
    ----------
    temperature : float
        The temperature is the total energy kinetic energy divided by the number of disks.

    Returns
    -------
    List[List[float]]
        The list of velocities for the four hard disks.
    """
    random_numbers = []
    normalizer = 0
    vel = [[] for _ in range(4)]
    for _ in range(8):
        random_numbers.append(random.gauss(0, 1))
        normalizer += random_numbers[-1] ** 2
    normalizer = math.sqrt(normalizer / 8 / temperature)
    for i in range(4):
        vel[i] = [random_numbers[2 * i] / normalizer,
                  random_numbers[2 * i + 1] / normalizer]
    return vel
