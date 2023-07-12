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
Executable Python script that samples the positions of hard disks in a periodic box using event-driven molecular
dynamics. The number of disks, density, and box aspect ratio are set by command-line arguments.

The script computes the minimum over all pair collision times for all pairs of disks. The code then updates all the
positions and the velocities of the colliding disks.

The number of samples and the time between two samples can also be set by the command-line arguments. By default, the
interval between two samples are 15.0, and 1000 samples are produced.

For more information about the command-line arguments, use the -h (or --help) command-line argument of this script.
An exemplary run can be started via "python3 molecular_dynamics.py 2 2 0.28 crystal --sample_time 15.0 --n_samples 10".

This script samples the positions of all hard disks in a given time interval and prints them to stdout. The
(2 * k)th and (2 * k + 1)th floats in the output are the x- and y-positions of the kth disk, respectively.
"""
import argparse
import math
import random
from typing import Sequence
from common import correct_periodic_position, separation_vector, create_packed, create_crystal, sample_vel
random.seed(1)

parser = argparse.ArgumentParser()
parser.add_argument("n_x", help="number of disks per row", type=int)
parser.add_argument("n_y", help="number of rows", type=int)
parser.add_argument("eta", help="packing fraction", type=float)
parser.add_argument("shape", choices=["square", "rectangle", "crystal"],
                    help="the shape of the box: square for aspect ratio 1, rectangle for aspect ratio sqrt(3)/2, "
                         "and crystal for aspect ratio compatible with a triangular lattice specified by n_x and n_y",
                    type=str)
parser.add_argument("-t", "--sample_time", help="time between two samples (default=15.0)", default=15.0, type=float)
parser.add_argument("-n", "--n_samples", help="number of samples (default=1000)", default=1000, type=int)
args = parser.parse_args()

n = args.n_x * args.n_y
sigma = math.sqrt(args.eta / (n * math.pi))

if args.shape == "square":
    aspect_ratio = 1.0
    box = (1.0 / math.sqrt(aspect_ratio), math.sqrt(aspect_ratio))
    pos = create_packed(n, sigma, box)
elif args.shape == "rectangle":
    aspect_ratio = math.sqrt(3.0) / 2.0
    box = (1.0 / math.sqrt(aspect_ratio), math.sqrt(aspect_ratio))
    pos = create_packed(n, sigma, box)
else:
    assert args.shape == "crystal"
    aspect_ratio = math.sqrt(3.0) / 2.0 * args.n_y / args.n_x
    box = (1.0 / math.sqrt(aspect_ratio), math.sqrt(aspect_ratio))
    pos = create_crystal(args.n_x, args.n_y, sigma, box)

vel = sample_vel(n)
mean_vel = [0.0, 0.0]
for vel_disk in vel:
    mean_vel[0] += vel_disk[0]
    mean_vel[1] += vel_disk[1]
mean_vel = [m / n for m in mean_vel]
for i in range(n):
    vel[i][0] -= mean_vel[0]
    vel[i][1] -= mean_vel[1]

cutoff = min(box[0], box[1]) / 2.0 - 2.0 * sigma


def find_event(pos_i: Sequence[float], vel_i: Sequence[float], pos_j: Sequence[float],
               vel_j: Sequence[float], sigma: float, box: Sequence[float]) -> float:
    """
    Compute the time when the two hard disks of radius sigma at the given positions with the given velocities in the
    given simulation box with periodic boundary conditions.

    Parameters
    ----------
    pos_i : Sequence[float]
        The position of the first hard disk.
    pos_j : Sequence[float]
        The position of the second hard disk.
    vel_i : Sequence[float]
        The velocity of the first hard disk.
    vel_j : Sequence[float]
        The velocity of the second hard disk.
    sigma : float
        The radius of the hard disks.
    box : Sequence[float]
        The geometry of the box.

    Returns
    -------
    float
        The time of the collision of the two disks or None if the disks never collide.
    """
    vel_rel = [vi - vj for vi, vj in zip(vel_i, vel_j)]
    vel_rel_sq = vel_rel[0] ** 2 + vel_rel[1] ** 2
    if vel_rel_sq == 0.0:
        return math.inf
    pos_rel = separation_vector(pos_i, pos_j, box)
    pos_rel_sq = pos_rel[0] ** 2 + pos_rel[1] ** 2
    scal = vel_rel[0] * pos_rel[0] + vel_rel[1] * pos_rel[1]
    upsilon = scal ** 2 - vel_rel_sq * (pos_rel_sq - 4.0 * sigma ** 2)
    if upsilon > 0.0 > scal:
        return -(scal + math.sqrt(upsilon)) / vel_rel_sq
    else:
        return math.inf


sample_count = args.n_samples
time_to_sample = args.sample_time
vel_max = max(math.sqrt(v[0] ** 2 + v[1] ** 2) for v in vel)

while sample_count > 0:
    events = [(i, j, find_event(pos[i], vel[i], pos[j], vel[j], sigma, box)) for i in range(n) for j in range(i + 1, n)]
    # As required by periodic boundary conditions, collisions with images of the disks are considered.
    # For each disk, only the collisions happening in a region centered at the disk are considered.
    # The region has the same geometry as the box, and it moves together with the disk.
    # To prevent the collision between disks out of this region, a "horizon" in time is implemented.
    # The horizon is obtained by dividing the minimum possible distance between a disk and another disk out of its
    # region, by twice (an upper bound to) the maximum velocity of the disks. The factor two is due to the possible
    # head-on collision.
    events.append((-1, -1, cutoff / vel_max / 2))
    events.append((-2, -2, time_to_sample))
    first_event = min(events, key=lambda m: m[2])
    i, j, t = first_event
    time_to_sample -= t
    for m in range(len(pos)):
        pos[m] = correct_periodic_position([p + t * v for p, v in zip(pos[m], vel[m])], box)
    if i != j:
        delta_x = separation_vector(pos[j], pos[i], box)
        delta_v = [vj - vi for vj, vi in zip(vel[j], vel[i])]
        delta_x_norm = math.sqrt(delta_x[0] ** 2 + delta_x[1] ** 2)
        direction = [dx / delta_x_norm for dx in delta_x]
        delta_v_dot_direction = delta_v[0] * direction[0] + delta_v[1] * direction[1]
        vel[i] = [vi + d * delta_v_dot_direction for vi, d in zip(vel[i], direction)]
        vel[j] = [vj - d * delta_v_dot_direction for vj, d in zip(vel[j], direction)]
        vel_max = max(vel_max, math.sqrt(vel[i][0] ** 2 + vel[i][1] ** 2),
                      math.sqrt(vel[j][0] ** 2 + vel[j][1] ** 2))
    if i == -2:
        time_to_sample = args.sample_time
        sample_count -= 1
        # If vel_max is only updated in line 154, it is a monotonic increasing function of time.
        # We recalculate vel_max at a resampling to prevent it becoming too large.
        vel_max = max(vel, key=lambda v: math.sqrt(v[0] ** 2 + v[1] ** 2))
        vel_max = math.sqrt(vel_max[0] ** 2 + vel_max[1] ** 2)
        print(*iter(comp for s in pos for comp in s))
