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
Executable Python script that samples the positions of hard disks in a periodic box using the straight event-chain Monte
Carlo algorithm. The number of disks, density, and box aspect ratio are set by command-line arguments.

The script computes the minimum over all pair collision times for the active disk and all other disks. The code then
updates the position of the active disk, and transfers its velocity to the other colliding disk. The velocity of the
active disk is restricted to (1, 0) and (0, 1), resulting in an easy implementation of periodic boundary conditions.

The number of samples, the number of chains between samplings, and the chain time can also be set by the command-line
arguments. By default, each chain has a chain time of 0.24, and there are 1000 chains between two samples. In total 1000
samples are produced by default.

For more information about the command-line arguments, use the -h (or --help) command-line argument of this script.
An exemplary run can be started via
"python3 ECMC_straight.py 2 2 0.28 crystal --chain_time 0.24 --n_chains 1000 --n_samples 10".

This script samples the positions of all hard disks in a given time interval and prints them to stdout. The
(2 * k)th and (2 * k + 1)th floats in the output are the x- and y-positions of the kth disk, respectively. The pressure
in x and in y direction, computed by Eq. 20, can also be printed to stdout. The output of the pressure is done at line
149 and 151, which are commented out by default.
"""
import argparse
import math
import random
from typing import Sequence
from common import create_packed, create_crystal
random.seed(1)

parser = argparse.ArgumentParser()
parser.add_argument("n_x", help="number of disks per row", type=int)
parser.add_argument("n_y", help="number of rows", type=int)
parser.add_argument("eta", help="packing fraction", type=float)
parser.add_argument("shape", choices=["square", "rectangle", "crystal"],
                    help="the shape of the box: square for aspect ratio 1, rectangle for aspect ratio sqrt(3)/2, "
                         "and crystal for aspect ratio compatible with a triangular lattice specified by n_x and n_y",
                    type=str)
parser.add_argument("-t", "--chain_time", help="length for each chain (default=0.24)", default=0.24, type=float)
parser.add_argument("-c", "--n_chains", help="number of chains between sampling (default=1000)", default=1000, type=int)
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

sum_delta_x = [0.0, 0.0]
sum_chain_time = [0.0, 0.0]
direction = random.randint(0, 1)


def find_event(pos_active: Sequence[float], pos_target: Sequence[float], direction: int,
               sigma: float, box: Sequence[float]) -> (float, float):
    """
    Compute the time when the active hard disk with a unit velocity in the given direction collides with the target
    disk. Also, return the distance between the two disks at the collision.

    The distance at the collision is returned as zero of the disks never collide, i.e., if the returned collision time
    is infinite.

    Parameters
    ----------
    pos_active : Sequence[float]
        The position of the active hard disk.
    pos_target : Sequence[float]
        The position of the target hard disk.
    direction : int
        The direction of the unit velocity of the active disk (0 and 1 correspond to a velocities parallel the x- and
        y-axes, respectively).
    sigma : float
        The radius of the hard disks.
    box : Sequence[float]
        The geometry of the box.

    Returns
    -------
    (float, float)
        (The time of the collision of the disks, the distance of the two disks at the collision.)
    """
    distance_perp = abs(pos_target[1 - direction] - pos_active[1 - direction])
    distance_perp = min(distance_perp, box[1 - direction] - distance_perp)
    if distance_perp >= 2.0 * sigma:
        return math.inf, 0.0
    else:
        distance_para = pos_target[direction] - pos_active[direction]
        if distance_para < 0.0:
            distance_para += box[direction]
        elif distance_para == 0.0:
            return math.inf, 0.0
        delta_x = math.sqrt(4.0 * sigma ** 2 - distance_perp ** 2)
        time_of_flight = distance_para - delta_x
        return time_of_flight, delta_x


for sample in range(args.n_samples * args.n_chains):
    active = random.randint(0, n - 1)
    chain_time = args.chain_time
    sum_chain_time[direction] += chain_time
    while chain_time > 0.0:
        events = [(target, *find_event(pos[active], pos[target], direction, sigma, box))
                  for target in range(n) if target != active]
        events.append((active, chain_time, 0.0))
        first_event = min(events, key=lambda t: t[1])
        target, event_time, delta_x = first_event
        # The event time could be slightly negative due to the rounding error of the trigonometry calculation.
        # If the event time is negative, it is set to 0.0 in order to prevent the active disk moving backwards.
        pos[active][direction] += max(event_time, 0.0)
        while pos[active][direction] > box[direction]:
            pos[active][direction] -= box[direction]
        sum_delta_x[direction] += delta_x
        active = target
        chain_time -= event_time
    if (sample + 1) % args.n_chains == 0:
        # P_x calculated using Eq. 20
        # print(n * (1 + sum_delta_x[0] / sum_chain_time[0]))
        # P_y calculated using Eq. 20
        # print(n * (1 + sum_delta_x[1] / sum_chain_time[1]))
        sum_delta_x = [0, 0]
        sum_chain_time = [0, 0]
        print(*iter(comp for s in pos for comp in s))
    direction = random.randint(0, 1)
