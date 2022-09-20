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
# Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, and Werner Krauth,
# Hard-disk computer simulations---a historic perspective,
# arXiv e-prints: 2207.07715 (2022), https://arxiv.org/abs/2207.07715.
#
"""
Executable Python script that samples the positions of hard disks in a periodic box using the Newtonian event-chain
Monte Carlo algorithm. The number of disks, density, and box aspect ratio are set by command-line arguments.

The script computes the minimum over all pair collision times for the active disk and all other disks. The code then
updates the position of the active disk, as well as the velocities of both colliding disks.

The number of samples, the number of chains between samplings, and the chain time can also be set by the command-line
arguments. By default, each chain has a chain time of 80.0, and there is a single chain between two samples. In total
1000 samples are produced by default.

For more information about the command-line arguments, use the -h (or --help) command-line argument of this script.
An exemplary run can be started via
"python3 ECMC_Newtonian.py 2 2 0.28 crystal --chain_time 80.0 --n_chains 1 --n_samples 10".

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
parser.add_argument("-t", "--chain_time", help="length for each chain (default=80.0)", default=80.0, type=float)
parser.add_argument("-c", "--n_chains", help="number of chains between sampling (default=1)", default=1, type=int)
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

# Only the collisions within a box centered at the active disk are considered. This cutoff prevents the active disk from
# interacting with the disks out of the box.
cutoff = min(box[0], box[1]) / 2.0 - 2.0 * sigma


def find_event(pos_active: Sequence[float], pos_target: Sequence[float], vel_active: Sequence[float],
               sigma: float, box: Sequence[float]) -> float:
    """
    Compute the time when the active hard disk with the given velocity collides with the target disk.

    Parameters
    ----------
    pos_active : Sequence[float]
        The position of the active hard disk.
    pos_target : Sequence[float]
        The position of the target hard disk.
    vel_active : Sequence[float]
        The velocity of the active disk.
    sigma : float
        The radius of the hard disks.
    box : Sequence[float]
        The geometry of the box.

    Returns
    -------
    float
        The time of the collision of the disks.
    """
    vel_sq = vel_active[0] ** 2 + vel_active[1] ** 2
    pos_rel = separation_vector(pos_active, pos_target, box)
    dist_sq = pos_rel[0] ** 2 + pos_rel[1] ** 2
    scal = vel_active[0] * pos_rel[0] + vel_active[1] * pos_rel[1]
    upsilon = scal ** 2 - vel_sq * (dist_sq - 4.0 * sigma ** 2)
    if upsilon > 0.0 > scal:
        return -(scal + math.sqrt(upsilon)) / vel_sq
    else:
        return math.inf


for sample in range(args.n_samples * args.n_chains):
    chain_time = args.chain_time
    active = random.randint(0, n - 1)
    vel = sample_vel(n)
    while chain_time > 0.0:
        time_cutoff = min(chain_time, cutoff / math.sqrt(vel[active][0] ** 2 + vel[active][1] ** 2))
        events = [(target, find_event(pos[active], pos[target], vel[active], sigma, box))
                  for target in range(n) if target != active]
        events.append((active, time_cutoff))
        first_event = min(events, key=lambda t: t[1])
        target, event_time = first_event
        for d in range(2):
            pos[active][d] += event_time * vel[active][d]
        pos[active] = correct_periodic_position(pos[active], box)
        chain_time -= event_time
        if active != target:
            sep = separation_vector(pos[target], pos[active], box)
            e_parallel = [c / 2.0 / sigma for c in sep]
            dot = (vel[target][0] - vel[active][0]) * e_parallel[0] + (vel[target][1] - vel[active][1]) * e_parallel[1]
            vel[active][0] += e_parallel[0] * dot
            vel[active][1] += e_parallel[1] * dot
            vel[target][0] -= e_parallel[0] * dot
            vel[target][1] -= e_parallel[1] * dot
            active = target
    if (sample + 1) % args.n_chains == 0:
        print(*iter(comp for s in pos for comp in s))
