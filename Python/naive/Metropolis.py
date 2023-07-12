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
Executable Python script that samples the positions of hard disks in a periodic box using the Metropolis algorithm.
The number of disks, density, and box aspect ratio are set by command-line arguments.

The script proposes the displacement of a uniformly sampled disk. The proposed position is uniformly sampled in a square
region around the disk center. Only if the proposed position does not introduce any overlap, the proposed position is
accepted.

The number of samples and the moves between two samples can also be set by the command-line arguments. By default, the
number of moves between two samples are 1000, and 1000 samples are produced.

For more information about the command-line arguments, use the -h (or --help) command-line argument of this script.
An exemplary run can be started via
"python3 Metropolis.py 2 2 0.28 crystal --sample_move 1000 --n_samples 10".

This script samples the positions of all hard disks in a given time interval and prints them to stdout. The
(2 * k)th and (2 * k + 1)th floats in the output are the x- and y-positions of the kth disk, respectively.
"""
import argparse
import math
import random
from typing import Sequence
from common import correct_periodic_position, create_packed, create_crystal, separation_vector
random.seed(1)

parser = argparse.ArgumentParser()
parser.add_argument("n_x", help="number of disks per row", type=int)
parser.add_argument("n_y", help="number of rows", type=int)
parser.add_argument("eta", help="packing fraction", type=float)
parser.add_argument("shape", choices=["square", "rectangle", "crystal"],
                    help="the shape of the box: square for aspect ratio 1, rectangle for aspect ratio sqrt(3)/2, "
                         "and crystal for aspect ratio compatible with a triangular lattice specified by n_x and n_y",
                    type=str)
parser.add_argument("-m", "--sample_move", help="number of moves between two samples (default=1000)", default=1000,
                    type=int)
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

four_sigma_sq = 4.0 * sigma ** 2
delta = (math.sqrt(1.0 / n / math.pi) - sigma) / 2.0
for sample in range(args.n_samples * args.sample_move):
    a = random.randint(0, n - 1)
    b = correct_periodic_position([(pos[a][d] + random.uniform(-delta, delta)) % box[d] for d in range(2)], box)
    reject = False
    for c in range(n):
        if c != a:
            separation_vec = separation_vector(b, pos[c], box)
            if separation_vec[0] ** 2 + separation_vec[1] ** 2 < four_sigma_sq:
                reject = True
                break
    if not reject:
        pos[a][:] = b
    if (sample + 1) % args.sample_move == 0:
        print(*iter(comp for s in pos for comp in s))
