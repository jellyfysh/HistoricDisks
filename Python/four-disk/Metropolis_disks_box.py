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
Executable Python script that samples the positions of four hard-disks in a non-periodic square box of side length 1.0
using the Metropolis algorithm.

In each Metropolis trial, a move is uniformly sampled from a square of side length 2 * delta. This move is proposed for
a random hard disk and only accepted if it leads to a legal configuration where no hard disks overlap with each other or
with the wall.

This script samples the positions of all four hard disks after a given number of (proposed) Metropolis moves and prints
them to stdout. The (2 * k)th and (2 * k + 1)th floats in the output are the x- and y-positions of the kth disk,
respectively.

The modifiable parameters are contained in a single code block below.
"""
import random

###### Start of modifiable parameters ######
sigma = 0.15  # The radius of the hard disks.
delta = 0.1  # The range of the proposed Metropolis move.
n_samples = 10000  # The total number of samples that will be taken before this script finishes.
sampling_interval = 5000  # The number of (proposed) Metropolis moves between two samples.
######  End of modifiable parameters  ######

pos = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
sigma_sq = sigma ** 2
for sample in range(n_samples * sampling_interval):
    a = random.choice(pos)
    b = [a[0] + random.uniform(-delta, delta), a[1] + random.uniform(-delta, delta)]
    min_dist = min((b[0] - c[0]) ** 2 + (b[1] - c[1]) ** 2 for c in pos if c != a)
    box_cond = min(b[0], b[1]) < sigma or max(b[0], b[1]) > 1.0 - sigma
    if not (box_cond or min_dist < 4.0 * sigma ** 2):
        a[:] = b
    if (sample + 1) % sampling_interval == 0:
        print(*iter(comp for s in pos for comp in s))
