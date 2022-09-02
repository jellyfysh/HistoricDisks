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
Executable Python script that samples the positions of four hard-disks in a non-periodic square box of side length 1.0
using Newtonian event-chain Monte Carlo.

Every hard disk has a velocity, but only single active disk actually moves. The script computes the next event as the
collision of the active disk with another disk or the wall. In a collision with a wall, the velocity is updated
according to a Newtonian collision and the active index remains unchanged. The velocity update in a disk collision also
uses Newtonian's collision rules by considering both involved velocities, and the active index is transferred to the
collision disk. The active disk and all velocities are resampled after each sampling.

This script samples the positions of all four hard disks in a given time interval and prints them to stdout. The
(2 * k)th and (2 * k + 1)th floats in the output are the x- and y-positions of the kth disk, respectively.

The modifiable parameters are contained in a single code block below.
"""
import math
import random
from common import wall_time, pair_time, sample_vel

###### Start of modifiable parameters ######
sigma = 0.15  # The radius of the hard disks.
n_samples = 10000  # The total number of samples that will be taken before this script finishes.
sample_time = 80.0  # The time between two samples.
######  End of modifiable parameters  ######

pos = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
n = len(pos)
pairs = [tuple([j for j in range(n) if not j == i]) for i in range(n)]

for sample in range(n_samples):
    t_current_chain = 0.0
    active = random.randint(0, n - 1)
    vel = sample_vel()
    while True:
        wall_times = [wall_time(pos[active][i], vel[active][i], sigma) for i in range(2)]
        pair_times = [pair_time(pos[active], vel[active], pos[i], [0, 0], sigma) for i in pairs[active]]
        next_event = min(wall_times + pair_times)
        if t_current_chain + next_event > sample_time:
            for i in range(2):
                pos[active][i] += vel[active][i] * (sample_time - t_current_chain)
            break
        t_current_chain += next_event
        for i in range(2):
            pos[active][i] += vel[active][i] * next_event
        if min(wall_times) < min(pair_times):
            vel[active][wall_times.index(next_event)] *= -1.0
        else:
            target = pairs[active][pair_times.index(next_event)]
            sep = [pos[target][0] - pos[active][0], pos[target][1] - pos[active][1]]
            abs_sep = math.sqrt(sep[0] ** 2 + sep[1] ** 2)
            e_parallel = [c / abs_sep for c in sep]
            dot = (vel[target][0] - vel[active][0]) * e_parallel[0] + (vel[target][1] - vel[active][1]) * e_parallel[1]
            vel[active][0] += e_parallel[0] * dot
            vel[active][1] += e_parallel[1] * dot
            vel[target][0] -= e_parallel[0] * dot
            vel[target][1] -= e_parallel[1] * dot
            active = target
    print(*iter(comp for s in pos for comp in s))
