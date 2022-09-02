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
using straight event-chain Monte Carlo.

The script computes the next event as the collision of the active disk with another disk or the wall. In a collision
with a wall, the velocity is updated according to a Newtonian collision and the active index remains unchanged. In a
disk collision, the velocity and the active index is transferred to the collision disk. The active disk and its velocity
are resampled after each event chain. Here, the length of an event chain is sampled from an exponential distribution.
The velocities are always parallel to the coordinate axes.

This script samples the positions of all four hard disks in a given time interval and prints them to stdout. The
(2 * k)th and (2 * k + 1)th floats in the output are the x- and y-positions of the kth disk, respectively.

This script can also output the pressure calculated between two samples. It can use the estimators in Eqs (14)
and (20) in [Li2022]. The square box allows to calculate P_x and P_y together. The pressure output in the lines
88 and 90 is commented out by default.

The modifiable parameters are contained in a single code block below.
"""
import random
from common import wall_time, pair_time

###### Start of modifiable parameters ######
sigma = 0.15  # The radius of the hard disks.
n_samples = 10000  # The total number of samples that will be taken before this script finishes.
chain_length = 0.24  # The average time between two resamplings.
sample_chain = 500  # The number of chains between two samples.
######  End of modifiable parameters  ######

pos = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
n = len(pos)
velocities = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]
pairs = [tuple([j for j in range(n) if not j == i]) for i in range(n)]

sum_delta_x = 0.0
sum_t_sim = 0.0
wall_collision_count = 0

for sample in range(n_samples * sample_chain):
    t_current_chain = 0.0
    active = random.randint(0, n - 1)
    vel = random.choice(velocities)
    chain_time = random.expovariate(chain_length)
    sum_t_sim += chain_time
    while True:
        wall_times = [wall_time(pos[active][i], vel[i], sigma) for i in range(2)]
        pair_times = [pair_time(pos[active], vel, pos[i], [0, 0], sigma) for i in pairs[active]]
        next_event = min(wall_times + pair_times)
        if t_current_chain + next_event > chain_time:
            for i in range(2):
                pos[active][i] += vel[i] * (chain_time - t_current_chain)
            break
        t_current_chain += next_event
        for i in range(2):
            pos[active][i] += vel[i] * next_event
        if min(wall_times) < min(pair_times):
            wall_collision_count += 1
            vel[wall_times.index(next_event)] *= -1.0
        else:
            target = pairs[active][pair_times.index(next_event)]
            sum_delta_x += vel[0] * (pos[target][0] - pos[active][0]) + vel[1] * (pos[target][1] - pos[active][1])
            active = pairs[active][pair_times.index(next_event)]
    if (sample + 1) % sample_chain == 0:
        print(*iter(comp for s in pos for comp in s))
        # Pressure as (P_x + P_y) / 2 calculated using 14.
        # print(n * wall_collision_count / sum_t_sim)
        # Pressure as (P_x + P_y) / 2 calculated using 20.
        # print(n + n * (2 * sigma * wall_collision_count + sum_delta_x) / sum_t_sim)
        sum_delta_x = 0.0
        sum_t_sim = 0.0
        wall_collision_count = 0
