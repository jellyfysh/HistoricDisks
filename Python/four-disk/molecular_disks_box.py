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
using event-driven molecular dynamics.

The script computes the minimum over all pair collision times for the six pairs of disks, and over the wall collision
times for the four disks in each direction. The code then updates all the positions and the velocities of the colliding
disks.

This script samples the positions of all four hard disks in a given time interval and prints them to stdout. The
(2 * k)th and (2 * k + 1)th floats in the output are the x- and y-positions of the kth disk, respectively.

This script can also output the pressure calculated between two samples. It can use the estimators in Eqs (13c)
and (19a) in [Li2022]. The pressure output in the lines 87, 88 and 90, 91 is commented out by default.

The modifiable parameters are contained in a single code block below.
"""
import math
from common import wall_time, pair_time, sample_vel

###### Start of modifiable parameters ######
sigma = 0.15  # The radius of the hard disks.
n_samples = 1000000  # The total number of samples that will be taken before this script finishes.
sample_time = 15.0  # The time between two samples.
######  End of modifiable parameters  ######

pos = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]  # Initial disk position
# vel = [[0.21, 0.12], [0.71, 0.18], [-0.23, -0.79], [0.78, 0.1177]]
vel = sample_vel()  # Initial velocity
r = math.sqrt(sum(comp ** 2 for v in vel for comp in v))  # Radius of velocity sphere
walls = [(disk, direction) for disk in range(len(pos)) for direction in range(2)]
pairs = [(disk_one, disk_two) for disk_one in range(len(pos)) for disk_two in range(disk_one + 1, len(pos))]
wall_collision_count = 0
pair_collision_count = 0

for sample in range(n_samples):
    time_after_last_sample = 0.0
    while True:
        wall_times = [wall_time(pos[disk][direction], vel[disk][direction], sigma) for disk, direction in walls]
        pair_times = [pair_time(pos[one], vel[one], pos[two], vel[two], sigma) for one, two in pairs]
        next_event = min(wall_times + pair_times)
        if time_after_last_sample + next_event > sample_time:
            for disk, direction in walls:
                pos[disk][direction] += vel[disk][direction] * (sample_time - time_after_last_sample)
            break
        time_after_last_sample += next_event
        for disk, direction in walls:
            pos[disk][direction] += vel[disk][direction] * next_event
        if min(wall_times) < min(pair_times):
            collision_disk, direction = walls[wall_times.index(next_event)]
            vel[collision_disk][direction] *= -1.0
            wall_collision_count += 1
        else:
            a, b = pairs[pair_times.index(next_event)]
            del_x = [pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]]
            abs_x = math.sqrt(del_x[0] ** 2 + del_x[1] ** 2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [vel[b][0] - vel[a][0], vel[b][1] - vel[a][1]]
            scal = del_v[0] * e_perp[0] + del_v[1] * e_perp[1]
            for k in range(2):
                vel[a][k] += e_perp[k] * scal
                vel[b][k] -= e_perp[k] * scal
            pair_collision_count += 1
    print(*iter(comp for s in pos for comp in s))
    # Pressure as (P_x + P_y) / 2 calculated using 13c.
    # print(math.sqrt(math.pi) / r * math.gamma(len(pos) + 0.5) / math.gamma(len(pos))
    #       * wall_collision_count / sample_time / 2)
    # Pressure calculated using 19a.
    # print(len(pos) + sigma * math.sqrt(math.pi) / r * math.gamma(len(pos) + 0.5) / math.gamma(len(pos))
    #       * (wall_collision_count + math.sqrt(2.0) * pair_collision_count) / sample_time)
    wall_collision_count = 0
    pair_collision_count = 0
