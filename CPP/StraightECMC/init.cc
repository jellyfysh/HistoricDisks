/*
 * HistoricDisks - Synopsis of pressure data, sampling algorithms and pressure
 * estimators for the hard-disk model of statistical physics
 * https://github.com/jellyfysh/HistoricDisks
 * Copyright (C) 2022 Botao Li, Yoshihiko Nishikawa, Philipp Höllmer,
 * Louis Carillo, A. C. Maggs, and Werner Krauth
 *
 * This file is part of HistoricDisks.
 *
 * HistoricDisks is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * HistoricDisks is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * HistoricDisks in the LICENSE file. If not, see
 * <https://www.gnu.org/licenses/>.
 *
 * If you use HistoricDisks in published work, please cite the following
 * reference (see [Li2022] in References.bib):
 * Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs,
 * Werner Krauth;
 * Hard-disk pressure computations—a historic perspective.
 * J. Chem. Phys. 21 December 2022; 157 (23): 234111.
 * https://doi.org/10.1063/5.0126437
 *
 */

#include "StraightEcmc.h"
#include <iostream>

using namespace std;

void StraightEcmc::init_cell() {
// Builds the neighbor list of the cell system, and stores the initial
// configuration in the cell system.
// For each cell, its neighbors are indexed as
// 6 7 8
// 3 4 5
// 0 1 2
// where the cell itself is has the index 4 in its neighbor list.
  int n_1, n_2, i_2, j_2;
  for (int j = 0; j < Param::number_cell[1]; j++) {
	for (int i = 0; i < Param::number_cell[0]; i++) {
	  for (int c_2 = -1; c_2 < 2; c_2++) {
		for (int c_1 = -1; c_1 < 2; c_1++) {
		  n_1 = i + j * Param::number_cell[0];
		  i_2 = (i + c_1 + Param::number_cell[0]) % Param::number_cell[0];
		  j_2 = (j + c_2 + Param::number_cell[1]) % Param::number_cell[1];
		  n_2 = i_2 + j_2 * Param::number_cell[0];
		  cell_neighbour[n_1][c_1 + 1 + (c_2 + 1) * 3] = n_2;
		}
	  }
	}
  }

  // puts the initial condition in the Cells
  int i_cell[2];
  for (unsigned int i = 0; i < Param::number_disks; i++) {
	for (unsigned int j = 0; j < 2; j++) {
	  i_cell[j] = int(floor(
		  (positions[i][j] + Param::box[j] / 2) / Param::cell_size[j]));
	  if (i_cell[j] == Param::number_cell[j]) {
		i_cell[j] = 0;
		positions[i][j] -= Param::box[j];
	  }
	}
	n_1 = i_cell[0] + i_cell[1] * Param::number_cell[0];
	if (i == 0) {
	  cell_cur = n_1;
	  k_cur = cell_ocp[size_t(n_1)];
	}
	for (int j = 0; j < 2; j++) {
	  cell[n_1][cell_ocp[size_t(n_1)]][j] = positions[i][j] + Param::box[j] / 2
		  - ((double)(i_cell[j]) + 0.5) * Param::cell_size[j];
	}
	cell_ocp[size_t(n_1)]++;
  }
}

void StraightEcmc::init_pos() {
// Creates an initial configuration according to the provided command-line
// parameters. Works with square or rectangular box.
// The origin of the coordinate system is at the center of the box.
  if (Param::shape == 1) {
	// The box is a square.
	auto n = static_cast<unsigned int>(int(sqrt(Param::number_disks)));
	// The initial configuration is almost fully packed.
	double dx[2] = {1.00001 * Param::sigma * 2., 0.};
	double dy[2] = {1.00001 * Param::sigma, 1.00001 * Param::sigma * sqrt(3.)};
	for (unsigned int i = 0; i < n; i++) {
	  for (unsigned int j = 0; j < n + 2; j++) {
		if (j * n + i + 1 > Param::number_disks) continue;
		positions[j * n + i][0] = fmod(i * dx[0] + j * dy[0], Param::box[0]);
		positions[j * n + i][1] = fmod(i * dx[1] + j * dy[1], Param::box[1]);
		for (int k = 0; k < 2; k++) {
		  if (positions[j * n + i][k] <= -Param::box[k] / 2)
			positions[j * n + i][k] += Param::box[k];
		  if (positions[j * n + i][k] >= Param::box[k] / 2)
			positions[j * n + i][k] -= Param::box[k];
		}
	  }
	}
  } else if (Param::shape == 2) {
	// The box is a rectangle with aspect ratio 1 : \sqrt{3} / 2
	auto n = static_cast<unsigned int>(int(sqrt(Param::number_disks)));
	// The initial configuration is almost fully packed.
	double dx[2] = {1.00001 * Param::sigma * 2, 0};
	double dy[2] = {1.00001 * Param::sigma, 1.00001 * Param::sigma * sqrt(3)};
	for (unsigned int i = 0; i < n; i++) {
	  for (unsigned int j = 0; j < n + 2; j++) {
		if (j * n + i + 1 > Param::number_disks) continue;
		positions[j * n + i][0] = fmod(i * dx[0] + j * dy[0], Param::box[0]);
		positions[j * n + i][1] = fmod(i * dx[1] + j * dy[1], Param::box[1]);
		for (int k = 0; k < 2; k++) {
		  if (positions[j * n + i][k] <= -Param::box[k] / 2)
			positions[j * n + i][k] += Param::box[k];
		  if (positions[j * n + i][k] >= Param::box[k] / 2)
			positions[j * n + i][k] -= Param::box[k];
		}
	  }
	}
  } else if (Param::shape == 0) {
	// The box allows fully-packed configuration.
	// The initial configuration is a perfect crystal if Param::slant is
	// zero.
	// Depending on Param::slant, the lattice can be slanted.
	unsigned int nx = Param::number_disks_x;
	unsigned int ny = Param::number_disks_y;
	double dx[2] =
		{1.0 * Param::box[0] / nx, Param::slant * Param::box[1] / ny / nx};
	double dy[2] = {0.5 * Param::box[0] / nx, 1.0 * Param::box[1] / ny};
	for (unsigned int i = 0; i < nx; i++) {
	  for (unsigned int j = 0; j < ny; j++) {
		if (j * nx + i + 1 > Param::number_disks) continue;
		positions[j * nx + i][0] =
			fmod(i * dx[0] + (j % 2) * dy[0], Param::box[0]);
		positions[j * nx + i][1] = fmod(i * dx[1] + j * dy[1]
											+ 0.5 * Param::slant * Param::box[1]
												/ ny / nx * (j % 2),
										Param::box[1]);
		for (int k = 0; k < 2; k++) {
		  if (positions[j * nx + i][k] <= -Param::box[k] / 2)
			positions[j * nx + i][k] += Param::box[k];
		  if (positions[j * nx + i][k] >= Param::box[k] / 2)
			positions[j * nx + i][k] -= Param::box[k];
		}
	  }
	}
  }
  cout << "Initial configuration created\n";
}

double StraightEcmc::distance_square(unsigned int i, unsigned int j) {
  // Finds the minimal distance between two disks in a periodic system
  // Not used anywhere other than check_overlap()
  double dist_x = abs(positions[i][0] - positions[j][0]);
  dist_x = min(dist_x, Param::box[0] - dist_x);
  double dist_y = abs(positions[i][1] - positions[j][1]);
  dist_y = min(dist_y, Param::box[1] - dist_y);
  return dist_x * dist_x + dist_y * dist_y;
}

bool StraightEcmc::check_overlap() {
  // Loops over all disk pairs to find overlap.
  double min_dist_square = 4 * Param::sigma * Param::sigma;
  for (unsigned int i = 0; i < Param::number_disks; i++) {
	for (unsigned int j = 0; j < i; j++)
	  if (distance_square(i, j) < min_dist_square + 1e-10) return true;
  }
  return false;
}
