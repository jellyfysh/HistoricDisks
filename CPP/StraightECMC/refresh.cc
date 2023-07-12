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

using namespace std;

void
StraightEcmc::refresh_cell_x(double dx) {
// Update the cell system with after-collision position of the active disk.
// The disk is moved in the x direction.
  int cell_new;
  if (dx > Param::cell_size[0] / 2) {
	// The active disk moves into a new cell.
	// x position in the new cell
	dx -= Param::cell_size[0];
	// Index of the new cell
	cell_new = cell_neighbour[cell_cur][5];
	// Refreshes # of disk in the new cell
	cell_ocp[size_t(cell_new)] = cell_ocp[size_t(cell_new)] + 1;
	// x position in the new cell
	cell[cell_new][cell_ocp[size_t(cell_new)] - 1][0] = dx;
	// y position in the new cell
	cell[cell_new][cell_ocp[size_t(cell_new)] - 1][1] = cell[cell_cur][k_cur][1];
	// Swaps the last disk in the old cell to take k_cur's vacant place
	cell[cell_cur][k_cur][0] = cell[cell_cur][cell_ocp[size_t(cell_cur)] - 1][0];
	cell[cell_cur][k_cur][1] = cell[cell_cur][cell_ocp[size_t(cell_cur)] - 1][1];
	// Refreshes # of disks in the old cell
	cell_ocp[size_t(cell_cur)] = cell_ocp[size_t(cell_cur)] - 1;
	cell_cur = cell_new;
	// The active disk is the last disk in the new cell.
	k_cur = cell_ocp[size_t(cell_new)] - 1;
  } else {
	// The active disk stays in cell_cur.
	cell[cell_cur][k_cur][0] = dx;
  }
}

void StraightEcmc::refresh_cell_y(double dy) {
// Update the cell system with after-collision position of the active disk.
// The disk is moved in the y direction.
  int cell_new;
  if (dy > Param::cell_size[1] / 2) {
	// The active disk moves into a new cell.
	// y position in the new cell
	dy -= Param::cell_size[1];
	// Index of the new cell
	cell_new = cell_neighbour[cell_cur][7];
	// Refreshes # of disk in the new cell
	cell_ocp[size_t(cell_new)] = cell_ocp[size_t(cell_new)] + 1;
	// y position in the new cell
	cell[cell_new][cell_ocp[size_t(cell_new)] - 1][1] = dy;
	// x position in the new cell
	cell[cell_new][cell_ocp[size_t(cell_new)] - 1][0] = cell[cell_cur][k_cur][0];
	// Swaps the last disk in the old cell to take k_cur's vacant place
	cell[cell_cur][k_cur][0] = cell[cell_cur][cell_ocp[size_t(cell_cur)] - 1][0];
	cell[cell_cur][k_cur][1] = cell[cell_cur][cell_ocp[size_t(cell_cur)] - 1][1];
	// Refreshes # of disks in the old cell
	cell_ocp[size_t(cell_cur)] = cell_ocp[size_t(cell_cur)] - 1;
	cell_cur = cell_new;
	// The active disk is the last disk in the new cell.
	k_cur = cell_ocp[size_t(cell_new)] - 1;
  } else {
	// The active disk stays in cell_cur.
	cell[cell_cur][k_cur][1] = dy;
  }
}
