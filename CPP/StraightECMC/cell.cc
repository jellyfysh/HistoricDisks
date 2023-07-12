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
#include <cmath>

using namespace std;

void StraightEcmc::expl_cell_x(double &L_min,
							   int &k_min,
							   int &cell_min,
							   double &delta_x_ij) {
  // Checks the neighboring cells of the active disk (specified by cell_cur
  // and k_cur) for potential collisions when the active disk moves in the
  // x-direction. The cells to be explored are (5, 1, 7, 4, 2, 8).
  // Updates L_min, the minimal distance to go before the collision;
  // cell_min and k_min, the target of the collision; delta_x_ij, the distance
  // in x-direction when the active disk and the target are in contact.
  // Loop is flattened for better performance.
  double dx[2];
  double x_cur[2] = {cell[cell_cur][k_cur][0], cell[cell_cur][k_cur][1]};
  int cell_act;
  double L_coll;
  int j;
  delta_x_ij = 0;
  cell_act = cell_neighbour[cell_cur][5];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[1] = cell[cell_act][j][1] - x_cur[1];
	if (abs(dx[1]) < 2 * Param::sigma) {
	  dx[0] = cell[cell_act][j][0] - x_cur[0] + Param::cell_size[0];
	  L_coll =
		  dx[0] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
	  if (L_coll < L_min) {
		L_min = L_coll;
		k_min = j;
		cell_min = cell_act;
		delta_x_ij = sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][1];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[1] = x_cur[1] + Param::cell_size[1] - cell[cell_act][j][1];
	if (abs(dx[1]) < 2 * Param::sigma) {
	  dx[0] = cell[cell_act][j][0] - x_cur[0];
	  if (dx[0] > 1.e-16) {
		L_coll =
			dx[0] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		if (L_coll < L_min) {
		  delta_x_ij =
			  sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		  L_min = L_coll;
		  k_min = j;
		  cell_min = cell_act;
		}
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][7];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[1] = cell[cell_act][j][1] - x_cur[1] + Param::cell_size[1];
	if (abs(dx[1]) < 2 * Param::sigma) {
	  dx[0] = cell[cell_act][j][0] - x_cur[0];
	  if (dx[0] > 1.e-16) {
		L_coll =
			dx[0] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		if (L_coll < L_min) {
		  delta_x_ij =
			  sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		  L_min = L_coll;
		  k_min = j;
		  cell_min = cell_act;
		}
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][4];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[1] = cell[cell_act][j][1] - x_cur[1];
	if (abs(dx[1]) < 2 * Param::sigma) {
	  dx[0] = cell[cell_act][j][0] - x_cur[0];
	  if (dx[0] > 1.e-14) {
		L_coll =
			dx[0] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		if (L_coll < L_min) {
		  delta_x_ij =
			  sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		  L_min = L_coll;
		  k_min = j;
		  cell_min = cell_act;
		}
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][2];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[1] = x_cur[1] - cell[cell_act][j][1] + Param::cell_size[1];
	if (abs(dx[1]) < 2 * Param::sigma) {
	  dx[0] = cell[cell_act][j][0] - x_cur[0] + Param::cell_size[0];
	  L_coll =
		  dx[0] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
	  if (L_coll < L_min) {
		delta_x_ij = sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		L_min = L_coll;
		k_min = j;
		cell_min = cell_act;
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][8];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[1] = cell[cell_act][j][1] - x_cur[1] + Param::cell_size[1];
	if (abs(dx[1]) < 2 * Param::sigma) {
	  dx[0] = cell[cell_act][j][0] - x_cur[0] + Param::cell_size[0];
	  L_coll =
		  dx[0] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
	  if (L_coll < L_min) {
		delta_x_ij = sqrt(4 * (Param::sigma * Param::sigma) - (dx[1] * dx[1]));
		L_min = L_coll;
		k_min = j;
		cell_min = cell_act;
	  }
	}
  }
}

void StraightEcmc::expl_cell_y(double &L_min,
							   int &k_min,
							   int &cell_min,
							   double &delta_y_ij) {
  // Checks the neighboring cells of the active disk (specified by cell_cur
  // and k_cur) for potential collisions when the active disk moves in the
  // y-direction. The cells to be explored are (7, 5, 3, 4, 8, 6).
  // Updates L_min, the minimal distance to go before the collision;
  // cell_min and k_min, the target of the collision; delta_y_ij, the distance
  // in y-direction when the active disk and the target are in contact.
  // Loop is flattened for better performance.
  double dx[2];
  double x_cur[2] = {cell[cell_cur][k_cur][0], cell[cell_cur][k_cur][1]};
  int cell_act;
  double L_coll;
  int j;
  delta_y_ij = 0;
  cell_act = cell_neighbour[cell_cur][7];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[0] = cell[cell_act][j][0] - x_cur[0];
	if (abs(dx[0]) < 2 * Param::sigma) {
	  dx[1] = cell[cell_act][j][1] - x_cur[1] + Param::cell_size[1];
	  L_coll =
		  dx[1] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
	  if (L_coll < L_min) {
		delta_y_ij = sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		L_min = L_coll;
		k_min = j;
		cell_min = cell_act;
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][5];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[0] = -x_cur[0] + Param::cell_size[0] + cell[cell_act][j][0];
	if (abs(dx[0]) < 2 * Param::sigma) {
	  dx[1] = cell[cell_act][j][1] - x_cur[1];
	  if (dx[1] > 1.e-16) {
		L_coll =
			dx[1] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		if (L_coll < L_min) {
		  delta_y_ij =
			  sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		  L_min = L_coll;
		  k_min = j;
		  cell_min = cell_act;
		}
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][3];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[0] = cell[cell_act][j][0] - x_cur[0] - Param::cell_size[0];
	if (abs(dx[0]) < 2 * Param::sigma) {
	  dx[1] = cell[cell_act][j][1] - x_cur[1];
	  if (dx[1] > 1.e-16) {
		L_coll =
			dx[1] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		if (L_coll < L_min) {
		  delta_y_ij =
			  sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		  L_min = L_coll;
		  k_min = j;
		  cell_min = cell_act;
		}
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][4];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[0] = cell[cell_act][j][0] - x_cur[0];
	if (abs(dx[0]) < 2 * Param::sigma) {
	  dx[1] = cell[cell_act][j][1] - x_cur[1];
	  if (dx[1] > 1.e-14) {
		L_coll =
			dx[1] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		if (L_coll < L_min) {
		  delta_y_ij =
			  sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		  L_min = L_coll;
		  k_min = j;
		  cell_min = cell_act;
		}
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][8];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[0] = -x_cur[0] + cell[cell_act][j][0] + Param::cell_size[0];
	if (abs(dx[0]) < 2 * Param::sigma) {
	  dx[1] = cell[cell_act][j][1] - x_cur[1] + Param::cell_size[1];
	  L_coll =
		  dx[1] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
	  if (L_coll < L_min) {
		delta_y_ij = sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		L_min = L_coll;
		k_min = j;
		cell_min = cell_act;
	  }
	}
  }
  cell_act = cell_neighbour[cell_cur][6];
  for (j = 0; j < cell_ocp[size_t(cell_act)]; j++) {
	dx[0] = cell[cell_act][j][0] - x_cur[0] - Param::cell_size[0];
	if (abs(dx[0]) < 2 * Param::sigma) {
	  dx[1] = cell[cell_act][j][1] - x_cur[1] + Param::cell_size[1];
	  L_coll =
		  dx[1] - sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
	  if (L_coll < L_min) {
		delta_y_ij = sqrt(4 * (Param::sigma * Param::sigma) - (dx[0] * dx[0]));
		L_min = L_coll;
		k_min = j;
		cell_min = cell_act;
	  }
	}
  }
}
