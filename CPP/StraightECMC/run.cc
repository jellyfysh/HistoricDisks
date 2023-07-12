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

void StraightEcmc::run() {
  // The straight ECMC simulation runs here.
  // The whole run is divided into intervals, at the end of which the
  // configurations are sampled.
  // Each interval is further divided into multiple chains.
  // The direction is changed at the end of each chain, and the active disk is
  // resampled at the end of each chain.

  // Disk index and cell index of the first collision
  int k_first_event, cell_first_event;

  // Counts the number of output configuration
  int configuration_counter = 0;

  // Counts number of chains, reset at the end of each interval.
  long int chain_counter = 0;

  // Accumulated $\Delta x_{ij}$ during each interval when moving in x
  // direction. Used for pressure computation.
  double accumulated_delta_x_ij = 0;

  // Accumulated $\Delta y_{ij}$ during each interval when moving in y
  // direction. Used for pressure computation.
  double accumulated_delta_y_ij = 0;

  // Accumulated chain length during an interval when moving in x direction.
  double accumulated_length_x = 0;

  // Accumulated chain length during an interval when moving in y direction.
  double accumulated_length_y = 0;

  // Dummy variable to store $\Delta x_{ij} (\Delta y_{ij})$ of the
  // collision.
  double delta_ij;

  // Cutoff for searching events due to limited cell size.
  double displacement_max = min(min(Param::box[0], Param::box[1]) / 2,
								min(Param::cell_size[0], Param::cell_size[1]))
	  - 2 * Param::sigma;

  // Total chain length (time) of the run.
  double total_chain_length =
	  Param::lambda_0 * Param::factor * Param::extra_factor;
  double total_distance_to_go = total_chain_length;

  // Length of the chain.
  double chain_length = 3.125 * sqrt(Param::number_disks) * Param::lambda_0;
  // Number of chains in an interval.
  long int reset = (int)(total_chain_length / Param::n_samples / chain_length);

  double distance_to_go;
  short direction = 0;

  double V_relative =
	  Param::box[0] * Param::box[1] / 2 / sqrt(3) / double(Param::number_disks)
		  / Param::sigma / Param::sigma;
  cout << "V/V_0: " << V_relative << endl;
  cout << "total chain length " << total_chain_length << endl;
  cout << "chain length " << chain_length << endl;
  time_start = clock();
  k_first_event = 0;
  cell_first_event = 0;
  while (true) {
	direction = 1 - direction;
	distance_to_go = min(chain_length, total_distance_to_go);
	choose_active(cell_cur, k_cur);
	if (direction == 0) {
	  accumulated_length_x += distance_to_go;
	  while (true) {   //loop over collisions
		double displacement_first_event = distance_to_go;
		expl_cell_x(displacement_first_event,
					k_first_event,
					cell_first_event,
					delta_ij);

		displacement_first_event = max(displacement_first_event, .0);
		double new_position = cell[cell_cur][k_cur][0]
			+ min(min(displacement_first_event, distance_to_go),
				  displacement_max);

		refresh_cell_x(new_position);

		if (displacement_max < min(displacement_first_event, distance_to_go)) {
		  // Displacement limited by cell size.
		  // The displacement continues with the same disk and the same
		  // direction.
		  distance_to_go -= displacement_max;
		  continue;
		} else if (displacement_first_event < distance_to_go) {
		  // A regular collision
		  accumulated_delta_x_ij += delta_ij;
		  distance_to_go = distance_to_go - displacement_first_event;
		  k_cur = k_first_event;
		  cell_cur = cell_first_event;
		  number_collisions++;
		  continue;
		} else {
		  // End of chain reached before any collision.
		  break;
		}
	  }
	} else {
	  accumulated_length_y += distance_to_go;
	  while (true) {   //loop over collisions
		double displacement_first_event = distance_to_go;
		expl_cell_y(displacement_first_event,
					k_first_event,
					cell_first_event,
					delta_ij);

		displacement_first_event = max(displacement_first_event, .0);
		double new_position = cell[cell_cur][k_cur][1]
			+ min(min(displacement_first_event, distance_to_go),
				  displacement_max);

		refresh_cell_y(new_position);

		if (displacement_max < min(displacement_first_event, distance_to_go)) {
		  distance_to_go -= displacement_max;
		  continue;
		} else if (displacement_first_event < distance_to_go) {
		  accumulated_delta_y_ij += delta_ij;
		  distance_to_go = distance_to_go - displacement_first_event;
		  k_cur = k_first_event;
		  cell_cur = cell_first_event;
		  number_collisions++;
		  continue;
		} else {
		  break;
		}
	  }
	}
	chain_counter++;
	if (chain_counter >= reset) {
	  // accumulated_length_x = accumulated_length_y,
	  // and pressure = (pressure_x + pressure_y) / 2
	  double pressure = ((accumulated_delta_x_ij + accumulated_delta_y_ij)
		  / (accumulated_length_x + accumulated_length_y) + 1) / V_relative;
	  double pressure_x =
		  (accumulated_delta_x_ij / accumulated_length_x + 1) / V_relative;
	  double pressure_y =
		  (accumulated_delta_y_ij / accumulated_length_y + 1) / V_relative;

	  out_series_h5(static_cast<unsigned int>(configuration_counter),
					pressure,
					"/pressure");
	  out_series_h5(static_cast<unsigned int>(configuration_counter),
					pressure_x,
					"/pressure_x");
	  out_series_h5(static_cast<unsigned int>(configuration_counter),
					pressure_y,
					"/pressure_y");
	  out_configuration_h5(configuration_counter);
	  accumulated_length_x = 0;
	  accumulated_length_y = 0;
	  chain_counter = 0;
	  accumulated_delta_x_ij = 0;
	  accumulated_delta_y_ij = 0;
	  configuration_counter++;
	}
	total_distance_to_go -= chain_length;
	if (total_distance_to_go < 0) break;
  }
  time_end = clock();
}

void StraightEcmc::choose_active(int &cell_chosen, int &k_chosen) {
  // Uniformly samples a disk from the cell system.
  int number_ocp = 0;
  uniform_int_distribution<int>
	  cell_distribution(0, Param::total_number_cell - 1);
  uniform_int_distribution<int> k_distribution(0, Param::n_cell_max - 1);
  while (true) {
	cell_chosen = cell_distribution(Param::random_generator);
	number_ocp = cell_ocp[size_t(cell_chosen)];
	k_chosen = k_distribution(Param::random_generator);
	if (k_chosen < number_ocp) break;
  }
}
