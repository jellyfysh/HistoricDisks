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
#include <fstream>

using namespace std;

StraightEcmc::StraightEcmc() {
  cout.precision(11);
  number_collisions = 0;
  cell_ocp.resize(size_t(Param::total_number_cell), 0.);
  cell_neighbour.resize(boost::extents[Param::total_number_cell][9]);
  cell.resize(boost::extents[Param::total_number_cell][Param::n_cell_max][2]);
  positions.resize(boost::extents[Param::number_disks][2]);

  cout << "Number of disks = " << Param::number_disks << "\n";
  cout << "Radius = " << Param::sigma << "\n";
  cout << "Box size = [" << Param::box[0] << ", " << Param::box[1] << "]\n";
  cout << "Number of cells = " << Param::number_cell[0] << ", "
	   << Param::number_cell[1] << "\n";

  ifstream out_file(Param::out_string);
  ifstream in_file(Param::in_string);
  if (in_file.is_open()) {
	cout << "Detected initial configuration\n";
	read_h5();
  } else {
	cout << "Start a fresh run with generated initial configuration\n";
	init_pos();
  }
  // Checking overlap
  // For large number_disks, it is strongly recommended that the check is
  // disabled due to its O(number_disks^2) complexity of check_overlap().
  if (check_overlap()) cerr << "Overlap in initial configuration!" << endl;
  init_cell();
  create_h5();
  // Writes initial configuration.
  out_configuration_h5(-1);
}

StraightEcmc::~StraightEcmc() {
  cout.precision(5);
  cout << number_collisions << " collisions, "
	   << double(time_end - time_start) / (double)CLOCKS_PER_SEC
	   << " seconds\n";
  double eph = (double)number_collisions * (double)CLOCKS_PER_SEC
	  / double(time_end - time_start) * 3600.;
  cout << "Estimated number of events per hour: " << eph << endl;
  finish_h5(eph);
  cout.precision(11);
}
