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

// Straight ECMC with a single active disk
// The algorithm is identical to Etienne's fortran program

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
#include <chrono>
#include <random>
#include "StraightEcmc.h"

using namespace std;

namespace Param {
unsigned int number_disks = 0;
unsigned int number_disks_x = 0;
unsigned int number_disks_y = 0;
double eta = 0;
double box[2] = {1.0, 1.0};
int number_cell[2] = {0, 0};
int total_number_cell = number_cell[0] * number_cell[1];
double lambda_0 = 0;
double sigma = 0;
double cell_size[2] =
	{Param::box[0] / number_cell[0], Param::box[1] / number_cell[1]};
int slant = 0;
int extra_factor = 1;
std::string out_string = "output.h5";
std::string in_string = "dummy.h5";
int shape = 0;
// By default, the seed for the RNG is 0.
mt19937 random_generator = mt19937(0);
}

int main(int argc, char *argv[]) {
  if (argc != 9) {
	cout << "Bad number of parameters, quitting.\n";
	return 0;
  }
  cout << "Initializing with command-line arguments.\n";
  Param::number_disks_x =
	  static_cast<unsigned int>(strtol(argv[1], nullptr, 10));
  Param::number_disks_y =
	  static_cast<unsigned int>(strtol(argv[2], nullptr, 10));
  Param::number_disks = Param::number_disks_x * Param::number_disks_y;
  Param::eta = strtod(argv[3], nullptr);
  Param::slant = (int)strtol(argv[4], nullptr, 10);
  Param::extra_factor = (int)strtol(argv[5], nullptr, 10);
  Param::out_string = "";
  Param::out_string += argv[7];
  Param::in_string = "";
  Param::in_string += argv[8];
  // The box volume is always 1
  if (strcmp(argv[6], "square") == 0) {
	Param::shape = 1;
	Param::number_cell[0] = (int)(sqrt(Param::number_disks) * 7.0 / 8.0);
	Param::number_cell[1] = (int)(sqrt(Param::number_disks) * 7.0 / 8.0);
  } else if (strcmp(argv[6], "rectangle") == 0) {
	Param::shape = 2;
	Param::box[0] = 1.0 / sqrt(sqrt(3.0) / 2.0);
	Param::box[1] = sqrt(sqrt(3.0) / 2.0);
	Param::number_cell[0] =
		(int)(sqrt(Param::number_disks) * Param::box[0] * 7.0 / 8.0);
	Param::number_cell[1] =
		(int)(sqrt(Param::number_disks) * Param::box[1] * 7.0 / 8.0);
  } else if (strcmp(argv[6], "crystal") == 0) {
	Param::shape = 0;
	Param::box[0] = 1.0
		/ sqrt(sqrt(3.0) / 2.0 * Param::number_disks_y / Param::number_disks_x);
	Param::box[1] =
		sqrt(sqrt(3.0) / 2.0 * Param::number_disks_y / Param::number_disks_x);
	Param::number_cell[0] = (int)(Param::number_disks_x * 7.0 / 8.0);
	Param::number_cell[1] = (int)(Param::number_disks_y * 7.0 / 8.0);
  } else {
	cout << "Bad shape, quitting.\n";
	return 0;
  }
  Param::sigma = sqrt(
	  Param::box[0] * Param::box[1] * Param::eta / 3.14159265358979323846
		  / double(Param::number_disks));
  //I n order to avoid self collision, there must be at least two cells in each
  // direction.
  Param::number_cell[0] = max(2, Param::number_cell[0]);
  Param::number_cell[1] = max(2, Param::number_cell[1]);
  Param::total_number_cell = Param::number_cell[0] * Param::number_cell[1];

  Param::lambda_0 = 0.07680 / sqrt(Param::number_disks);
  Param::cell_size[0] = Param::box[0] / Param::number_cell[0];
  Param::cell_size[1] = Param::box[1] / Param::number_cell[1];

  // Uncomment the following line to seed the RNG by the time the run starts.
  /* chrono::high_resolution_clock::time_point
	  t1 = chrono::high_resolution_clock::now();
  auto t_int =
	  chrono::duration_cast<chrono::nanoseconds>(t1.time_since_epoch()).count();
  cout << "seed: " << t_int << endl;
  Param::random_generator = mt19937(static_cast<unsigned int>(t_int)); */

  StraightEcmc c;
  cout << "Run created\n";
  c.run();
  return 0;
}
