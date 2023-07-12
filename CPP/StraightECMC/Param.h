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

#ifndef _Param_h
#define _Param_h

#include <random>

using namespace std;

namespace Param {
// Multiplier of standard chain length, controls the total run time
// 200000000 means roughly 20s on a laptop with Intel CORE i7 9th Gen.
const long int factor = 200000000;

// Number of sampled configurations
const int n_samples = 1000;

// Max number of disks in a cell
const int n_cell_max = 5;

// Path of the output hdf5 file
extern std::string out_string;

// Path of the hdf5 containing the initial configuration.
extern std::string in_string;

// Density (packing fraction)
extern double eta;

// Box's geometry
extern double box[2];

// Number of disks
extern unsigned int number_disks;

// For crystalline initial configuration, the number of disks in the x direction
// and y direction is specified
// The crystalline initial configuration is a "base" configuration, in which the
// edge of the lattice is aligned with the x-axis. For example:
//  o o o o o o o
// o o o o o o o
//  o o o o o o o
// o o o o o o o
// is a base configuration with number_disks_x = 7 and number_disks_y = 4

// Number of disks in a row in the lattice
extern unsigned int number_disks_x;

// Number of rows, the distance between two rows are \sqrt{3} / 2 * the distance
// between two disks in a row
// number_disks_y has to be an even number to be compatible with the
// periodic boundary condition
extern unsigned int number_disks_y;

// Number of cells in x and y directions
extern int number_cell[2];

// Total number of cells
extern int total_number_cell;

// Mean free path
extern double lambda_0;

// Radius of the disk
extern double sigma;

// Length of the edges of the cell in x and y directions
extern double cell_size[2];

// Multiplied with factor to control the run time externally
extern int extra_factor;

// Parameter controls slant for crystalline initial configuration.
// Slant means raise the disks in a row gradually such that the end of the row
// connects with the beginning another row at the periodic boundary.
// The example can be found in
// Statistical Mechanics: Algorithms and Computations, Fig. 2.26.
// Variable slant controls how rows are matched.
// For example, value 2 means the end of the first row matched the beginning of
// the third row.
// Slant could only be even number
extern int slant;

// Integer indicating the shape of the system.
// 0 means crystalline aspect ratio; 1 means square;
// 2 means aspect ratio \sqrt{3} / 2.
extern int shape;

// Random number generator
extern mt19937 random_generator;
}

#endif
