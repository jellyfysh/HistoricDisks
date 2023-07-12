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

#ifndef _StraightEcmc_h
#define _StraightEcmc_h

#include "boost/multi_array.hpp"
#include <unordered_map>
#include "Param.h" //set all parameters here

using namespace std;

class StraightEcmc {
 private:
  // total number of collisions
  long int number_collisions;

  // Wall time of the start and the end
  clock_t time_start{}, time_end{};

  // The cell containing the active disk
  int cell_cur{};

  // The index of the active disk in its cell
  int k_cur{};
  // cell_cur, k_cur also used as disk being treated when filling the cells
  // during initialization.

  // Number of disks in each cell
  vector<int> cell_ocp;

  // Disk positions
  // Only used before the cell system is built.
  boost::multi_array<double, 2> positions;

  // Cell system
  // Disk position is stored as relative position with respect to the center
  // of each cell.
  boost::multi_array<double, 3> cell;

  // Neighbour list
  // 2D array of size number_cells * 9
  // elements are indices of the cells
  boost::multi_array<int, 2> cell_neighbour;

  /// Generate initial configuration.
  void init_pos();

  /** Initializes the cell system.
  Called after init_pos() or reading in the initial configuration. */
  void init_cell();

  /** Schedules the collision when the active disk moves in the x direction.
   *
   * Searches six out of the nine neighboring cells to find possible collision
   * between the active disk and the static ones.
   *
   * @param[in, out] L_min the time before the collision takes place. The rest
   * of the chain length is passed into the function through this parameter. If
   * a collision can take place, its value is updated as the time before the
   * collision taking place.
   * @param[out] k_min the index of the disk at the collision in its cell
   * @param[out] cell_min the index of the cell containing the disk in the
   * coming collision
   * @param[out] delta_x_ij the $\ Delta x_{ij}$ of the disks at the
   * collision; used in pressure computation.
   * */
  void expl_cell_x(double &L_min,
				   int &k_min,
				   int &cell_min,
				   double &delta_x_ij);

  /** Schedules the collision when the active disk moves in the y direction.
   *
   * Searches six out of the nine neighboring cells to find possible collision
   * between the active disk and the static ones.
   *
   * @param[in, out] L_min the time before the collision takes place. The rest
   * of the chain length is passed into the function through this parameter. If
   * a collision can take place, its value is updated as the time before the
   * collision taking place.
   * @param[out] k_min the index of the disk at the collision in its cell
   * @param[out] cell_min the index of the cell containing the disk in the
   * coming collision
   * @param delta_y_ij the $\ Delta y_{ij}$ of the disks at the
   * collision; used in pressure computation.
   * */
  void expl_cell_y(double &L_min,
				   int &k_min,
				   int &cell_min,
				   double &delta_y_ij);

  /** Updates the cell system when the active disk moves in x direction.
   *
   * Update the cell-occupation numbers and the position of the active disk.
   * Called after the new position of the active disk found by expl_cell_x().
   *
   * \param dx the displacement of the active disk in x direction
   * */
  void refresh_cell_x(double dx);

  /** Updates the cell system when the active disk moves in y direction.
   *
   * Update the cell-occupation numbers and the position of the active disk.
   * Called after the new position of the active disk found by expl_cell_y().
   *
   * \param dy the displacement of the active disk in y direction
   * */
  void refresh_cell_y(double dy);

  /** Uniformly samples a disk from the cell system.
   *
   * This function is used to choose active disk when resampling.
   *
   * @param[out] cell the index of the cell containing the chosen disk
   * @param[out] k the index of the chosen disk in its cell
   * */
  void choose_active(int &cell, int &k);

  /// Initializes the output hdf5 file.
  void create_h5();

  /// Reads initial configuration stored in a hdf5 file.
  void read_h5();

  /** Writes the end time of the run and the estimated number of collisions per
   * hour into the output hdf5 file.
   *
   * @param eph the estimated number of collisions per hour
   * */
  void finish_h5(double eph);

  /** Writes the disk configuration into the output hdf5 file.
   * Also updates the configuration count "/count" in the hdf5 file.
   *
   * @param counter the number of sampled configurations minus 1
   * */
  void out_configuration_h5(int counter);

  /** Appends a new element to a list in the hdf5 file.
   *
   * @param size the size of the list
   * @param value the value to be added
   * @param name the name of the list
   * */
  void out_series_h5(unsigned int size, double value, string name);

  /** Finds the square of the distance between two disks factoring in the
   * periodic boundary condition.
   *
   * Used before the cell system is built.
   *
   * @param i the index of a disk
   * @param j the index of another disk
   *
   * @return the minimum distance between the disks.
   * */
  double distance_square(unsigned int i, unsigned int j);

 public:
  StraightEcmc();

  ~StraightEcmc();

  /// Runs the simulation.
  void run();

  /** Loops over all disk pairs to find overlap.
   *
   * @return boolean value indicating whether at least an overlap is found
   * */
  bool check_overlap();
};

#endif
