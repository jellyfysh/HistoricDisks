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

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSpace.hpp>
#include "StraightEcmc.h"

using namespace std;
using namespace HighFive;
using namespace H5Easy;

void StraightEcmc::create_h5() {
  // Create the output hdf5 file and write into it the basic information of the run.
  string string_h5 = Param::out_string;
  File fd(string_h5, File::ReadWrite | File::Create | File::Truncate);
  dump(fd, "/parameters/Nx", Param::number_disks_x);
  dump(fd, "/parameters/Ny", Param::number_disks_y);
  dump(fd, "/parameters/N", Param::number_disks);
  dump(fd, "/parameters/eta", Param::eta);
  DataSet fd_L = fd.createDataSet<double>("/parameters/L", DataSpace(2));
  fd_L.write(Param::box);
  dump(fd, "/parameters/sigma", Param::sigma);
  dump(fd, "/parameters/slant", Param::slant);
  dump(fd, "/parameters/shape", Param::shape);
  auto t = time(nullptr);//now write the time as a string
  auto tm = *localtime(&t);
  ostringstream oss;
  oss << put_time(&tm, "%d-%m-%Y %H:%M:%S");
  dump(fd, "/stats/start_time", oss.str());
  DataSpace ds = DataSpace({1}, {DataSpace::UNLIMITED});
  DataSetCreateProps props;
  props.add(Chunking(vector<hsize_t>{1}));
  DataSet
	  pressure = fd.createDataSet("/pressure", ds, AtomicType<double>(), props);
  DataSet pressure_x =
	  fd.createDataSet("/pressure_x", ds, AtomicType<double>(), props);
  DataSet pressure_y =
	  fd.createDataSet("/pressure_y", ds, AtomicType<double>(), props);
}

void StraightEcmc::read_h5() {
  // Reads initial configuration stored in a hdf5 file.
  string string_h5 = Param::in_string;
  File fd(string_h5, File::ReadOnly);
  vector<vector<double>>
	  configuration = load<vector<vector<double>>>(fd, "/config-init");
  for (unsigned int i = 0; i < Param::number_disks; i++) {
	positions[i][0] = configuration[i][0];
	positions[i][1] = configuration[i][1];
	if (positions[i][0] > Param::box[0] / 2) positions[i][0] -= Param::box[0];
	if (positions[i][1] > Param::box[1] / 2) positions[i][1] -= Param::box[1];
  }
  cout << "Finish reading config\n";
}

void StraightEcmc::finish_h5(double eph) {
  // Writes the end time of the run and the estimated number of collisions per
  // hour into the output hdf5 file.
  string string_h5 = Param::out_string;
  File fd(string_h5, File::ReadWrite);
  auto t = time(nullptr);
  auto tm = *localtime(&t);
  ostringstream oss;
  oss << put_time(&tm, "%d-%m-%Y %H:%M:%S");
  dump(fd, "/stats/end_time", oss.str(), DumpMode::Overwrite);
  dump(fd, "/stats/EPH", eph, DumpMode::Overwrite);
}

void StraightEcmc::out_configuration_h5(int counter) {
  // Writes the configuration into the output hdf5 file
  // Updates the configuration count "/count" in the output hdf5 file.
  // "index" is the configuration count during the run.
  // When index < -1, the configuration is written into the output hdf5 file
  // bearing the name of the initial configuration "config-init".
  string string_h5 = Param::out_string;
  File fd(string_h5, File::ReadWrite);
  string str_con = "config-" + to_string(counter);
  if (counter < 0) str_con = "config-init";
  vector<size_t> dims{Param::number_disks, 2};
  DataSet config = fd.createDataSet<double>(str_con, DataSpace(dims));
  vector<vector<double>> pos(Param::number_disks, vector<double>(2));
  unsigned int count = 0;
  for (int i = 0; i < Param::total_number_cell; i++)
	for (int k = 0; k < cell_ocp[size_t(i)]; k++) {
	  pos[count][0] = cell[i][k][0]
		  + (i % Param::number_cell[0] + 0.5) * Param::cell_size[0];
	  pos[count][1] = cell[i][k][1]
		  + (floor(i / Param::number_cell[0]) + 0.5) * Param::cell_size[1];
	  count++;
	}
  config.write(pos);
  dump(fd, "/count", counter, DumpMode::Overwrite);
  dump(fd, "/stats/collisions", number_collisions, DumpMode::Overwrite);
}

void StraightEcmc::out_series_h5(unsigned int size,
								 double value,
								 string name) {
  string string_h5 = Param::out_string;
  File fd(string_h5, File::ReadWrite);
  DataSet dataset = fd.getDataSet(name);
  dataset.resize({size + 1});
  dataset.select({size}, {1}).write(value);
}
