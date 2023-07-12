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
Executable python script that analyzes a simulation stored in a hdf5 file.
If psi6 is not in the file, this script calculates psi6 and adds it to the file.
If psi6 is in the file, this script reads psi6 from the file.
The evolution of psi6 is visualized in "psi6.pdf", and the correlation between the angle of psi6 and pressure is
visualized in "pressure_vs_psi6.pdf"
The error bar of the pressure is estimated using stationary bootstrap.
Usage example: "python3 analysis_h5.py example.h5"
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
from h5_read import extract_item
from psi6 import global_psi6
from error_bar import error_bar


if __name__ == '__main__':
    f = h5py.File(sys.argv[1], 'r+')
    box = extract_item(f, 'L')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Abandoning the beginning of the run. By default, the first 30% is abandoned.
    start_from = int(0.3 * (extract_item(f, 'count') + 1))
    pressure = extract_item(f, 'pressure')[start_from:]
    psi6_snapshot = []
    # If psi6 is already in the hdf5 data file, this script reads in the psi6.
    if 'psi6' in f:
        psi6_snapshot = extract_item(f, 'psi6')[start_from:]
        print('Psi^6 is already calculated for this run.')
    # If psi6 is not in the data file, this script calculates psi6 and adds it into the data file.
    else:
        for i in range(extract_item(f, 'count') + 1):
            configuration = extract_item(f, 'config-' + str(i))
            psi6_snapshot.append(global_psi6(configuration, box))
        psi6_snapshot = np.array(psi6_snapshot)
        f.create_dataset("psi6", data=psi6_snapshot)
        psi6_snapshot = psi6_snapshot[start_from:]

    print("Average Psi^6 during the run: {}".format(np.average(psi6_snapshot)))

    # Visualization of the evolution of psi6.
    ax.scatter(psi6_snapshot.real, psi6_snapshot.imag, marker='o', s=2, color='C0')
    ax.plot(psi6_snapshot.real, psi6_snapshot.imag, alpha=0.2, color='C0')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal')
    plt.savefig('psi6.pdf')
    plt.close()

    # Visualization of the correlation between psi6 and pressure.
    plt.scatter(np.angle(psi6_snapshot) / np.pi, pressure, marker='o', s=2)
    pressure_average, pressure_se = error_bar(pressure)
    print("Pressure: {} +/- {}".format(pressure_average, pressure_se))
    plt.plot([-1.5, 1.5], [pressure_average, pressure_average], color='k')
    plt.xlabel(r'arg($\Psi_6$) / $\pi$')
    plt.ylabel(r'$\beta P V_0 / N$')
    plt.xlim(-1, 1)
    plt.savefig('pressure_vs_psi6.pdf')
    plt.close()
