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
Functions involved in extracting data from a hdf5 file.
"""
import h5py
import sys
from functools import partial
from typing import Union, Any


def find_item(name: str, target: str) -> Union[str, None]:
    """
    Function used in the visit() method of the hdf5 file object.
    Checks whether the name of a dataset is contained in a path.

    Parameters
    ----------
    name : str
        Trial path.
    target :
        Name of the dataset to be found.

    Returns
    -------
    Union[str, None]
        Trial path if it contains the target. In other cases, None.
    """
    if target in name:
        return name


def extract_item(file: h5py.File, target: str) -> Any:
    """
    Extract the content of a dataset from a hdf5 file.

    Parameters
    ----------
    file: h5py.File
        The hdf5 file.
    target: str
        Name of the dataset to be found.

    Returns
    -------
    Any
        The content of the dataset.
    """
    path = file.visit(partial(find_item, target=target))
    output = file[path]
    output = output[()]
    return output


if __name__ == '__main__':
    # Examples of extracting information from a hdf5 data file generated by the C++ ECMC program.
    # Example usage: "python3 h5_read.py example.h5"
    f = h5py.File(sys.argv[1], 'r')
    N = extract_item(f, 'N')
    print('N', N)
    Nx = extract_item(f, 'Nx')
    print('Nx', Nx)
    Ny = extract_item(f, 'Ny')
    print('Ny', Ny)
    eta = extract_item(f, 'eta')
    print('eta', eta)
    L = extract_item(f, 'L')
    print('L', L)
    sigma = extract_item(f, 'sigma')
    print('sigma', sigma)
    slant = extract_item(f, 'slant')
    print('slant', slant)
    shape = extract_item(f, 'shape')
    print('shape', shape)
    start_time = extract_item(f, 'start_time')
    print('start_time', start_time)
    end_time = extract_item(f, 'end_time')
    print('end_time', end_time)
    EPH = extract_item(f, 'EPH')
    print('EPH', EPH)
    pressure = extract_item(f, 'pressure')
    print('pressure', pressure)
    config = extract_item(f, 'config-0')
    print('config', config)
