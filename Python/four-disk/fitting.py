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
# Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, and Werner Krauth,
# Hard-disk computer simulations---a historic perspective,
# arXiv e-prints: 2207.07715 (2022), https://arxiv.org/abs/2207.07715.
#
"""
Executable Python script that calculates the pressure from hard-disk configurations using Eqs (12) and (27a) in
[Li2022].

The script reads in hard-disk configurations in the output format of the sampling programs. The filename should be the
only positional argument of this script (i.e., one should use the command 'python3 fitting.py configuration.txt'). This
script relies on NumPy as an external dependency.

The script considers the distances between all disk pairs and the distances between disks and the walls where the
distance is small enough, i.e., the shifted distance such that the minimum distance is 0 is smaller than the fit
interval. Here, The size of the fit interval is 0.1 * sigma. The fit interval is then divided into 100 equal sized bins
that serve as the bins for the histogram. The pair-distance histogram is then processed to be the pair-correlation
function. The pair-correlation function is fitted by a fourth-order polynomial. The value of the pair-correlation
function at 2 * sigma is obtained by an extrapolation of the polynomial fit. The value of the rescaled line density is
obtained with the same procedure.

This script prints the pressures and the corresponding error bars calculated from Eqs (12) and (27a) in [Li2022].

The fitting procedure is defined as a class, initialized with configurations, fit interval, bin size, number of disks,
box geometry, and radius. The operations on the configurations are defined as methods in the class.
"""
import sys
from typing import List, Sequence
import numpy as np


class Fitting:
    """
    Class that provides methods to calculate the extrapolated pair-correlation function and rescaled line density from
    hard-disk configurations.

    Attributes
    ----------
    n : int
        The number of hard disks.
    sigma : float
        The radius of the hard disks.
    box : Sequence[float]
        The geometry of the simulation box, i.e., the side lengths L_x and L_y.
    bin_size : float
        The bin size for the histogram.
    fit_interval : float
        The maximum considered (shifted) distance in the wall and pair distances.
    wall_distances : List[List[float]]
        The list of wall distances for both directions shifted by sigma such that the minimum distance is 0.
        Any shifted distance larger than fit_interval is excluded.
    wall_sample_size : int
        Number of all wall distances in either of the directions
    pair_distances_sq : List[float]
        The list of squared pair distances shifted by (2 * sigma) ** 2 such that the minimum distance is 0.
        Any shifted distance larger than fit_interval is excluded.
    pair_sample_size : int
        Number of all pair distances.
    """
    def __init__(self, fit_interval: float, bin_size: float, n: int, sigma: float, box: Sequence[float]):
        """
        Initialize the instance by storing the relevant parameters.

        Parameters
        ----------
        fit_interval : float
            The maximum considered (shifted) distance in the wall and pair distances.
        bin_size : float
            The bin size for the histogram.
        n : int
            The number of hard disks.
        sigma : float
            The radius of the hard disk.
        box : Sequence[float]
            The geometry of the simulation box, i.e., the side lengths L_x and L_y.
        """
        self.n = n
        self.sigma = sigma
        self.box = box
        self.bin_size = bin_size
        self.fit_interval = fit_interval
        self.wall_distances = [[], []]
        self.wall_sample_size = 0
        self.pair_distances_sq = []
        self.pair_sample_size = 0

    @staticmethod
    def load_configurations(filename: str) -> List[List[float]]:
        """
        Load the hard-disk configurations from the given file.

        Each line of the file contains a single hard-disk configuration. The (2 * k)th and (2 * k + 1)th floats in the
        line should be the x- and y-positions of the kth disk, respectively.

        Parameters
        ----------
        filename : str
            The name of the file that stores the hard-disk configurations.

        Returns
        -------
        List[List[float]]
            The hard-disk configurations.
        """
        configurations = []
        with open(filename, "r") as file:
            for line in file:
                configurations.append(list(map(float, line.split())))
        return configurations

    def compute_wall_distances(self, configurations: List[List[float]]) -> None:
        """
        Compute and store the wall distances shifted by sigma from the given hard-disk configurations. Only shifted
        distances smaller than self.fit_interval are included.

        The required format of the hard-disk configurations is documented in the static self.load_configurations method.

        Parameters
        ----------
        configurations : List[List[float]]
            The hard-disk configurations.
        """
        for configuration in configurations:
            for i in range(self.n):
                for d in range(2):
                    position = configuration[2 * i + d]
                    if position < self.fit_interval + self.sigma:
                        self.wall_distances[d].append(position - self.sigma)
                    elif position > self.box[d] - self.sigma - self.fit_interval:
                        self.wall_distances[d].append(self.box[d] - position - self.sigma)
        self.wall_sample_size = len(configurations) * self.n

    @staticmethod
    def distance_sq(disk_one: Sequence[float], disk_two: Sequence[float]) -> float:
        """
        Return the squared distance between two hard disks.

        Parameters
        ----------
        disk_one : Sequence[float]
            The (two-dimensional) position of the first hard disk.
        disk_two : Sequence[float]
            The (two-dimensional) position of the second hard disk.

        Returns
        -------
        float
            The squared distance.
        """
        return (disk_one[0] - disk_two[0]) ** 2 + (disk_one[1] - disk_two[1]) ** 2

    def compute_distances_sq(self, configurations: List[List[float]]) -> None:
        """
         Compute and store the squared pair distances shifted by (2 * sigma) ** 2 from the given hard-disk
         configurations. Only shifted distances smaller than self.fit_interval are included.

        The required format of the hard-disk configurations is documented in the static self.load_configurations method.

        Parameters
        ----------
        configurations : List[List[float]]
            The hard-disk configurations.
        """
        criterion = (2. * self.sigma + self.fit_interval) ** 2
        for configuration in configurations:
            for i in range(self.n):
                for j in range(i):
                    square_distance = self.distance_sq([configuration[2 * i], configuration[2 * i + 1]],
                                                       [configuration[2 * j], configuration[2 * j + 1]])
                    if square_distance < criterion:
                        self.pair_distances_sq.append(square_distance)
        self.pair_sample_size = len(configurations) * self.n * (self.n - 1) / 2

    def fit_rho(self, direction: int) -> float:
        """
        Create the histogram of the shifted wall distances in the given direction, and obtain the approximate rescaled
        line density. Then, fit a fourth-order polynomial for the final extrapolation to the shifted distance 0.

        Parameters
        ----------
        direction : int
            The direction (0 for x- and 1 for y-direction).

        Returns
        -------
        float
            The rescaled line density extrapolated to the shifted distance 0.

        Raises
        ------
        AssertionError
            If the direction is not 0 or 1.
        """
        assert direction == 0 or direction == 1
        hist, bins = np.histogram(self.wall_distances[direction], np.arange(0, self.fit_interval, self.bin_size))
        pdf = [h / self.wall_sample_size / self.bin_size * self.box[direction] for h in hist]
        r = np.arange(bins[0] + self.bin_size / 2, bins[-1] + self.bin_size / 2, self.bin_size)
        poly = np.polyfit(r, pdf, 4)
        p = np.poly1d(poly)
        return p(0) / 2.0

    def fit_g(self) -> float:
        """
        Create the histogram of the shifted squared pair distances, and obtain the approximate pair-correlation
        function. Then, fit a fourth-order polynomial for final extrapolation to the shifted distance 0.

        Returns
        -------
        float
            The pair-correlation function extrapolated to the shifted distance 0.
        """
        bins = np.arange(0, self.fit_interval, self.bin_size)
        bins_sq = (bins + 2 * self.sigma) ** 2
        hist, _ = np.histogram(self.pair_distances_sq, bins_sq)
        pdf = [h / self.pair_sample_size / self.bin_size * self.box[0] * self.box[1] for h in hist]
        pdf = [pdf[i] / 2 / np.pi / (bins[i] + 2 * self.sigma + self.bin_size / 2) for i in range(len(pdf))]
        r = np.arange(bins[0] + self.bin_size / 2, bins[-1] + self.bin_size / 2, self.bin_size)
        poly = np.polyfit(r, pdf, 4)
        p = np.poly1d(poly)
        return p(0)


def main() -> None:
    """
    Read the hard-disk configurations from the file given by the first positional argument to this script, and compute
    the pressures and the corresponding error bars calculated from Eqs (12) and (27a) in [Li2022]. The error bars are
    estimated from computing a pressure estimate for batches of the hard-disk configurations.
    """
    configurations = Fitting.load_configurations(sys.argv[1])
    n = 4
    number_batch = 100
    box = [1.0, 1.0]
    sigma = 0.15
    fit_interval = 0.1 * sigma
    bin_size = 0.01 * fit_interval
    batch_size = len(configurations) // number_batch
    pressure_wall = []
    pressure_homothetic = []
    for i in range(number_batch):
        selected_configurations = configurations[batch_size * i: batch_size * (i + 1)][:]
        fit = Fitting(fit_interval, bin_size, n, sigma, box)
        fit.compute_wall_distances(selected_configurations)
        fit.compute_distances_sq(selected_configurations)
        rho_x = fit.fit_rho(0)
        rho_y = fit.fit_rho(1)
        g = fit.fit_g()
        pressure_wall.append(n * (rho_x + rho_y) / 2.0)
        pressure_homothetic.append(n + n * (2 * np.pi * (n - 1) * sigma ** 2 * g + sigma * (rho_x + rho_y)))
    pressure_wall = np.array(pressure_wall)
    pressure_homothetic = np.array(pressure_homothetic)
    print(r'Wall pressure: {:+.6f} \pm {:+.6f}'.format(
        pressure_wall.mean(), pressure_wall.std(ddof=1) / np.sqrt(len(pressure_wall))))
    print(r'Homothetic pressure: {:+.6f} \pm {:+.6f}'.format(
        pressure_homothetic.mean(), pressure_homothetic.std(ddof=1) / np.sqrt(len(pressure_homothetic))))


if __name__ == '__main__':
    main()
