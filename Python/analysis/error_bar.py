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
Use stationary bootstrap to calculate the error bar of the average of a time series.
The stationary bootstrap is performed using the implementation https://github.com/YoshihikoNishikawa/StationaryBootstrap.
"""
import sys
from typing import Tuple
import numpy as np
from stresampling import stationary_bootstrap as sbm


def error_bar(time_series: np.ndarray, alpha: float = 0.68) -> Tuple[float, float]:
    """
    Find the error bar of a time series.

    Parameters
    ----------
    time_series : np.ndarray
        The series whose error bar is calculated.
    alpha : float
        Confidence level. Set to 0.68 by default.

    Returns
    -------
    Tuple[float, float]
        The average of the series and the standard error of the mean value of the series.
    """
    stat = sbm.conf_int(time_series, np.mean, alpha)
    return stat.mean, stat.se


if __name__ == '__main__':
    # Example usage: "python3 error_bar.py pressure.dat"
    series = np.loadtxt(sys.argv[1])
    mean, se = error_bar(series)
    print("{} +/- {}".format(mean, se))
