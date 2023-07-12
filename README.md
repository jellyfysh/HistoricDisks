[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

# HistoricDisks

This repository accompanies the work
[Hard-disk pressure computations—a historic perspective](https://doi.org/10.1063/5.0126437). It provides the
pressure data of the hard-disk model extracted from the literature since 1953, and also the set of high-precision
pressures that are presented in the work above. Furthermore, this repository contains naive Markov-chain Monte Carlo
and molecular dynamics implementations (with or without pressure estimation) for the hard-disk model, as well as a
state-of-the-art implementation of the hard-disk event-chain Monte Carlo algorithm.

## List of programs

This repository contains:

* Pressure data, equations of states
    * Pressure data files (CSV, see the [DigitizedData](DigitizedData) directory)
    * Equations of states visualization (Python, requires `matplotlib`, see the [synopsis](Python/synopsis) directory)

* Four-disk non-periodic-box programs (Python)
    * Sampling program using Metropolis algorithm (Python, see the
      [Python/four-disk/Metropolis_disks_box.py](Python/four-disk/Metropolis_disks_box.py) script)
    * Sampling program using Molecular dynamics with pressure estimators (Python, see the
      [Python/four-disk/molecular_disks_box.py](Python/four-disk/molecular_disks_box.py) script)
    * Sampling program using straight ECMC with pressure estimators (Python, see the
      [Python/four-disk/ECMC_straight_disks.py](Python/four-disk/ECMC_straight_disks.py) script)
    * Sampling program using reflective ECMC (Python, see the
      [Python/four-disk/ECMC_reflective_disks_box.py](Python/four-disk/ECMC_reflective_disks_box.py) script)
    * Sampling program using forward ECMC (Python, see the
      [Python/four-disk/ECMC_forward_disks_box.py](Python/four-disk/ECMC_forward_disks_box.py) script)
    * Sampling program using Newtonian ECMC (Python, see the
      [Python/four-disk/ECMC_Newtonian_disks_box.py](Python/four-disk/ECMC_Newtonian_disks_box.py) script)

* Naive periodic-box programs
    * Sampling program using Metropolis algorithm (Python, see
      the [Python/naive/Metropolis.py](Python/naive/Metropolis.py) script)
    * Sampling program using Molecular dynamics (Python, see
      the [Python/naive/molecular_dynamics.py](Python/naive/molecular_dynamics.py) script)
    * Sampling program using straight ECMC with pressure estimators (Python, see
      the [Python/naive/ECMC_straight.py](Python/naive/ECMC_straight.py) script)
    * Sampling program using reflective ECMC (Python, see
      the [Python/naive/ECMC_reflective.py](Python/naive/ECMC_reflective.py) script)
    * Sampling program using forward ECMC (Python, see the [Python/naive/ECMC_forward.py](Python/naive/ECMC_forward.py)
      script)
    * Sampling program using Newtonian ECMC (Python, see
      the [Python/naive/ECMC_Newtonian.py](Python/naive/ECMC_Newtonian.py) script)

* State-of-the-art hard-disk programs
    * Sampling program using straight ECMC with pressure estimators (C++, see the [CPP/StraightECMC](CPP/StraightECMC)
      directory)

* Analysis
    * Pressure calculation using the fitting formula (Python, see the
      [Python/four-disk/fitting.py](Python/four-disk/fitting.py) script)
    * Global Orientational order parameter calculation (Python, requires `numpy`, `scipy`, see
      the [Python/analysis/psi6.py](Python/analysis/psi6.py) script)
    * Confidence-interval calculation (Python, requires `numpy`, `stresampling`, see
      the [Python/analysis/error_bar.py](Python/analysis/error_bar.py) script)
    * Analysis for straight ECMC runs stored in the HDF5 format (Python, requires `numpy`, `stresampling`, `h5py`, see
      the [Python/analysis/analysis_h5.py](Python/analysis/analysis_h5.py) script)

## Installing

All Python scripts can be executed with any Python3 implementation (e.g., standard [CPython](https://www.python.org) or
[PyPy3](https://www.pypy.org)).

The sampling programs do not have any further requirements. The script for the
pressure calculation using the fitting formulas (see [Python/four-disk/fitting.py](Python/four-disk/fitting.py)) relies
on [NumPy](https://numpy.org).

In addition to NumPy, the analysis programs rely on [h5py](https://www.h5py.org) to read and store data in the HDF5 data
format. [SciPy](https://scipy.org) is used to compute the orientational order parameter of a hard-disk configuration,
and [matplotlib](https://matplotlib.org) is used for visualization. The
[stresampling](https://pypi.org/project/stresampling/) package estimates error bars by using stationary bootstrap.

The requirements can also be found in the [List of Programs](#list-of-programs). All external dependencies for the
Python programs are included in [requirements.txt](requirements.txt). 

We recommend setting up a virtual environment and installing the correct versions of the external dependencies by
running the following commands (replace `python3` by the Python3 interpreter of your choice):

```shell
python3 -m venv historic_venv
source historic_venv/bin/activate
python3 -m pip install -U pip setuptools
python3 -m pip install -r requirements.txt
```

Please refer to
the [README.md](CPP/StraightECMC/README.md) of the C++ program for its compilation requirements.

## Authors

Check the [AUTHORS.md](AUTHORS.md) file to see who participated in this project.

## License

This project is licensed under the GNU General Public License, version 3 (see the [LICENSE](LICENSE) file).

## Contact

If you have questions regarding the HistoricDisks software package, just raise an issue here on GitHub or contact us via
mail (see the [AUTHORS.md](AUTHORS.md) file). We are happy to help you!

## Citation

If you use (parts of the) the programs in published work, please cite the following reference (see
[[Li2022]](https://doi.org/10.1063/5.0126437) in [References.bib](References.bib)):

Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, Werner Krauth;\
Hard-disk pressure computations&mdash;a historic perspective.\
J. Chem. Phys. 21 December 2022; 157 (23): 234111. https://doi.org/10.1063/5.0126437
