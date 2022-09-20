[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

# HistoricDisks&mdash;Preview

This repository accompanies the work 
[Hard-disk computer simulations&mdash;a historic perspective](https://arxiv.org/abs/2207.07715). It provides the 
pressure data of the hard-disk model extracted from the literature since 1953, and also the set of high-precision 
pressures that are presented in the work above. Furthermore, this repository contains naive Markov-chain Monte Carlo 
and molecular dynamics implementations for the hard-disk model, several pressure estimators, as well as a 
state-of-the-art implementations of the hard-disk event-chain Monte Carlo algorithm.

## Version 0.2

In its current version 0.2, this repository serves as a preview for the complete version 1 that is expected to be 
finished in September 2022. In order to demonstrate the overall style of the programs, this preview version contains 
sampling programs for four hard disks in a non-periodic box that use the Metropolis algorithm, event-driven 
molecular dynamics, or four variants (straight, reflective, forward, and Newtonian) of the event-chain Monte 
Carlo (ECMC) algorithm. The preview version also contains naive sampling programs for an arbitrary number of hard 
disks in a periodic box using all algorithms implemented in the four-disk programs. We also provide a Python script 
that computes the pressure from hard-disk configurations by fitting a polynomial to the estimated pair-correlation 
function and the rescaled line density.

## List of programs
The full repository will contain:

- [ ] Pressure data, equations of states
   - [ ] Pressure data files (CSV)
   - [ ] Equations of states visualization (Python)

- [x] Four-disk non-periodic-box programs (Python)
   - [x] Sampling program using Metropolis algorithm (Python, see the 
         [Python/four-disk/Metropolis_disks_box.py](Python/four-disk/Metropolis_disks_box.py) script)
   - [x] Sampling program using Molecular dynamics with pressure estimators (Python, see the 
         [Python/four-disk/molecular_disks_box.py](Python/four-disk/molecular_disks_box.py) script)
   - [x] Sampling program using straight ECMC with pressure estimators (Python, see the 
         [Python/four-disk/ECMC_straight_disks.py](Python/four-disk/ECMC_straight_disks.py) script)
   - [x] Sampling program using reflective ECMC (Python, see the
         [Python/four-disk/ECMC_reflective_disks_box.py](Python/four-disk/ECMC_reflective_disks_box.py) script)
   - [x] Sampling program using forward ECMC (Python, see the
         [Python/four-disk/ECMC_forward_disks_box.py](Python/four-disk/ECMC_forward_disks_box.py) script)
   - [x] Sampling program using Newtonian ECMC (Python, see the
         [Python/four-disk/ECMC_Newtonian_disks_box.py](Python/four-disk/ECMC_Newtonian_disks_box.py) script)
   
- [x] Naive periodic-box programs
   - [x] Sampling program using Metropolis algorithm (Python)
   - [x] Sampling program using Molecular dynamics (Python)
   - [x] Sampling program using straight ECMC with pressure estimators (Python)
   - [x] Sampling program using reflective ECMC (Python)
   - [x] Sampling program using forward ECMC (Python)
   - [x] Sampling program using Newtonian ECMC (Python)

- [ ] State-of-the-art hard-disk programs
   - [ ] Sampling program using straight ECMC with pressure estimators (C++)

- [ ] Analysis
   - [x] Pressure calculation using the fitting formula (Python, see the 
         [Python/four-disk/fitting.py](Python/four-disk/fitting.py) script)
   - [ ] Global Orientational order parameter calculation (Python)
   - [ ] Confidence-interval calculation (Python)

Completed tasks that are marked in the list above are already contained in this version of the repository.

## Installing

All Python scripts can be executed with any Python3 implementation (e.g., standard [CPython](https://www.python.org) or 
[PyPy3](https://www.pypy.org)). The sampling programs do not have any further requirements. Only the script for the 
pressure calculation using the fitting formulas (see [Python/four-disk/fitting.py](Python/four-disk/fitting.py)) relies 
on [NumPy](https://numpy.org) as an external dependency (see [requirements.txt](requirements.txt)). We recommend setting 
up a virtual environment and installing the correct NumPy version by running the following commands (replace `python3` 
by the Python3 interpreter of your choice):

```shell
python3 -m venv historic_venv
source historic_venv/bin/activate
python3 -m pip install -U pip setuptools
python3 -m pip install -r requirements.txt
```

## Authors 

Check the [AUTHORS.md](AUTHORS.md) file to see who participated in this project.

## License

This project is licensed under the GNU General Public License, version 3 (see the [LICENSE](LICENSE) file).

## Contact

If you have questions regarding the HistoricDisks software package, just raise an issue here on GitHub or contact us via 
mail (see the [AUTHORS.md](AUTHORS.md) file). We are happy to help you!

## Citation

If you use (parts of the) the programs in published work, please cite the following reference (see
[[Li2022]](https://arxiv.org/abs/2207.07715) in [References.bib](References.bib)):

Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, and Werner Krauth,\
Hard-disk computer simulations&mdash;a historic perspective,\
arXiv e-prints: 2207.07715 (2022), https://arxiv.org/abs/2207.07715.
