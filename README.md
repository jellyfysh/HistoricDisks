[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

# HistoricDisks&mdash;Teaser

This repository accompanies the work 
[Hard-disk computer simulations&mdash;a historic perspective](https://arxiv.org/abs/2207.07715). It provides the 
pressure data of the hard-disk model extracted from the literature since 1953, and also the set of high-precision 
pressures that are presented in the work above. Furthermore, this repository contains naive Markov-chain Monte Carlo 
and molecular dynamics implementations for the hard-disk model, several pressure estimators, as well as a 
state-of-the-art implementations of the hard-disk event-chain Monte Carlo algorithm.

## Version 0

In its current version 0, this repository serves as a teaser for the complete version 1 that is expected to be finished 
in September 2022. In order to demonstrate the overall style of the programs, this teaser version contains a naive 
sampling program for four hard disks in a non-periodic box that uses the Metropolis algorithm. This repository will
be updated continuously until all programs are online.

## List of programs
The full repository will contain:

- [ ] Pressure data, equations of states
   - [ ] Pressure data files (CSV)
   - [ ] Equations of states visualization (Python)

- [ ] Four-disk non-periodic-box programs (Python)
   - [x] Sampling program using Metropolis algorithm (Python, see the 
         [Python/four-disk/Metropolis_disks_box.py](Python/four-disk/Metropolis_disks_box.py) script)
   - [ ] Sampling program using Molecular dynamics with pressure estimators (Python)
   - [ ] Sampling program using straight ECMC with pressure estimators (Python)
   - [ ] Sampling program using reflective ECMC (Python)
   - [ ] Sampling program using forward ECMC (Python)
   - [ ] Sampling program using Newtonian ECMC (Python)
   
- [ ] Naive periodic-box programs
   - [ ] Sampling program using Metropolis algorithm (Python)
   - [ ] Sampling program using Molecular dynamics (Python)
   - [ ] Sampling program using straight ECMC with pressure estimators (Python)
   - [ ] Sampling program using reflective ECMC (Python)
   - [ ] Sampling program using forward ECMC (Python)
   - [ ] Sampling program using Newtonian ECMC (Python)

- [ ] State-of-the-art hard-disk programs
   - [ ] Sampling program using straight ECMC with pressure estimators (C++)

- [ ] Analysis
   - [ ] Pressure calculation using the fitting formula (Python)
   - [ ] Global Orientational order parameter calculation (Python)
   - [ ] Confidence-interval calculation (Python)

Completed tasks that are marked in the list above are already contained in this version of the repository.

## Installing

The currently contained Python scripts can be executed with any Python3 implementation without any further requirements.
The docstrings in the Python scripts provide further information on their parameters and their output.

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
