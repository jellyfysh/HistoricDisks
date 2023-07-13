## Straight ECMC simulation of hard-disk system implementing cell-based data structure in C++.

This code is based on a `C++` translation of `Fortran90` code used
in [Bernard2011](https://link.aps.org/doi/10.1103/PhysRevLett.107.155704) (more specifically, `GraphValidateCellECMC.cc`
described in [Li2020](https://doi.org/10.1016/j.cpc.2020.107702)). Unlike the original `Fortran90` program, the pressure
can be calculated during the run using equation (21) in [Li2022](https://doi.org/10.1063/5.0126437). The active disk,
resampled after a fixed time interval (a chain), moves in both x and y directions. The direction is updated sequentially
at each resampling.

Physical parameters are set with command-line arguments. The program is able to generate initial configuration from the
command-line arguments, and reading in initial configuration is also an option.

Multidimensional arrays are managed using `boost::multi_array`, see [here](https://www.boost.org/). All input and output
are in the format of `HDF5`, see [here](https://www.hdfgroup.org/solutions/hdf5/).

### Dependencies

The program is tested with `cmake 3.26.4`, `boost 1.71.0` and [`HighFive 2.7.1`](https://github.com/BlueBrain/HighFive).
The `HighFive`
repository is imported as a submodule of this repository. Please refer to the [readme](HighFive/README.md) of HighFive
for its dependencies.

### Usage

* After cloning the repository, use `git submodule update --init` to initialize the `HighFive` submodule.

* Compile the code in the `build` directory using `cmake`:
  ```console
  cmake . -B build/
  cmake --build build
  ```
  This creates the executable `bin/StraightEcmc`. On any changes to C++ code, only the second line has to be
  run again. Note that this procedure uses the `-O3 -DNDEBUG` compilation flags.
* For a (slower) debug build, use the following cmake commands:
  ```console
  cmake . -B debug_build/ -D CMAKE_BUILD_TYPE=Debug
  cmake --build debug_build
  ```
* Run by `./bin/StraightEcmc Nx Ny eta slant length_factor shape out_file in_file`, where `Nx`, `Ny`, ..., `out_file`
  are command-line arguments. Their meanings are explained in the [Input](#input) section. As an example, a run
  that can be featured in one of the blue points in FIG. 9(b) is created by
  executing `./bin/StraightEcmc 9 8 0.708 0 1 square output.h5 dummy.h5`. However, when producing the results, we used a
  much longer run.

### Conventions

The following conventions is adopted in this program.

* The number of disks is not directly given. As we need to create box and initial configuration that forms a perfect
  crystal, the disk number is expressed as a multiplication of Nx, the number of disks in a row in a triangular lattice,
  and Ny, the number of rows. Ny has to be an even number to comply with the periodic boundary condition. The Nx is
  associated with the length of the box and the Ny is associated with the width of the box when the box allows perfect
  crystal. If the box geometry does not allow perfect crystal to form, only Nx * Ny is relevant.
  ``` 
   ____________________________
    /\  /\  /\  /\  /\  /\  /\  /  | -> 2 rows
   /__\/__\/__\/__\/__\/__\/__\/   | -> direction of a row, 8 disks in a row
  ```
* The system volume is always one. The program supports three box shapes: "crystal", "square", and "rectangle". The shape
  of the system controls not only the aspect ratio of the box but also the generated initial configuration. For a
  "crystal" system, the
  aspect ratio of the box is the aspect ratio of the triangular lattice, namely $(1 : Ny / Nx * \sqrt{3} / 2)$. For a
  "square" system the aspect ratio is $(1 : 1)$, while for a"rectangle" system it is $(1 : \sqrt{3} / 2)$. For a "crystal"
  system, the initial configuration is a perfect crystal, while for the
  other two cases, the initial configuration is almost fully packed in a part of the box and empty in the other part.
* If the number of disks does not allow for a perfect crystal to form, one can set the Nx to the number of disks and
  Ny to one. For this case, the shape of the system has to be set to "square" or "rectangle", or a very slim system
  will be formed.
* The density eta is given in the form of packing fraction. The radius of the disk is derived from eta.
* In the slanted crystalline configuration, the y coordinate of the disks in a row is slightly increased from left to
  right. A n-slant configuration connects the right-most disk of the first row with the
  left-most disk of the (n+1)th row. Here, n has to be an even integer.
* The center of the box is the origin of the coordinate system.
* The pressure is given in the unit of $\beta P V_0 / N$, i.e. the legacy unit.

### Input

The program takes eight command-line parameters. Some of them are used in both initial configuration creation and
simulation, while some of them only participate in initial configuration creation.

The parameters featured in both are:

* `Nx`: number of disks in each row, positive integer.
* `Ny`: number of rows, positive even integer.
* `eta`: density expressed in packing fraction, float number.
* `slant`: the number of slant, not used if the shape is square or rectangle, positive even integer. This parameter is
  only used during initial configuration creation.
* `length_factor`: proportional to the length of run, positive integer. By default, the `length_factor` is set to one, 
  and the program runs for roughly 200000000 collisions, corresponding to roughly 20 seconds on a laptop
  with Intel CORE i7 9th Gen. When setting to a value other than one, the length of the run is multiplied by the value 
  of `length_factor`. It is recommended to do a test run with `length_factor` set to 1 and adjust it according to
  the time of the run. As there is no dependence of length of run on the number of disks, larger systems requires
  larger `length_factor` in general.
* `shape`: string representing the shape of the system, could be either "square", "rectangle", or "crystal".
* `out_file`: the path to the HDF5 file storing output of the run. The detail of its contents is
  described in the [Output](#output) section.
* `in_file` the path to the HDF5 file storing the initial configuration of the run. This file has to possess the initial
  configuration--dataset "config-init".

### Output

All the output of the run is stored in a single HDF5 file. This file contains:

* `config-init`: the initial configuration of the run.
* `config-i`: the ith sampled configuration.
* `count`: the configuration count.
* `pressure_x`: the pressure in x-direction calculated using equation (21).
* `pressure_y`: the pressure in y-direction calculated using equation (21).
* `pressure`: the pressure calculated by averaging `pressure_x` and `pressure_y`.
* `parameters/L`: the box geometry.
* `parameters/N`: the number of disks.
* `parameters/Nx`: Nx as described in section [Conventions](#conventions).
* `parameters/Ny`: Nx as described in section [Conventions](#conventions).
* `parameters/eta`: packing fraction.
* `parameters/shape`: integer indicating box shape. 0 indicates crystal; 1 indicates square; 2 indicates rectangle.
* `parameters/sigma`: disk radius.
* `parameters/slant`: slant as described in section [Input](#input).
* `stats/EPH`: estimated collisions per hour.
* `stats/collisions`: total number of collisions.
* `stats/start_time`: real-world time when the run starts.
* `stats/end_time`: real-world time when the run terminates.