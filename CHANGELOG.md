# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [unreleased]

### Added

* Using sitmo for PRNG generation to have parallel PRNG concistency.
* Data processors for 2D mutlri grid and multi block lattices.

### Removed

* Removed the use of `rand()` and `plbRandom()`.

### Changed

* Changed regression testing to account for sitmo PRNG:
* Changed 2D VTK to accept Box2D.

### Fixed

* Bugs in 2D multi grid refinement.

## [2.3.0]

### Added

* Added example of a 3D flow past a cylinder in a rectangular channel.
* Added example of 3D porous media with partial bounce back.
* Added partial bounce back boundary condition.
* Added 3D collision models: Raw (RM), Hermite (HM), Central (CM), Central Hermite (CHM) and Cumulant (K),  Gauss-Hermite formalism (GH), weighted version of populations in the HM formalism,  Regularized approaches: standard (REG-HM) and recursive (SRT-RR).
* Added dsl2d and dsl3d examples.

### Removed

### Changed

* Made TRTdynamics more consistent with original paper.
* Added a parameter to modify the bulk viscosity for 2D and 3D formulations of the RM, HM, CM, CHM, GH collision models.

### Fixed

* Fixed error in the He forcing scheme implementation.
* Fixed minor error in CI causing an exmaple not to be compiled.

## [2.2.1]

### Added

* The D2Q9 version of collision models: Raw (RM), Hermite (HM), Central (CM), Central Hermite (CHM) and Cumulant (K),  Gauss-Hermite formalism (GH), weighted version of populations in the HM formalism,  Regularized approaches: standard (REG-HM) and recursive (SRT-RR).
* Added test case: `isothermCoVo2D`
* Added interface to use magic parameter in TRTdynamics.

### Removed

### Changed

* Renames sMinus to omegaMinus for consistency in the notations in TRTdynamics.

### Fixed

## [2.2.0]

### Added

  * Added custom communicator for possible Cuda additions.
  * Added coupledSimulators folder for external couplings
  * Added npFEM library in coupledSimulators folder
  * Added bloodFlowDefoBodies example in examples/showCases for the simuation of cellular blood flow using npFEM (in coupledSimulators)
  * Added support for HDF5 support.
  * Added cylinder3d example.
  * Added multi block generation from MultiBlock in 2D.
  * Added `ccache` support for CI.
  * Added references for off-lattice BCs.
  * Added `FilippovaHaenelLocalModel3D` that implements the original Filippova-Haenel boundary condition.
  * Added `MeiLuoShyyModel3D` as a replacement fot the `FilippovaHaenelModel3D` off-lattice boundary condition which was wrongly named FH as it was the MLS version that was implemented. **This is a breaking change**.

### Removed

  * Removed `scons` compilation system to keep only `cmake`.
  * Removed the `FilippovaHaenelModel3D` for off-lattice boundary conditions. **This is a breaking change**. To continue using that former FH BC one should switch to `MeiLuoShyyModel3D`.

### Changed

  * Changed CMake system: compile all examples from palabos root directory
  * Changed the Filippova-Haenel BC to be consistent with the paper. **There is a breaking change here.**

### Fixed

  * Fixed bug in the `externalflowAroundObstacle` showCase.
  * Fixed bug with TRT dynamics.
  * Fixed bug with timer in the single threaded compilation.

## [2.1.0]

### Added

### Removed

  * Removed useless SMP compilation flag.

### Changed

  * Updated `scons` to 3.1.1.
  * Compilation for continuous integration made faster.

### Fixed



