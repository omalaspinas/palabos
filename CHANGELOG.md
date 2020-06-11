# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [unreleased]

### Added

  * Added custom communicator for possible Cuda additions.
  * Added coupledSimulators folder for external couplings
  * Added npFEM library in coupledSimulators folder
  * Added bloodFlowDefoBodies example in examples/showCases for the simuation of cellular blood flow using npFEM (in coupledSimulators)

### Removed

### Changed

  * Changed CMake system: compile all examples from palabos root directory

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



