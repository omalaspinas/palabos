# Lid driven cavity in 2D

This example is the simulation of the famous numerical benchamrk of the two dimensional lid-driven cavity.

![Lid driven cavity geometry. Side  walls have zero velocity.](figs/cavity2d.svg)

The top lid is driven with a velocity $`\vec u=(U, 0)`$ while
all the other boundaries have velocity $`\vec u=\vec 0`$.

## Concepts use in this example

1. Creation of a simluation domain (`MultiBlockLattice2D`) with the BGK collision model (`BGKdynamics()`).
2. Instatiation of finitie-difference regularized boundary conditions `createInterpBoundaryCondition2D` (Skordos-like).
3. Setting on-lattice velocity boundaries on all the domain boundaries (`setVelocityConditionOnBlockBoundaries()`).
4. Setting a constant velocity on a boundary (`setBoundaryVelocity()`).
5. Initializing populations with constant density and velocity (`initializeAtEquilibrium()`).
6. Use of `Box2D` objects to handle simulation subdomains.
7. Output `gif` images on simulation slices (`writeScaledGif()`).
8. Ouput `vtk` files (`ParallelVtkImageOutput2D`, `writeData()`).
9. Coputation of velocity, norm and vorticity (`computeVelocity()`, `computeVelocityNorm()`, `computeVorticity()`).
10. Use of units helper functions `IncomprFlowParam`.
11. Timestep progression `collideAndStream()`.
12. Computation of average quantities `getStoredAverageEnergy()` and `getStoredAverageDensity()`.
13. Time measurement `global::timer()`.