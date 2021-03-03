# Use TRT dynamics
In this example we simulate the lid driven cavity selecting different variants of the TRTdynamics.
The BGK dynamics can be considered a sub-set of the TRT, when the symmetric (omega+) and antisymmetric (omega-)
collision frequencies are equal. In the general case, when using the TRT collision model,
we are not interested in manually choosing both omega+ and omega-, but we want to set omega+ that
is related to viscosity and choose omega- do get the desired behaviour of the collision model.
The parametric behavior of the TRT is controlled by the magic parameter Lambda, that can be
set in the following way:
````dynamics->setParameter(dynamicParams::magicParameter,magic);````

In Palabos, there are different variants of the TRTdynamics. In this example we show how you can use them.
Check how function 

``IsoThermalBulkDynamics<T, DESCRIPTOR> *
  selectDynamics(const CollisionType collision_type, double tau_plus, double magic);``
  
selects different type of TRT and BGK dynamics. If you want to try this example with a different
dynamics, you can set `collision_model` (in main()) equal to the one of the value defined by the `CollisionType`
enum:

``enum class CollisionType {BGK,TRT,BGKma1,TRTma1};``

Ma1TRTdynamics referes to the TRTdynamics implemented using the linearized equilibrium, useful if you are 
interested in Stokes flows.

Feel free to try also other dynamics of Palabos.

## Lid driven cavity in 2D

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