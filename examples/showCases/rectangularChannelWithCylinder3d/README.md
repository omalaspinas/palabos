## Palabos Example : rectangularChannelWithCylinder3d

This example simulates a rectangular channel (walls everywhere except for the inlet/outlet) with a cylinder as an obstacle at the middle of the cross-section.

This case study combines essentially the rectangularChannel3d & cylinder3d examples, and its purpose is to investigate the impact of the walls on the formation of the Von Karman Vortex Street. If the walls are too close to the cylinder then the vortex street is suppressed.

The compilation of this example follows exactly the same steps as any other Palabos example.

Keep in mind that the command line arguments referring to the geometry (channel & cylinder) are in physical units, and subsequently they are automatically converted to dimensionless quantities using the cylinder diameter as the reference length scale. The resolution refers to the discretization of the reference length scale (diameter), e.g. a resolution 20 means that the cylinder diameter is discretized with 20 lattice sites.

Palabos output (e.g. gifs, VTKs) is in dimensionless units. The conversion to physical units is simply done by multiplying the lengths with the reference length scale and the velocities with the reference velocity (e.g. the one used to compute the Reynolds number). Furthermore, given a reference length and velocity, one can readily deduct a reference time scale.
