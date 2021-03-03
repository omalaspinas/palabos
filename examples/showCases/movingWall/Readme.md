## Waving Plate
This example shows you a rectangular plate that repeatly waving along the rotation axis as y-axis in the fluid domain. 
#### Definition of time
"T timeLB = iT + 1.0;" while "iT" is the time step starts from 0 to maximum iteration time.
#### Behavior of the plate
The rectangular plate will not move along with flow but having a behavior that defined by parameters such as "ampl", "freq" in "struct Param" and time t. 
Those parameters can be manually changed by the user.
By default the angle of this repeatly waving obeys a sinusoidal function in time.
The sequence of utilizing the velocity of the vertices on the plate is that:
1. calculate the angle of the plate from the "getAngle(T t)".
2. have the angular velocity from the "getAngularVelocity(T t)".
3. thus the velocity of every vertex can be obtained from the "class SurfaceVelocity" and then be used in "inamuroIteration()".
4. update of new positions of vertices during iterations will be achieved by looping over every vertex,
 and calculating coordinate z as r*cos(phi), coordinate x as r*sin(phi), while r is the distance between the vertex and the rotation axis.
 The coordinate y doesn't need to be changed.
#### Let rhoBar and j involved in iteration
To perform the immersed boundary method, it is required to have rhoBar and j for the interpolation and force spreading processes.
The code made below operations to fullfill this requirement.
1. when defining lattice, rhobar, and j, there is a pointer.
2. define a "rhoBarJarg" to hold for lattice, rhoBar, and j.
3. let "rhoBarJarg" being integrated by "ExternalRhoJcollideAndStream3D" and "BoxRhoBarJfunctional3D", 
 - The "ExternalRhoJcollideAndStream3D" will make variables inside to do a "collideAndStream". 
 So here it is for "rhoBarJarg"'s collideAndStream during iterations.
 - The "BoxRhoBarJfunctional3D" fundamentally will perform "computeRhoBarJ" for every cell of the domain and will return rhoBar and j.
 So here it is for "rhoBarJarg"'s update of rhoBar and j during iterations.
 
 The "ExternalRhoJcollideAndStream3D" was set as level 0, and "BoxRhoBarJfunctional3D" was set as level 2. 
 Then by "integrateProcessingFunctional", during every iteration those processors with non-negative level will be executed, starting from level 0.

4. after initialization of equilibrium functions by "initializeAtEquilibrium", implement once "BoxRhoBarJfunctional3D" by "applyProcessingFunctional".
This operation will return the rhoBar and j to the "rhoBarJarg".
5. then the "rhoBarJarg" will be ready for "lattice->executeInternalProcessors()" and "inamuroIteration()".
#### Create and place the immersed rectangle plate
The immersed surface can be analytically obtained or loaded from STL file. The idea is to place the triangle set to the desired location,
 then get the information of vertices and areas by the "DEFscaledMesh", and push those information to the defined "vertices" and "areas". Finally define a container for the vertices and areas for the implementation of the "inamuroIteration".
1. define vertices and areas
2. use "constructRectangle" by "param.xSideLB" and so on to analytically create a surface "rectangleTriangleSet", 
instantiateImmersedWallData(vertices, areas, container);
3. define the "mount", then use translate to remove the plate to the desired location by "rectangleTriangleSet.translate(mount);"
 We should be careful with the subtracted half "param.ySideLB" in the defining of the "mount",
 that is because the plate itself has length, we want to reach to the correct position by adding "param.mountPointLB" after the subtraction of y direction.
4. use "rectangleDef" from the placed "rectangleTriangleSet" to get the information of the plate, and assign values into vertices and areas.
5. define a "container" that holds for the vertices and areas.
6. then the "container" will be ready for "inamuroIteration".

#### The main part of one iteration
1. lattice->executeInternalProcessors(); 
 The "rhoBarJarg" will do the "collideAndStream()" and then return the rhoBar and j.
2. update the postions of the immersed plate.
3. instantiate the immersed data with the "container".
4. perform the immersed boundary iterations.

