# Instructions for Palabos-npFEM

For more information on the library (numerics & parameters), you can consult the publications below:

* https://arxiv.org/abs/1903.06479
* https://arxiv.org/abs/1911.03062

npFEM is a heavily modified version of ShapeOp (https://www.shapeop.org/). The different naming originates from the fact that the modifications make the solver more like an FEM solver instead of a computer graphics tool. In more details:

1. We have changed radically the original ShapeOp solver. Our approach follows the paper of Liu et al. (https://arxiv.org/abs/1604.07378), thus the solver uses a Quasi-Newton optimization technique, and not the standard local-global alternating approach of ShapeOp.
2. In Computer Graphics, the solvers are approximating Newton's equations, reducing the computational cost. Our solver is not approximating Newton's equations, but we provide a converged solution, focusing on accuracy and physically correct states.
3. We provide a CUDA implementation of the solver.
4. From ShapeOp, we are using mainly the code structure.

## Compilation

Palabos-npFEM has been tested in both UNIX and Windows (provided CMake File):

1. ```mkdir build```
2. ```cd build```
3. ```cmake ..```
4. ```make -j && make -j gpu``` (if ENABLE_CUDA ON)

In Windows, step 4 is executed by opening the Visual Studio solution file (.sln) in the build folder, and building the project/solution using the top **Build** tab. After finishing compilation in /bloodFlowDefoBodies/Release(or)Debug you can find the executable. Move it outside this folder and execute the application with the help of the xml files.

## Perform Cell Packing

Cell Packing is the software that randomly initializes RBCs & PLTs in the flow field, and after some iterations it resolves all the collisions/interpenetrations. The resolved positions per blood cell are stored in the folder CPs. Open/Explore the cellPacking_params.xml to decide the geometry, hematocrit and many other parameters.

**CPU-version**:
```
mpirun -n X ./bloodFlowDefoBodies cellPacking_params.xml
```

**GPU-version**:
```
mpirun -n X ./bloodFlowDefoBodies_gpu cellPacking_params.xml NumberOfNodes NumberOfGPUsPerNode
```
NumberOfNodes: 1 for a workstation, varies in a cluster
NumberOfGPUsPerNode: number of GPUs per node

In folder tmp you can find the vtks of RBCs & PLTs, along with profiling and other output.

After having the initial positions of the blood cells (stored in the CPs folder), we can perform simulations with a proper flow field.

## Cellular Blood Flow Simuations

These simulations either start from a given CPs folder (generated through Cell Packing or provided) or from a file (initialPlacing.pos) that stores the position (center of mass) and orientation of the blood cells. The initialPlacing.pos stores first the RBCs as ID (0: RBC, 1: PLT), position and orientation (in degrees) and after the PLTs. Cell Packing is the prefered way to generate complex flow fields with many parameters.

For a quick start just extract the provided CPs_X.tar.gz and perform the corresponding Cellular Blood Flow Simuation.

***It is very important to understand that without the right CPs folder (generated through Cell Packing or provided) or initialPlacing.pos file, the Cellular Blood Flow Simuation is not going to start (it crashes).***

**CPU-version**:
```
mpirun -n X ./bloodFlowDefoBodies shear_params/poiseuille_params.xml
```

**GPU-version**:
```
mpirun -n X ./bloodFlowDefoBodies_gpu shear_params/poiseuille_params.xml NumberOfNodes NumberOfGPUsPerNode
```

## Case study: Collision at an obstacle

We provide one case study that uses the initialPlacing.pos. It can be executed as normally with the obstacle_params.xml and it involves 1 RBC colliding with a cubic obstacle.
