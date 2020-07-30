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

### CUDA-capable workstations - important note

If you have a CUDA capable workstation, you can manually enable the CUDA support through the CMake file (the npFEM solver can solve bodies in NVIDIA GPUs). The npFEM-GPU branch supports GPUs with compute capability >6.0, due to the use of the ***atomicAdd*** function. For GPUs with lower compute capability, you could implement an alternative, but less efficient, atomicAdd as instructed in **CUDA toolkit documentation**. To check the requirement of the compute capability, go to CUDA samples (comes with the CUDA installation), compile the deviceQuery example, and by running it you will have a detailed list of your GPU specifications.

For optimization purposes (registers), we have used the CUDA-specific ```__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)``` in ```palabos_root_dir/coupledSimulators/npFEM/src_GPU/device_utilities.cu```. The  ```maxThreadsPerBlock ```  specifies the maximum number of threads per block with which the application will ever launch the targeted CUDA-kernel. The  ```minBlocksPerMultiprocessor ``` is optional and specifies the desired minimum number of resident blocks per multiprocessor. If you need to solve RBCs with more than 258 surface vertices (see ```bloodFlowDefoBodies/Case_Studies``` folder), you have to manually tune the **hardcoded** maxThreadsPerBlock from 258 (current value) to the desired one. For platelets, which always have less surface vertices than RBCs (due to size difference), the ```maxThreadsPerBlock ``` readily covers them.

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

## Cellular Blood Flow Simulations

These simulations either start from a given CPs folder (generated through Cell Packing or provided) or from a file (initialPlacing.pos) that stores the position (center of mass) and orientation of the blood cells. The initialPlacing.pos stores first the RBCs as ID (0: RBC, 1: PLT), position and orientation (in degrees) and after the PLTs. Cell Packing is the preferred way to generate complex flow fields with many parameters.

For a quick start just extract the provided CPs_X.tar.gz and perform the corresponding Cellular Blood Flow Simulation.

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
