# Instructions for npFEM (ShapeOp) with Rhino3D-Grasshopper

This is a stand-alone version of the library, in the sense that it is completely detached from Palabos.

Whatever presented here is intended for Windows only because Rhino3D & Grasshopper (https://www.rhino3d.com/) are available only for Windows (and Mac recently). Also, we are using Rhino version 5 and Grasshopper (GH) installed separately (in latest Rhino, GH is integrated in it). The last choice is because we want to keep this code compatible with the ShapeOp project (https://www.shapeop.org/).

npFEM is a heavily modified version of ShapeOp. The different naming originates from the fact that the modifications make the solver more like an FEM solver instead of a computer graphics tool. In more details:

1. We have changed radically the original ShapeOp solver. Our approach follows the paper of Liu et al. (https://arxiv.org/abs/1604.07378), thus the solver uses a Quasi-Newton optimization technique, and not the standard local-global alternating approach of ShapeOp.
2. In Computer Graphics, the solvers are approximating Newton's equations, reducing the computational cost. Our solver is not approximating Newton's equations, but we provide a converged solution, focusing on accuracy and physically correct states.
3. We provide a CUDA implementation of the solver.
4. From ShapeOp, we are using mainly the code structure.

For more information about ShapeOp, and how it can be integrated in Rhino-GH, you can consult the links below:

* https://www.shapeop.org/
* https://doi.org/10.1007/978-3-319-24208-8_42

Essentially, we are going to compile the library into a dynamic library (dll), and use Rhino-GH for simple tasks and visualization. By simple tasks, we mean one RBC/PLT only, as this is intended for building the materials and then feed them into the Palabos-npFEM library for more intricate simulations.

This part of the project is inspired by the project below:

* https://github.com/AndersDeleuran/ShapeOpGHPython

What we provide in the npFEM_RhinoGH folder adopts ideas and code from the project above, and it is tailored for blood cells.

## Compilation

npFEM has been tested in both UNIX and Windows. Before following the steps below, make sure to set the enviroment variable PALABOS_ROOT as the root folder of palabos.

Windows only (CMake & Visual Studio):

1. mkdir build (current folder)
2. cd build
3. cmake ..
4. make (through Visual studio)

This results in a dynamic library inside the build folder.

* Follow the instructions at ShapeOp page on how to install the Rhino-GH environment.
* Copy the produced dll in the libraries folder of Grasshopper.
* Open Rhino-GH and the Grasshopper file provided in the npFEM_RhinoGH folder.

## Rhino-GH

* Go to npFEM_RhinoGH folder
* Open the npFEM_RhinoGH.3dm Rhino File (this will open Rhino3D)
* In Rhino's command-line type: Grasshopper (this will open GH)
* In GH's environment open npFEM_RhinoGH.gh file
* The rest are similar to what is described in this project: https://github.com/AndersDeleuran/ShapeOpGHPython
