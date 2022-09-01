////////////////////////////////
// Settling Driven Convection //
////////////////////////////////

This example describes a Rayleigh-Taylor problem initiated by the settling of particles.

Assuming that particles are coupled with the fluid (small particles, high number ==> no inertia), a single-phase model is implemented. In that case, particles are considered as a continuum phase described by an advection-diffusion-sedimentation equation.
The model presented here is implemented through a hybrid algorithm where the fluid motion is solved using the Lattice Boltzmann Method while the advection-diffusion-sedimentation of particles is solved by a low diffusive finite difference scheme called Weighted Essentially Non-Oscillatory.

The initial goal of this program was to reproduce some experiments performed in laboratory in order to characterise the instability growing at the interface between two layers.
The configuration is a box divided in two layers initially separated by a barrier. The lower layer is a mixture of fresh water and sugar while the upper layer is a mixture of fresh water and particles (glass beads). For the sugar species, a classical 1st order finite difference scheme is used.
The initial particle volume fraction in the upper layer is taken in order to ensure that the bulk density density is less than the lower layer. The system being initially stable, the settling of particle will initiate the growth of an unstable layer underneath the initial interface. Then, any perturbation of this layer will trigger an instability which will grow as downward moving columns of particles.

In order to compile the code, type "cd build" then "cmake .."
Once the cmake is successful, type "make" in order to conclude.
Go back to the main directory ("cd ..") then run the program with "mpirun -np X settlingDrivenConvection" X being the number of processors you want to use.
The different outputs (.gif images of the density field, particle field, fluid velocity field and vorticity, .vtk for paraview) are saved in the /tmp folder.

In the example, the spatial resolution is 250 is order to allow the program to run on a classical computer. Increasing the resolution will provide more precision but will require more resources.

For more precisions, please refer to "Lemus J, Fries A, Jarvis PA, Bonadonna C, Chopard B and Lätt J (2021) Modelling Settling-Driven Gravitational Instabilities at the Base of Volcanic Clouds Using the Lattice Boltzmann Method. Front. Earth Sci. 9:713175. doi: 10.3389/feart.2021.713175"
