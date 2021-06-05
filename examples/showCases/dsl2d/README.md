# Stability of collision models: Double shear layer

## Introduction
The goal of this exercise is to understand how different collision models behave in the low viscosity regime, and for under-resolved conditions. Particular attention will be paid to particular choices of relaxation parameters, since their values are almost never discussed in the literature.

With this idea in mind, a 2D double shear layer will be simulated in under-resolved conditions since it allows to highlight dispersion and stability issues through very intuitive reasonings [1,2]. This test is carried out in a fully periodic domain without boundary conditions, and it consists in two shear layers evolving over time. By imposing a given (transverse velocity) perturbation at the initial time, these layers roll-up and form two counter-rotating vortices. Interestingly, when the Reynolds number is increased, while keeping the mesh resolution constant, one can observe the formation of secondary spurious vortices. The latter help us vizualizing dispersion issues that usually lead to the simulation blow-up, e.g., when the Reynolds or Mach numbers are increased, or if the mesh resolution if further decreased.   

Your mission is to try out different collisions models implemented in Palabos and observe the impact of the choice of model on the numerical results. In particular, you will quantify the impact of two aspects of the numerical model, namely: (1) the moment space in which the collision step is performed, and (2) the values of high-order and bulk-viscosity-related relaxation parameters. The most common collision models are considered in this exercice, namely, models based on moment spaces (raw (RM), Hermite (HM), central (CM), central Hermite (CHM) and cumulant (K)), Gauss-Hermite formalism (GH), and recursive regularization (RR). It is worth noting that for the D2Q9 lattice, HM- and GH-LBMs are mathematically equivalent.

## Guide
### Exercice 1: Test the compilation and run the code in parallel
Open a terminal in the directory `dsl2d` (for 2D double shear layer) and type the following commands to compile the code with cmake:
```
    cd build
    cmake .. && make -j
    cd ..
```
At this point you should have an executable in the current directory named `dsl2d`. Then run the code in parallel via
```
    mpirun -np 2 ./dsl2d config.xml 10000 128 0.1
```
to run the simulation for Re=10'000, Ma=0.1, N=128 points in each direction, with the collision model and other parameters being defined in `config.xml`. All parameters are provided in either the .xml file and the command line arguments so that the user need not to recompile the code each time parameters are modified. 

The performance of each collision model can be evaluated through
```
    time mpirun -np 2 ./dsl2d config.xml 10000 128 0.1
```
However, to do proper performance measurements, do not forget to comment in `dsl2d.cpp`
```
    ///// Initial state is outputed
    writeVTK(lattice, param, 0);    
    writeGif(lattice, 0);
```
and to fix, e.g., 
```
<!-- .vti output frequency (in terms of caracteristic time) -->
    <vtkT> 100 </vtkT>
```
in the .xml file so that the performance will not be biased by the time required to output data.

### Explanations: Code structure
In the file `utility_dsl2d_param.h`, you will find the class `Param` which is used to read data from the .xml file, and the command line arguments. All parameters are then written in a log file through the function `writeLogFile()`, so that the user can easily reproduce previous runs. It is worth noting that the computation of the time step is based on the physical speed of sound (acoustic scaling) that is given in the .xml file. Doing so, one can impose the proper speed of sound even with standard (isothermal) LBMs.

As far as `dsl2d.cpp` is concerned, `getDynamics()` is used to choose the collision model via the name entered in the .xml file. `DoubleShearLayerInitialVelocityField` initializes macroscopic variables (density and velocity components). `simulationSetup()` sets periodic boundaries, initialize populations at equilibrium via density and velocity fields. `writeGif()` and `writeVTK()` are used to output data either in a .gif or .vti format.

At the end of the day, only the values of dynName, hoOmega, bulk (.xml file) and the Reynolds number (argument of the command line) will be of interest in this exercice.

### Exercice 2: Improving stability without modifying relaxation parameters
Run the code with the standard BGK collision model (BGK_Ma2) for Re=10'000, Ma=0.1 and N=128. The simulation remains stable but spurious vortices appear as the shear layers become thinner and thinner. As this collision model is more dispersive than dissipative [3], small oscillations appear close to strong gradients and in under-resolved conditions. 

Try to change the collision model while keeping a single relaxation time (hoOmega=SRT) and without modifying the bulk viscosity (bulk=false). Only the cumulant (K) and recursive regularized (RR) approaches should remain stable, even though, the K-LBM shows severe dispersion issues. The latter issues are even worse when increasing the Reynolds number to 1'000'000.

### Exercice 3: Improving stability and accuracy by increasing the bulk viscosity
Keep Re=1'000'000, and try all collision models with bulk=true. Doing so, the relaxation frequency related to the trace of second-order moments is fixed to 1, hence imposing nu_bulk = cs2/2 = 1/6 in lattice units. 

Hermite-based formulations (HM, GH and CHM) are now stable but show dispersion issues. Drastically increasing the bulk viscosity definitely reduces dispersion issues observed with the K-LBM, even though, they remain visible. As previously observed, the RR-LBM is the only model that does not show dispersion issues.

### Exercice 4: Further improving stability and accuracy by regularizing high-order moments
In the literature, it is common to regularize high-order moments, i.e., to fix the corresponding relaxation frequency to 1. From the physical viewpoint, this amounts to filtering out non-equilibrium contributions of high-order moments. Mathematically speaking, the rank of the system is decreased, from 9 to 6 for D2Q9-LBMs, since populations only depend on rho, ux, uy, Pixx, Pixy and Piyy. 

To quantify the impact of the regularization of high-order moments, keep Re=1'000'000 and bulk=true, and this time, fix hoOmega=REG in the config file. Now, all collision models lead to stable simulations. While RM, HM, K and RR show dispersion issues, CM and CHM lead to the most accurate results. More precisely, standard and Hermite formalisms are now equivalent (since omega_bulk=omega4) which means that RM is equivalent to HM whereas CM is equivalent to CHM. Interestingly, results obtained with RR are worse for this particular configuration, and this is explained by the fact that imposing hoOmega=REG boils down to discarding the recursive computation of high-order moments, i.e., RR-REG and HM-REG are equivalent. Eventually, CHM-REG recovers the results that were obtained with RR-SRT. In-depth discussions on equivalencies between collision models can be found in Refs [4,5].      

### Exercice 5: Improving stability and accuracy without sacrifying acoustic waves
Even if increasing the bulk viscosity helps stabilizing simulations, it is also done at the cost of over-damped acoustic waves. To understand which collision models can be used for computational aeroacoustic simulations in the low-viscosity regime, impose bulk=false while keeping hoOmega=REG.

Interestingly, all collision models are stable but the CM approach. Once again, RR-REG boils down to HM-REG, and CHM-REG (RR-SRT) shows the best results, closely followed by the cumulant approach.

## Take home ideas
For pure aerodynamics, regularizing high-order moments and increasing the bulk viscosity are good ways to improve the stability of most collision models. Only the RR approach should be used in an SRT formalism. In case of computational aeroacoustic simulations, K-REG, CHM-REG and RR-SRT should be preferred. 

Eventually, it is worth noting that in all these tests, the collision models were used without accounting for turbulence modeling. In practice, collision models should be systematically combined with subgrid scale models when simulating industrial configurations, so that the mesh density remains reasonable [6].

## Bibliography
[1] Brown, D. L., & Minion, M. L. (1995). Performance of under-resolved two-dimensional incompressible flow simulations. Journal of Computational Physics, 122(1). http://www.sciencedirect.com/science/article/pii/S0021999185712053

[2] Minion, M. L., & Brown, D. L. (1997). Performance of under-resolved two-dimensional incompressible flow simulations, II. Journal of Computational Physics, 138(2), 734-765. http://www.sciencedirect.com/science/article/pii/S0021999197958435

[3] Marié, S., Ricot, D., & Sagaut, P. (2009). Comparison between lattice Boltzmann method and Navier–Stokes high order schemes for computational aeroacoustics. Journal of Computational Physics, 228(4), 1056-1070. http://www.sciencedirect.com/science/article/pii/S002199910800538X

[4] Coreixas, C., Chopard, B., & Latt, J. (2019). Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. Physical Review E, 100(3), 033305. https://www.researchgate.net/publication/335727817_Comprehensive_comparison_of_collision_models_in_the_lattice_Boltzmann_framework_Theoretical_investigations

[5] Coreixas, C., Wissocq, G., Chopard, B., & Latt, J. (2020). Impact of collision models on the physical properties and the stability of lattice Boltzmann methods. Philosophical Transactions of the Royal Society A, 378(2175), 20190397. https://www.researchgate.net/publication/339136823_Impact_of_collision_models_on_the_physical_properties_and_the_stability_of_lattice_Boltzmann_methods

[6] Coreixas, C. (2021). On the use of LBM in the industry: Subgrid scale model, wall law and grid refinement. The First International Workshop on Lattice Boltzmann for Wind Energy. https://www.researchgate.net/publication/349681779_On_the_use_of_LBM_in_the_industry_Subgrid_scale_model_wall_law_and_grid_refinement
