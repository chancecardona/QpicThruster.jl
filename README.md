# Quantum Particle-in-Cell (Q-PIC) Thruster simulation in Julia

-[![Build Status](https://github.com/chancecardona/QpicThruster.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/chancecardona/QpicThruster.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Background (Plasma Simulation)

Fluid Parcel:
    Constant mass infinitesimal volume of fluid.
    Volume may change (compressible) or be const (isochoric).
    These describe the avg velocity and other properties of fluid particles when averaged over a (length >> the mean free path). 
        This length is the "characteristic length" L_c. Typically = Volume / SurfaceArea.
        Or in rockets it's the combustion chamber's volume / nozzle throat area.
            On deLebal nozzles the cross section of the chamber is bigger than the throat so the L_c is > physical length of chamber.
Knudsen number (Kn) is the ratio of the mean-free-path to the L_c. 
    determines if you should use continuum mechanics (Kn << 1) or statistical mechanics (Kn >~ 1) for a fluid.
        Kn < 0.01: continuum flow (airflow around an aircraft)
        0.01 < Kn < 0.1: Slip flow
        0.1 < Kn < 10: transitional flow
        Kn > 10: free molecular flow (e.g. a dust particle thru the lower atmosphere, or a satellite in the exosphere)
    Depends on the density of the gas (and thus the temperature) due to the mean free path.
    Can use to determine the rarefaction of a flow.
    Can seperate gasses of different molecular masses (e.g. isotopes) by sending the mixed gas thru small holes in a thin wall,
        since the pressure of the gas is inv propertional to its molecular mass.
U_infinity = The freestream speed [L/T units]
    freestream: the air itself that is far upstream of an aerodynamic boody.
    so this is the air before the body can deflect slow down or compress the air.
Flow Velocity:
    flow velocity "u" (continuum), macroscopic velocity (statistical), drift velocity (E&M), describes motion of a continuum of particles.
    u = u(x,t) (flow velocity is a time varying vector field, where at a specific x,t it gives the velocity at that point of spacetime.
    bulk velocity: average of u. = volume flow rate / cross sectional area.

Boltzmann Relationship:
    n_e = n_0 * exp(- (q_e * phi) / (k*T_e))
    where:
        n_e is the electron density at a given potential (phi)
        n_0 is the equilibrium electron density
        q_e is the elementary charge of the electron
        phi is the electric potential
        k is the boltzmann constant
        T_e is the temperatue of the electrons.
    Relates how density of electrons changes, and simplifies it for calculation. 
    This assumes the electrons are much less massive and much faster than the ions and thus move "instantaneously".


Mesh Sims:
    e.g. dx,dy,dz meshes over a discrete space. each spatially connected "node" of simulation is connected in the simulation domain too.
    Just a mesh of data points, w/ a predefined number of neighbors but these can become degenerate or tangled during simulation if the 
    sim material is moving like in plasma physics.
    
Meshfree Sims:
    interaction of each node with its neighbors. Can use different shapes.
    Can do Lagrangian simulations so that the "node" can move according to the velocity field!

Flow Field Specifications:
    Eulerian (specification of the flow field) focuses on which specific locations the fluid passes as time passes.
        Eulerian sims employ a fixed mesh, which discretizes x (and or y, z, etc) and t.
    Lagrangian (specification of the flow field) follows a fluid parcel as it moves through spacetime (to get its pathline)
        Lagranian sims use a mesh representing not spacetime but the fluid-parcels, and instead needs to calculate their x position
        at a specific time.

Particle-In-Cell:
    PIC methods use the langrangian specification to track a fluid parcel in (continuos) phase-space by:
    1. Calculate Plasma parameters and values such as the Charge density. 
       This is done by summing up each "particle" for each "cell" and then distributing each particle's effects to the appropriate mesh point and apply boundaries.
       The cell is the area (if 2D) or the volume (if 3D), etc formed by the mesh-points.
    2. Solve Maxwell's Equations (or an approximation / simplification of them) for the Electric (and Magnetic, etc) potential.
       - Finite Difference Method:
         Do calculation for the field only on the discrete set of mesh-grid points. 
         Can then utilize a solver for the PDE on the mesh points, such as the [Gauss-Seidel method](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method) used to solve Poisson's Equation in 2D for our geometry.
       - Finite Element Method:
         Also discretize space into mesh points, then take the PDE and optimize error of a trial solution (using basis-functions) until convergence.
       - Spectral Method:
         Use the FFT and then again optimize solutions using basis functions (but now the eigenvalue problem is high oerder and defined over the whole domain usually).
    3. Next, we can move the particles by weighing the fields (summed over each mesh point that encapsulates the cell) for each particle to calculate the force.
    4. Can then solve the particle velocity and position by using something like the [Leap-frog method](https://en.wikipedia.org/wiki/Leapfrog_method) (a 2nd order explicit solver) (or the [Boris Method](https://www.particleincell.com/2011/vxb-rotation/)) to compute the position and velocity at the appropriate timesteps after the fields are updated.
       4.5 Boundary conditions are important to keep track to and can range from moving the particle, specular diffraction, deleting the particle, etc).
    5. We can then Apply Boundary Conditions (if not done in the previous step).
    6. Generate Particles (such as the Birdsall approximation of a Maxwellian electron sampling distribution).
    7. And whatever other effects to account for (like collisions).
    
Simulating Collisions:
    Direct Simulation Monte Carlo (DSMC):
        A method for simulating interparticle collisions. This is essentially PIC for neutral gasses.
        Another approach is to solve the fluid conservation equations. (CFD)(https://en.wikipedia.org/wiki/Computational_fluid_dynamics) techniques generally assume a Maxwellian distribution function for the gas molecule velocities (and thus isn't valid for non thermalized fluid flows).
        DSMC doesn't assume any distribution and instead, at each time step, pairs up each particle that's in the same cell and then collides them according to a specified probability.
    Monte Carlo Collisions (MCC):
        If density of the collision target >> density of the collision source, and the collision frequency is low, we can treat the target particles as a "cloud".
        This is much faster but doesn't conserve momentum.
    

## Particle-In-Cell Julia

This repo uses the Particle In Cell method on a fixed mesh (Eulerian specification) to simulate various thrusters and plasma objects.
Oftentimes COMSOL is used with the AC-DC module to model the magnetic field, which can be exported 

### Flow Around Plate
The first, and working is a simple 2D conducting plate (that serves as an absorbing boundary) accounting for electrostatics.
See https://www.particleincell.com/2010/es-pic-method/ for more explanation.

This plotss the electron Density and the E-field across the geometry.
To run this, try:
```
julia src/flow_around_plate.jl
```

TODO: Use [Gridap.jl](https://github.com/gridap/Tutorials/blob/master/docs/src/index.md) to model the mesh and solve the PDE's using FEM.
Can also implement something like [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) for CFD.

### Hall Effect Sim
(Reference: *Particle in Cell Simulation of Stationary Plasma Thruster - Taccogna*).
Code references for this already exist, see:
https://docs.juliahub.com/General/HallThruster/stable/
called as an example in:
```
julia src/hall_thruster/hall_thruster.jl
```


Would be cool to eventually model this in a consistent framework and verify against the model listed above (TODO).

### Q Thruster Sim (TODO)
I am currently in the process of converting this repo into julia.
https://github.com/AndrewChap/Q-ThrusterSimulation/blob/master/README.md
Code Setup (since it doesn't have much of a readme):
- QThrusterMain.cu
  - Defines macro constants (like CPUrun, toggling if particles are created on CPU or GPU, etc)
  - Defines physical constants, like pi
  - Create E&B Fields:
    - Loads `input_deck_2_7.txt`, a COMSOL file.
    - Reads the electric field file's array's rows and col's, and other COMSOL consts from the file.
    - Uses `EB_LOAD.cuh` to load the `complex/b_real.txt` `complex/a_imag.txt` etc field files to arrays.
    - Then starts calculating lagrangian packet parameters like real particles per macroparticle, calculates Cone (E-M cavity shape) geometry and those parameters.
  - Original CPU Calculation:
    - Defined last in `QThrusterMain.cu`. 
    - Original derivation, has the most comments on the physical effects being taken into consideration and the branching simulation strategy, explained in QThrusterMemoryOrganization.pptx
  - CPU:
    - Modern CPU kernel is then defined if not using GPU.
    - Manually creates electrons/positrons as normal, but can see a little clearer than in the GPU, but with optimizations.
  - GPU:
    - Starts managing branches as outlined in the QThrusterMemoryOrganization.pptx, transfers particles to GPU, or can create them on GPU.
    - ThrustCalculationGPU 1-3 are then called.
    - `kernel.cu` is called here as simulation is run as normal using full branching strategy described in powerpoint.
    - GPU memory is then freed.

Process:
- Convert .sln file and .vcxproj files to (cmake)/ makefiles (MakeItSo.jl for the project?)
- Convert ForceCurveMatlabReader.m to jl
- Figure out how to load comsol data into julia --> should already be done by the cuda.
  - `complex/b_imag.txt` and etc are literally just `sprintf`'d and indexed into arrays.
- Need to save data to VTK file, then view it.
  - WriteVTK
  - [VTKView](https://juliapackages.com/p/vtkview) [examples](https://github.com/j-fu/VTKView.jl/blob/master/src/Examples/Examples.jl)
    - Can use this to verify against original cone_open.vtk and cone_closed.vtk files.
- Call kernel and .cuh files from julia if need be. Call QThrusterMain.cu in jl or redefine it.
