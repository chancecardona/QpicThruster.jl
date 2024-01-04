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

## Particle-In-Cell Julia

This repo uses the Particle In Cell method on a fixed mesh (Eulerian specification) to simulate various thrusters and plasma objects.

### Flow Around Plate
The first, and working is a simple conducting plate (that serves as an absorbing boundary) accounting for electrostatics.
See https://www.particleincell.com/2010/es-pic-method/ for more explanation.

This plotss the electron Density and the E-field across the geometry.
To run this, try:
```
julia src/flow_around_plate.jl
```

### Hall Effect Sim
Todo (Particle in Cell Simulation of Stationary Plasma Thruster - Taccogna).
https://docs.juliahub.com/General/HallThruster/stable/

### Q Thruster Sim
Todo https://github.com/AndrewChap/Q-ThrusterSimulation/blob/master/README.md
