# Thruster Simulation (Plasma and Vacuum) in Julia!

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

Mesh Sims:
    e.g. dx,dy,dz meshes over a discrete space. each spatially connected "node" of simulation is connected in the simulation domain too.
    Just a mesh of data points, w/ a predefined number of neighbors but these can become degenerate or tangled during simulation if the 
    sim    material is moving like in plasma physics.
    
Meshfree Sims:
    interaction of each node with its neighbors. Can use different shapes.
    Can do Lagrangian simulations so that the "node" can move according to the velocity field!


## Particle-In-Cell Julia

Solves the (what?) PDE for charged particles in an e-field by tracking individual particles in a Lagrangian frame.
Eulerian (specification of the flow field) focuses on which specific locations the fluid passes as time passes.
    Eulerian sims employ a fixed mesh.
Lagrangian (specification of the flow field) follows a fluid parcel (def: a constant mass 3D infinitesimal volume of fluid) as it moves through spacetime (to get its pathline)
    Lagranian sims use simul

https://github.com/AndrewChap/Q-ThrusterSimulation/blob/master/README.md

https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.033312
