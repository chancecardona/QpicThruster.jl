# Quantum Particle-in-Cell (Q-PIC) Thruster simulation in Julia

-[![Build Status](https://github.com/chancecardona/QpicThruster.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/chancecardona/QpicThruster.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Particle-In-Cell Julia

This repo uses the Particle In Cell method on a fixed mesh (Eulerian specification) to simulate various thrusters and plasma objects.
See [BACKGROUND](BACKGROUND.md) for more info.
Oftentimes COMSOL is used with the AC-DC module to model the magnetic field, which can be exported and loaded like in the QThruster repo listed below.

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
Code references for this already exist, see: https://docs.juliahub.com/General/HallThruster/stable/
called as an example in:
```
julia src/hall_thruster/hall_thruster.jl
```

Would be cool to eventually model this in a consistent framework and verify against the model listed above.

### Q Thruster Sim
This repo (credit Andrew Chap): 
https://github.com/AndrewChap/Q-ThrusterSimulation/blob/master/README.md
has been converted to build on linux with cuda and use Julia instead of matlab:
https://github.com/chancecardona/Q-ThrusterSimulation/tree/linux-cmake

In the future it would be nice to simulate the dynamical QFT vacuum (including more than electrodynamics even hypothetically using lattice gauge theory) using Trixi.jl or another more coherent framework, and possibly take advantage of acclerated hardware and algorithms next. 
For more on modelling the vacuum as a virtual fermionic plasma see:
- Brady, White, March: 'Dynamics of the Vacuum and Casimir Analogs to the Hydrogen Atom'
- Urban: 'The quantum vacuum as the origin of the speed of light'.
- Bilson-Thompson, Leinweber, Williams: 'Highly-improved lattice eld-strength tensor' for full lattice calculation potentailly.

### Simulating the Casimir Effect
- See "Worldlines" for Casimir force approximation techniques, such as White: 'Worldline numerics applied to custom Casimir geometry generated unanticipated intersection with Alcubierre warp metric'.
