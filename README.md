# Rusty WAM

Rusty WAM is currently a prototype in very early development. Many features have yet to be implemented.

## Overview

Rusty WAM is a 1D Euler equation solver intended to model flow within internal combustion engines, inspired by the open source [OpenWAM]([url](https://github.com/EliDeCo/OpenWAM-Refactored)) repository. Features for the minimum viable product include:

- [x] Interior pipe solver for calculating flow
- [ ] Junction model and solver
- [ ] 0D Cylinder and Wiebe combustion model
- [ ] Boundary conditions or inlets and outputs to the atmosphere
- [ ] Visual editor for constructing, removing, and joining nodes to represent an internal combustion engine
- [ ] Ability to simulate naturally aspirated 4 stroke combustion engine within 5% of empirical data

Future Planned Features Include
- [ ] Options to simulate certain regions in 2D or 3D
- [ ] Support for 3D geometry to influence flow calculations
- [ ] Support for more engine types such as 2 Stroke, Wankel, and Opposed Piston
- [ ] Emissions tracking
- [ ] Full species based combustion for simple fuels like hydrogen and natural gas
