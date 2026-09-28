# Rusty WAM

Rusty WAM is currently a prototype in very early development. Many features have yet to be implemented.

## Overview

Rusty WAM is a 1D Euler equation solver intended to model flow within internal combustion engines, inspired by the open source [OpenWAM]([url](https://github.com/EliDeCo/OpenWAM-Refactored)) repository. 

## Current State
Pipes, junctions, and boundary conditions implimented. Can choose between Roe, RoeM and Muscl + Roe interior methods, as well as Euler, Ssp2, Ssp3, and Ssp43 time integration schemes.

## Validation
Every implimented method is validated against all the tests in the paper it was derived from, explained in detail in the corresponding file in the `validation` directory. The convergence table and Sod's problem output for the 3rd order Muscl + RoeM method is shown below:

| Cells | L1 error | Observed order |
|---|---|---|
| 100 | 9.3723e-4 | - |
| 200 | 1.2060e-4 | 2.96 |
| 400 | 1.5140e-5 | 2.99 |
| 800 | 1.8961e-6 | 3.00 |
| 1600 | 2.3842e-7 | 2.99 |

![Sod's Problem 3rd order Muscl + RoeM](validation\Sods.png)

## Planned Features:
Pipe Solver

- [x] Interior pipe solver for calculating flow
- [x] Junction model and solver
- [x] Boundary conditions or inlets and outputs to the atmosphere
- [ ] Visual editor for constructing, removing, and joining nodes

Minimum Viable Engine Solver
- [ ] 0D Cylinder and Wiebe combustion model
- [ ] Ability to simulate naturally aspirated 4 stroke combustion engine within 5% of empirical data

Long Term Additional Features
- [ ] Options to simulate certain regions in 2D or 3D
- [ ] Support for 3D geometry to influence flow calculations
- [ ] Support for more engine types such as 2 Stroke, Wankel, and Opposed Piston
- [ ] Emissions tracking
- [ ] Full species based combustion for simple fuels like hydrogen and natural gas