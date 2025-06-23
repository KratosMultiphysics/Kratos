---
title: Overview
keywords:
tags: [Overview.md]
sidebar: cosimulation_application
summary:
---
## Overview
This application handles the interaction between two or more solver in a CoSimulation. Generally, the user can control which coupling strategy to be used, the order of solver to be executed, and how the data are transfered between the solvers. Each solver will have their own mesh file and their json file containing the solver's settings
![MPMApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/workshop_2019_tutorials/fsi3.jpg)

A detailed explaination of the json structure can be found [here](../General/JSON_Structure.html)

The following solver wrappers are available:
- [Convection-Diffusion](../Solver_Wrappers/Convection-Diffusion.html)
- [DEM](../Solver_Wrappers/DEM.html)
- [Fluid Dynamics](../Solver_Wrappers/Fluid_Dynamics.html)
- [MPM Dirichlet](../Solver_Wrappers/MPM_Dirichlet.html)
- [MPM Neumann](../Solver_Wrappers/MPM_Neumann.html)
- [PFEM Fluid Dynamics](../Solver_Wrappers/PFEM_Fluid_Dynamics.html)
- [Potential Flow](../Solver_Wrappers/Potential_Flow.html)
- [Rigid Body Solver](../Solver_Wrappers/Rigid_Body_Solver.html)
- [Structural Mechanics](../Solver_Wrappers/Structural_Mechanics.html)
- [SDOF Solver](../Solver_Wrappers/SDOF_Solver.html)
