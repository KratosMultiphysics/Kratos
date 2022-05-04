---
title: Solving strategies
keywords: 
tags: [Solving-strategies.md]
sidebar: kratos_for_developers
summary: 
---

# Basic structure of the Python layer

## 1. Introduction
This tutorial explains the top-level python layer of Kratos and how to use it. It is based on the previous tutorial of the high-rise building (structural case). The final files can be retrieved [here](https://github.com/KratosMultiphysics/Documentation/tree/master/Workshops_files/Kratos_Workshop_2019/Sources/4_solving_strategies).

## 2. The AnalysisStage
The AnalysisStage has replaced what was previously done in the MainKratos-Script: Plugging together and combining the individual components that are necessary for performing a simulation in a particular application.

This means that **this is the place to do the user-scripting**, which was previously done in the MainScript. This is achieved by deriving a user-defined AnalysisStage based on the AnalysisStage from the application that is used,  e.g. `StructuralMechanicsAnalysis` or `FluidDynamicsAnalysis`.

The big advantage is that the only the functions, which are needed for the user-scripting, will be overridden. This way the if something changes in the baseclass, the changes are automatically used.

More information can be found in the [description of the AnalysisStage](Common-Python-Interface-of-Applications-for-Users#analysisstage).

**The tasks for this part of the tutorial are:**
* Modifying the thickness of the structure
* Printing the reaction forces on the fixed end of the structure.
* Creating point-load conditions, which are used later to change the way the loading is applied on the structure.

Start by checking the functions of the [AnalysisStage](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/analysis_stage.py) for suitable places to perform the above tasks.

## 3. The Processes
For a description of what a process is and how we are going to use it to write a custom process, please read [this entry in the Wiki](Using-processes-to-customize-a-simulation) first.

As opposed to defining a custom AnalysisStage, another way of printing the reaction-forces is to define a custom process. The link above contains a template for a process defined in python.

**Tasks:**
* use the AnalysisStage from above to change the loading from line-loads to point-loads using the following process: `assign_vector_variable_process`
* create a custom process for printing the reactions and add it to the simulation.

## 4. The PythonSolver
The [`PythonSolver`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/python_solver.py) is used inside the AnalysisStage, it is responsible for solving the "physics" of the problem. A detailed list of functionalities can be found [in the Wiki](Common-Python-Interface-of-Applications-for-Users#pythonsolver).

Writing a custom solver is an advanced usecase, which is not covered in this tutorial. An example can be found [in the Wiki](Implementing-thermal-solver).

This part of the tutorial aims to modify the settings of the solver, which can be found in `ProjectParametrs.json` under "solver_settings".

_Helpful advice_: using "wrong" settings on purpose will print hte admissible settings that are available.