---
title: Common Python Interface of Applications for Users
keywords: 
tags: [Common-Python-Interface-of-Applications-for-Users.md]
sidebar: kratos_for_developers
summary: 
---

# Overview
1. [Introduction](Common-Python-Interface-of-Applications-for-Users#introduction)
2. [AnalysisStage](Common-Python-Interface-of-Applications-for-Users#analysisstage)
    1. [AnalysisStage: Overview](Common-Python-Interface-of-Applications-for-Users#analysisstage-overview)
    2. [AnalysisStage: Responsibilities and provided Functionalities](Common-Python-Interface-of-Applications-for-Users#analysisstage-responsibilities-and-provided-functionalities)
    3. [AnalysisStage: Usage](Common-Python-Interface-of-Applications-for-Users#analysisstage-usage)
3. [PythonSolver](Common-Python-Interface-of-Applications-for-Users#pythonsolver)
    1. [PythonSolver: Overview](Common-Python-Interface-of-Applications-for-Users#pythonsolver-overview)
    2. [PythonSolver: Responsibilities and provided Functionalities](Common-Python-Interface-of-Applications-for-Users#pythonsolver-responsibilities-and-provided-functionalities)
    3. [PythonSolver: Usage](Common-Python-Interface-of-Applications-for-Users#pythonsolver-usage)
4. [Future Outlook](Common-Python-Interface-of-Applications-for-Users#outlook-kratos-project-multi-stage-simulation)

# Introduction
Solving a problem with Kratos is divided into two Python-objects : The **AnalysisStage** and the **PythonSolver**. 

The `PythonSolver` is responsible for everything related to the physics of the problem (e.g. how to setup the system of equations), whereas the `AnalysisStage` is related to everything that is not related to the physics (e.g. when and what output to write).

This means that coupling of physics is done on solver level. A coupled solver would contain the solvers of the involved physics as well as e.g. coupling logic, data exchange, ... .
The coupling of things not related to physics is done in the `AnalysisStage`, e.g. combined output.

The `PythonSolver` is a member of the `AnalysisStage`, therefore the `AnalysisStage` can be seen as "outer" layer and the `PythonSolver` as "inner" layer.

# AnalysisStage
## AnalysisStage: Overview
The baseclass of the `AnalysisStage` is located in the KratosCore ([here](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/analysis_stage.py)). It provides a set of functionalities needed to perform a simulation. Applications should derive from this object to implement things specific to the application.

In coupled simulations everything that is not related to the physics of a problem is done in the `AnalysisStage`. This can be e.g. special user scripting or combined output.

These derived classes replace what was formerly done in `MainKratos.py`. This means that also the **user scripting should be in these classes: The idea is that instead of having a custom `MainKratos.py` the user derives a class from the `AnalysisStage` of the application to be used**. Only the functions require modifications are being overridden, the remaining implementation is used from the baseclass. This way updates to the baseclass are automatically being used in the users custom `AnalysisStage`.

## AnalysisStage: Responsibilities and provided Functionalities
The `AnalysisStage` handles everything not related to the physics of the problem. This includes e.g.
- Managing and calling the `PythonSolver`
- Construction and handling of the Processes
- Managing the output (post-processing in GiD/h5, saving restart, ...)

The main **public** functions are listed together with a brief explanation in the following. For a more detailed explanation it is referred to the docstrings of the respective functions.

- **Run**: this function executes the entire simulation
- **Initialize**: this function initializes the `AnalyisStage`, i.e. it performs all the operations necessary before the solution loop
- **RunSolutionLoop**: this function runs the solution loop
- **Finalize**: this function finalizes the `AnalyisStage`, i.e. it performs all the operations necessary after the solution loop

The main **protected** functions that are supposed to be used in derived classes are:
- **_GetSolver** : This function returns the `PythonSolver`. It also internally creates it if it does not exist yet
- **_GetListOfProcesses** : This function returns the list of processes. It will throw an error if the processes have not yet been created, since the creation happens after the ModelParts were read.
- **_GetListOfOutputProcesses** : This function returns the list of output processes. It will throw an error if the processes have not yet been created, since the creation happens after the ModelParts were read.

## AnalysisStage: Usage
In order to use the `AnalysisStage` it has to be constructed with specific objects:
- `KratosMultiphysics.Model`: The model containing all the modelparts involved in a simulation
- `KratosMultiphysics.Parameters`: The settings for the simulation. They are expecting that the following settings are present:
    * `problem_data` : general settings for the simulation
    * `solver_settings` : settings for the `PythonSolver`
    * `processes` : regular processes, e.g. for the boundary conditions. _Note_: also [user-defined processes can be included here](User-defined-python-processes)
    * `output_processes` : processes that write the output

```
{
"problem_data" : {
    "echo_level"    : 0
    "parallel_type" : "OpenMP" # or "MPI"
    "start_time"    : 0.0,
    "end_time"      : 1.0
},
"solver_settings" : {
...
settings for the PythonSolver
...
},
"processes" : {
    "my_processes" : [
    list of Kratos Processes
    ],
    "list_initial_processes" : [
    list of Kratos Processes
    ],
    "list_boundary_processes" : [
    list of Kratos Processes
    ],
    "list_custom_processes" : [
    list of Kratos Processes
    ]
},
"output_processes" : {
    "all_output_processes" : [
    list of Kratos Output Processes
    ]
}
```
***
Objects deriving from the `AnalysisStage` have to implement the `_CreateSolver` function which creates and returns the specific `PythonSolver`

**Note**: If the order in which the processes-blocks are initialized matters (if e.g. some processes would overwrite settings of other processes), then the function `_GetOrderOfProcessesInitialization` (resp. `_GetOrderOfOutputProcessesInitialization`) has to be overridden in the derived class. This function returns a list with the order in which the processes will be initialized.

As example we consider the settings above: We want the processes "list_initial_processes" to be constructed first and "list_custom_processes" to be constructed second. The order in which the other processes are initialized does not matter.
In this case we have to override the `_GetOrderOfProcessesInitialization` function to return `["list_initial_processes", "list_custom_processes"]`. With this we achieve the desired behavior.

# PythonSolver
## PythonSolver: Overview
The baseclass of the `PythonSolver` is located in the KratosCore ([here](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/python_solver.py)). It provides a set of functionalities that are needed for solving a physical problem. Applications should derive from this object to implement the application-specific tasks.

If physics are being coupled (e.g. for Fluid-Structure Interaction) then this should be implemented on solver-level. The coupled solver used the solvers of the involved physics and does also other tasks such as coupling logic or data exchange. E.g. an FSISolver would have a fluid and a structural solver.

## PythonSolver: Responsibilities and provided Functionalities
The `PythonSolver` is responsible for everything related to the physics of a problem. This includes e.g. 
- Settings up and solving of the system of equations
- Importing and preparing the ModelPart
- Advancing in time

The main **public** functions are listed together with a brief explanation in the following. For a more detailed explanation it is referred to the docstrings of the respective functions.
- **AddVariables** : this function adds the variables needed in the solution to the ModelPart
- **AddDofs** : this function adds the dofs needed in the solution to the ModelPart
- **ImportModelPart** : this function imports the ModelPart used by the solver (e.g. form an mdpa- or a restart-file)
- **PrepareModelPart** : this function prepares the ModelPart to be used by the solver (e.g. create SubModelParts necessary for the solution)
- **AdvanceInTime**: this function advances the `PythonSolver` in time
- **Initialize**: this function initializes the `PythonSolver`
- **Predict**: this function predicts the new solution
- **InitializeSolutionStep**: this function prepares solving a solutionstep
- **SolveSolutionStep**: this function solves a solutionstep
- **FinalizeSolutionStep**: this function finalizes solving a solutionstep
- **Finalize**: this function finalizes the `PythonSolver`

## PythonSolver: Usage
In order to use the `PythonSolver` it has to be constructed with specific objects:
- `KratosMultiphysics.Model`: The model to be used by the `PythonSolver`
- `KratosMultiphysics.Parameters`: The settings for the `PythonSolver`. They are expecting that the following settings are present:
    * `echo_level` : echo_level for printing informations
    * `model_import_settings` : settings for importing the modelpart
```
{
"echo_level" : 0,
"model_import_settings" : {
    "input_type"     : "mdpa" # or "rest"
    "input_filename" : "input_file_name"
}
```
For importing the ModelPart it can also be necessary (depending on the details of the solver) to pass the name of the ModelPart such that it can interact correctly with the Model

# Outlook (Kratos-Project, Multi-Stage Simulation)
**Note** This is a collection of ideas. Please note that the following is in a very early design phase.

In the future the objects presented here can be used in a larger context, e.g. a Multi-Stage Analysis.
This means that e.g. a FormFinding Analysis can be performed with doing a FSI-simulation afterwards and finally a fatigue-analysis

The `Model` as container of `ModelPart`s plays an important role here, since it is being used in every stage.
This is how data can be passed from one Stage to the next, i.e. if the `ModelPart` is the same then it will also contain e.g. the nodal-results from previous stages which can be reused then.

This could look like this:
```python
import KratosMultiphysics
# if needed import the apps containing the stage

model = KratosMultiphysics.Model() # only ONE for all stages!
# reading ProjectParameters

# construct all the stages
stage_1 = FormFindingStage(model, project_params_formfinding)
stage_2 = FSIStage(model, project_params_fsi)
stage_3 = FatigueAnalysisStage(model, project_params_fatigue)

list_of_stages = [stage_1, stage_2, stage_3]

for stage in list_of_stages:
    stage.Run()
```