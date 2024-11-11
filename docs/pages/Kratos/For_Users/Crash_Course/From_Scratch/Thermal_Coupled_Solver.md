---
title: Thermal coupled solver
keywords: 
tags: [Example Simple Fluid Thermal Coupled Solver]
sidebar: kratos_for_users
summary: 
---

Let's consider we wish to solve a simple fluid-thermal problem in which the velocity field obtained by a **CFD** analysis is to be employed as a basis in the related *thermal convection-diffusion* problem. You can download the example [here](https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Example:-A-simple-Fluid-Thermal-solver/coupled_thermal_cylinder.gid.zip).

Conceptually, the problem is a very simple one-way coupled problem, in which the temperature *convection-diffusion* problem is governed by the fluid flow. Provided that the same mesh is employed in the same domain, we would expect the solution to be equivalent to doing:

```python
fluid_solver.Solve() # Here we compute VELOCITY
thermal_solver.Solve() # Here we use VELOCITY to convect the temperature.
```

Current tutorial is about implementing such behaviour within the Kratos framework.

# Implementation

In the previous part of the tutorial we discussed the existance of two objecs: a __Stage__ and a __Solver__ .
The implementation of the proposed coupled problem involves dealing with the physics of the problem. The time stepping involved as well as the application of the Boundary Conditions is equivalent to that of a "normal" **CFD** problem, we will thus use the standard **Stage** employed in the **CFD** solution and implement a new `CoupledThermalSolver`

An example of implementation of such a solver can be found [here](https://github.com/KratosMultiphysics/Kratos/blob/Release-6.0/applications/FluidDynamicsApplication/python_scripts/coupled_fluid_thermal_solver.py)

however let's go step by step and consider the design of such a solver. 

The first step is to load in the memory the libraries needed. This is achieved by doing

```python
# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
```

Note that the libraries should be imported first in the main script, we verify this by calling the `CheckRegisteredApplications` to verify that the library is already loaded before entering in the current script.

## Solver Construction
The next step is to define a class implementing the behaviour we need. We do so by defining the class we want, as well as a `CreateSolver ` auxiliary function to be used in the object construction

```python
def CreateSolver(main_model_part, custom_settings):
    return CoupledFluidThermalSolver(main_model_part, custom_settings)

class CoupledFluidThermalSolver(object):
```

We now need to decide which solving technology shall be employed in solving the **CFD** and **Thermal** problems. 
Many possibilities exist to this end, since in principle any combination of transient solvers is viable. 

The best way to go is thus to postpone the decision on the choice of the solvers to runtime, so that the user is allowed to pick the choice out of the existing possibilities.
The concept here is simple: from the point of view of the coupling the two solvers are a black-box. We shall configure them as we would if we were to solve an uncoupled problem, by simply nesting the corresponding configuration file in the config being passed to the coupled solver. This is done in the class constructor by doing

```python
def __init__(self, model, custom_settings):
    
    self.model = model 
    
    default_settings = KratosMultiphysics.Parameters("""
    {
        "model_part_name" : "MainModelPart",
        "solver_type" : "ThermallyCoupled",
        "fluid_solver_settings": {
            "solver_type": "navier_stokes_solver_vmsmonolithic",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            }
        },
        "thermal_solver_settings": {
            "solver_type": "Transient",
            "analysis_type": "linear",
            "model_import_settings": {
                "input_type": "use_input_model_part",
                "input_filename": "unknown_name"
            },
            "computing_model_part_name": "Thermal",
            "material_import_settings": {
                "materials_filename": "ThermicMaterials.json"
            },
            "convection_diffusion_variables": {
                "density_variable": "DENSITY",
                "diffusion_variable": "CONDUCTIVITY",
                "unknown_variable": "TEMPERATURE",
                "volume_source_variable": "HEAT_FLUX",
                "surface_source_variable": "FACE_HEAT_FLUX",
                "projection_variable": "PROJECTED_SCALAR1",
                "convection_variable": "CONVECTION_VELOCITY",
                "mesh_velocity_variable": "MESH_VELOCITY",
                "transfer_coefficient_variable": "",
                "velocity_variable": "VELOCITY",
                "specific_heat_variable": "SPECIFIC_HEAT",
                "reaction_variable": "REACTION_FLUX"
            }
        }
    }""")

    ## Overwrite the default settings with user-provided parameters
    self.settings = custom_settings
    self.settings.ValidateAndAssignDefaults(default_settings)

    model_part_name = self.settings["model_part_name"].GetString()
    self.main_model_part = self.model[model_part_name]
    import python_solvers_wrapper_fluid
    self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(
        self.main_model_part, 
        self.settings["fluid_solver_settings"],"OpenMP")
    

    self.thermal_model_part = self.model.CreateModelPart("thermal_model_part")
    modeler = KratosMultiphysics.ConnectivityPreserveModeler()
    modeler.GenerateModelPart(
        self.main_model_part, 
        self.thermal_model_part, 
        "Element2D3N", 
        "Condition2D2N")
    
    import python_solvers_wrapper_convection_diffusion
    self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(
        self.thermal_model_part,
        self.settings["thermal_solver_settings"],
        "OpenMP")
```

If we take a detailed look, we observe that we added a variable `default_settings`, which tells which will be the default configuration involved in the solution phase. In this specific case, the default settings has a structure of the type

```python
default_settings = {
    "solver_type" : ...here the name of the current solver ...,
    "fluid_solver_settings": { ... CFD solver settings ... },
    "thermal_solver_settings": {... thermal solver settings ... }
}
```

Such default settings are overridden by the user-provided options, passed in the constructor as `settings`.
The user-provided settings may however contain wrong or unexpected input, or simply miss some information (for example not making a choice for the thermal solver).
We thus want to "validate" them against the defaults, which is achieved by the call:

```python
self.settings.ValidateAndAssignDefaults(default_settings)
```

In this specific case this function verifies that the inputs are `solver_type`,`fluid_solver_settings`,`thermal_solver_settings` and that they are respectively a string and two "nested parameters objects". The user-provided settings cannot provide any other entry. This allows to throw an error if a misspelled entry is provided as input.
In the case the input does not contain one of the given entries (for example the thermal_solver_settings) such entry is taken from the defaults.
An important feature of the  `ValidateAndAssignDefaults` is that it acts only on the outmost level. It does **NOT** do any validation on the subparameters.

The idea is that we can now construct the `fluid_solver` and the `thermal_solver` by passing the parameters. 
This is done by the lines

```python
import python_solvers_wrapper_fluid
self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(
    self.main_model_part, 
    self.settings["fluid_solver_settings"],
    "OpenMP")

import python_solvers_wrapper_convection_diffusion
self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(
    self.thermal_model_part,
    self.settings["thermal_solver_settings"],
    "OpenMP")
```

:warning: The syntax of these functions will probably change in the next future for something equivalent.

note that here we pass to the factory functions the subparameters as taken from the self.settings .
Validation of the relevant settings is thus demanded to each of the solvers, which will use hierarchically a mechanism similar to the one described here.

Aside of the details in any case the new solver now defined two inner solvers called 

```python
self.fluid_solver
self.thermal_solver
```

which are now ready for use. The rest of the class implementation now consists in simply calling the relevant functions for the new solver.

## AddVariables, AddDofs, ImportModelPart
The new solver implements the solver API, and thus defines a number of functions. Let's take a look at them:

The AddVariables function allocated memory for the variables needed. Essentially any solver tells to the corresponding modelpart which variables need to be stored in the database. For example, the **CFD** solver will require `VELOCITY `while the thermal solver will require `TEMPERATURE`.
:warning:  A hack is needed as of now, however this will be cleaned up in the next future

```python
def AddVariables(self):
    self.fluid_solver.AddVariables()
    
    #this is a HACK: TODO: cleanup so that this is not needed in the future 
    #the problem is that variables are not added to the fluid_solver.MainModelPart 
    #who is in charge of creating the nodes
    self.tmp = self.thermal_solver.main_model_part
    self.thermal_solver.main_model_part = self.fluid_solver.main_model_part
    
    self.thermal_solver.AddVariables()

    #this is a HACK: TODO: cleanup so that this is not needed in the future 
    self.thermal_solver.main_model_part = self.tmp
 ```

The `ImportModelPart `function is in charge of importing the modelparts involved.  
by taking a look on the default settings, we can observe that the fluid settings contain

```json
"fluid_solver_settings": {
    "model_import_settings": {
        "input_type": "mdpa",
        "input_filename": "unknown_name"
    }
}
```

while the thermal settings are slightly different:

```json
"thermal_solver_settings": {
    "model_import_settings": {
        "input_type": "use_input_model_part",
    }
}
```

this expresses that the `fluid_model_part` will be read from a file (named `input_filename`), while the `thermal_solver` expects the termal_solver to be already available. This is so as instead of reading the `thermal_model_part` from a file, we construct it on the basis of the `fluid_solver`, **so that it shares the same nodes**.
This is achieved by the class 
    
```python
modeler = KratosMultiphysics.ConnectivityPreserveModeler()
```

which fills the thermal modelpart with new elements and conditions, preserving however the same nodes (Pointers) as for the source modelpart (the fluid) as well as the same connectivity.

The overall importing function thus looks 

```python
def ImportModelPart(self):
    self.fluid_solver.ImportModelPart()
    
    #here cloning the fluid modelpart to thermal_model_part so that the nodes are shared
    modeler = KratosMultiphysics.ConnectivityPreserveModeler()
    modeler.GenerateModelPart(
        self.main_model_part,
        self.thermal_model_part,
        "Element2D3N",
            "Condition2D2N")

    self.thermal_solver.ImportModelPart()
```

The `AddDofs` function simply adds the Dofs needed in the solution problem and presents no difficulty
        
```python
def AddDofs(self):
    self.fluid_solver.AddDofs()
    self.thermal_solver.AddDofs()
```

## Solution Process
Even though conceptually the solution process consists in calling consecutively the two solvers, the solve function is split in several steps within kratos.

The Initialize function in particular is designed to be called exactly once, same as the `Check()` function.

```python
def Initialize(self):
    self.fluid_solver.Initialize()
    self.thermal_solver.Initialize()
```

Each complete solution step is organized as

```python
self.InitializeSolutionStep()
self.Predict()
self.SolveSolutionStep()
self.FinalizeSolutionStep()
```

where each of the substeps take into account the different solvers as in:

```python
def InitializeSolutionStep(self):
    self.fluid_solver.InitializeSolutionStep()
    self.thermal_solver.InitializeSolutionStep()

def Predict(self):
    self.fluid_solver.Predict()
    self.thermal_solver.Predict()

def SolveSolutionStep(self):
    self.fluid_solver.SolveSolutionStep()
    self.thermal_solver.SolveSolutionStep()
    

def FinalizeSolutionStep(self):
    self.fluid_solver.FinalizeSolutionStep()
    self.thermal_solver.FinalizeSolutionStep()
```

Where the functions:

- `InitializeSolutionStep`
- `FinalizeSolutionStep`

are to be called **EXACTLY ONCE** per time step.

Indeed for a simple one-way coupled problem as the one at hand, this is exactly equivalent to calling the `Solve` function for both solvers. The `Solve` functions are however split into several substeps as for  complex solvers two way coupled solvers one would do something equivalent to

```python
def Solve(self):
    # Do once per time step
    self.InitializeSolutionStep()
    self.Predict()

    # Non linear loop
    converged = False
    while(not converged):
        self.SolveSolutionStep()
        ...maybe here change something to make the problem two way coupled ...
        converged = ...here a convergence_check

    # Do once per time step
    self.FinalizeSolutionStep()
```

this subdivision is essential since by guaranteeing that `FinalizeSolutionStep` is called just once, we provide a safe place for updating internal variables.

## Advancing in time
The time stepping is dealt into the **Stage**. The time advancement from one step to the next is however done in the **Solver**. The releavant functions are:

The funciton `GetMinimumBufferSize` telling how many steps in the past need to be preserved

```python
    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)
```

and `AdvanceInTime` which actually makes the time to progress.
it is important to make a remark here: **since the nodes are shared between the modelparts the `CloneTimeStep` function should be done only on one of the modelparts**

```python
def AdvanceInTime(self, current_time):
    dt = self.ComputeDeltaTime()
    new_time = current_time + dt

    #NOTE: the cloning is done ONLY ONCE since the nodes are shared
    self.main_model_part.CloneTimeStep(new_time)
    self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

    return new_time
```
