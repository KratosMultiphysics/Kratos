---
title: 5 - Analysis stage and Simulation Loop
keywords:
tags: [Kratos Crash Course analysis Stage Simulation Loop]
sidebar: kratos_for_users
summary:
---

# 1. Introduction

Up to this point, we've explored the principal components of Kratos: the essential files, data structures, and main components. However, we have yet to see how these elements come together to create a fully functional simulation workflow.

At the beginning of this course, you implemented a custom `AnalysisStage` class. This class is the core of Kratos' execution workflow, managing the main simulation loop, interactions with solvers, and the orchestration of processes and modelers. In this chapter, you will delve deeper into the components of the `AnalysisStage`, `Solvers` and `Processes`, examining how it integrates with Kratos solvers and auxiliary processes to execute a simulation.

# 2. The Analysis Stage

The `AnalysisStage` class is the principal component in a Kratos simulation. This class encapsulates the essential functions and method calls necessary for executing a simulation, providing a structured, adaptable interface to the Kratos framework. As previously mentioned, Kratos provides a generic `AnalysisStage`, which you can use directly or customize to meet more complex requirements.

It's important to note that the AnalysisStage itself does not contain any information about the specific physics of the problem. This responsibility lies within the `Elements` and `Solver` components selected for the simulation, allowing the `AnalysisStage` to remain flexible and applicable across different physical domains.

The AnalysisStage divides its execution into three distinct phases, which are sequentially invoked when the `Run` method is called. Each phase has a specific role in the simulation lifecycle.

You can find the complete sequence diagram for this class in our wiki, but here we will focus on the most relevant parts

Let's take a closer look at these three diferent stages using the generic analysis stage as a guide:

## 2.1 Initialize

This phase prepares the simulation environment. It loads model parts, initializes data structures, and prepares solver settings based on parameters from the configuration file. This setup phase is critical, as it ensures that the simulation begins in a well-defined state with all necessary components initialized:

*Note: we have ommited some code that is not interesting until you touch advanced topics*

```python
def Initialize(self):
    """This function initializes the AnalysisStage
    Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
    This function has to be implemented in deriving classes!
    """

    # Modelers:
    self._CreateModelers()
    self._ModelersSetupGeometryModel()
    self._ModelersPrepareGeometryModel()
    self._ModelersSetupModelPart()

    # Solver
    self._GetSolver().ImportModelPart()
    self._GetSolver().PrepareModelPart()
    self._GetSolver().AddDofs()

    ## Processes & Solver Initialization
    self.__CreateListOfProcesses()
    for process in self._GetListOfProcesses():
        process.ExecuteInitialize()

    self._GetSolver().Initialize()
    self.Check()

    for process in self._GetListOfProcesses():
        process.ExecuteBeforeSolutionLoop()

    # Get stepping and time settings
    self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

    # Get the time and the solver
    self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
    self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
```
{: data-lang="Python"}

This may look intimidating at first, but lets go step by step

- **Modelers**: The first block we encounter is `Modelers` block. We will se what a modeler is in section 4, and in this block we are creating the ones we will need.

```python
# Modelers:
self._CreateModelers()
self._ModelersSetupGeometryModel()
self._ModelersPrepareGeometryModel()
self._ModelersSetupModelPart()
```
{: data-lang="Python"}

- **Solver**: Next is the `Solver` block. As with `Modelers`, we will present a preview of solvers in section 3, but we can start to see some interesting things that should ring a bell on us: It Imports and prepares the modelpart.

    If you remember, we have stated several times that the solver was in charge of setting up the modelparts because it has de list of variables.

    Aside from setting up, it will also prepare the modelpart making the necessary changes to addapt it to a specific solver, and finally add the correct degrees of freedom.

```python
# Solver
self._GetSolver().ImportModelPart()
self._GetSolver().PrepareModelPart()
self._GetSolver().AddDofs()
```
{: data-lang="Python"}

- **Processess and Solver Initialization**: Is the section in charge of creating the user provided processes, which will se in section 4 and execute some of their parts. In this particular block we execute the `ExecuteInitialize` and `ExecuteBeforeSolutionLoop` methods, and in the middle we initialize the solver we created in the previous section. We do it this way because some processes need to be executed before the solver is initialized, and some need to be executed after we have some of the information that the solved provides once initialized. This way, processes give enouch flexibility to approach both situations:

```python
## Processes
self.__CreateListOfProcesses()
for process in self._GetListOfProcesses():
    process.ExecuteInitialize()

self._GetSolver().Initialize()
self.Check()

for process in self._GetListOfProcesses():
    process.ExecuteBeforeSolutionLoop()
```
{: data-lang="Python"}

- **Time and Timestep**: Finaly, as the analysis stage is the class in control of the time loop, the last task of the initialization stage is to provide the initial, final and time step of our simulation:

```python
# Get stepping and time settings
self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

# Get the time and the solver
self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
```
{: data-lang="Python"}

## 2.2 Execution

The main simulation loop occurs here. During this phase, the AnalysisStage manages iterative solution steps, invoking solvers and applying boundary conditions or other processes as needed. This phase continues until the simulation reaches the predefined end conditions or maximum time steps:

```python
def RunSolutionLoop(self):
    """This function executes the solution loop of the AnalysisStage
    It can be overridden by derived classes
    """
    while self.KeepAdvancingSolutionLoop():
        self.time = self._AdvanceTime()
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        is_converged = self._GetSolver().SolveSolutionStep()
        self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()
```
{: data-lang="Python"}

As you can see at first glance the logic is much more simpler

- 1) **Advance the simulation time**
- 2) **Initialize the solution step**
- 3) **Predict a initial solution**
- 4) **Perform the solution of the system**
- 5) **Finalize the solution step**
- 6) **Print and output**

This loop continues to execute until the KeepAdvancingSolutionLoop method returns False. This method governs the loop’s termination, which may occur for several reasons:

- The final simulation time has been reached.
- The solution has converged.
- Other user-defined criteria have been met.

Two additional methods are particularly important within this loop: `InitializeSolutionStep` and `FinalizeSolutionStep`. Similar in concept to the `Initialize` and `Finalize` methods, these functions prepare and clean up data specific to each time step:

```python
def InitializeSolutionStep(self):
    """This function performs all the required operations that should be executed
    (for each step) BEFORE solving the solution step.
    """
    self.PrintAnalysisStageProgressInformation()

    for process in self._GetListOfProcesses():
        process.ExecuteFinalizeSolutionStep()

    self.ChangeMaterialProperties() # This is normally empty
    self._GetSolver().InitializeSolutionStep()
```
{: data-lang="Python"}

As you can see, the `InitializeSolutionStep` prints info to let the user knwon that we are advancing in time, calls the `ExecuteInitializeSolutionStep` from our processes (for historical reasons, this may appear as `ApplyBoundaryConditions` in some analysis stages), and then call same method from the solver.

```python
def FinalizeSolutionStep(self):
    """This function performs all the required operations that should be executed
    (for each step) AFTER solving the solution step.
    """
    self._GetSolver().FinalizeSolutionStep()

    for process in self._GetListOfProcesses():
        process.ExecuteFinalizeSolutionStep()
```
{: data-lang="Python"}

Anologously, the `FinalizeSolutionStep` perform the same steps (except for changing the material properties)

## 2.3 Finalize

Once the simulation loop concludes, the `Finalize` phase handles any post-processing or data cleanup tasks. In the generic stage this consist on calling the `ExecuteFinalize` method from the processes and solver but may also include saving results, releasing memory, or performing additional actions required for analysis after the simulation run.

```python
def Finalize(self):
    """This function finalizes the AnalysisStage
    Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
    """
    for process in self._GetListOfProcesses():
        process.ExecuteFinalize()

    self._GetSolver().Finalize()
```
{: data-lang="Python"}

# 3. The Solver

The solver is the responsible of providing the simulation physics and solving the problem.

Hence, it can be easily guessed that there are as many solvers as physics and ways to solve such physics.

Some well-established examples are the `structural_mechanics_solver.py`, `fluid_solver` and `convection_diffusion_solver.py` as well as their derived classes (e.g., `structural_mechanics_static_solver.py` or `structural_mechanics_dynamic_implicit_solver.py`).

At this point, we should remark the following points:

- Solvers are always implemented in Python level for the sake of flexibility.

- All solvers in Kratos have `PythonSolver` class, which can be found in `python_solver.py`, as root class. This solver has no physics and its main purpose is to define the solver class API.

- It is a common, and very good by the way, practice to have a base solver for one kind of physics and then deriving several solvers from it to implement the particularities. For instance, there is a unique fluid solver but then there are several derived ones according to the resolution scheme or particular physics (e.g., segregated or monolithic, single phase or two-phase, etc.).


It has been already stated that the main purpose of the solver is to define and solve the physics of the simulation.

In a nutshell, this boils down to the following

- Validating the `solver_settings` in the `ProjectParameters.json` with the corresponding defaults in the solver.

- Adding the required nodal historical variables to the database. Note that this is done according to the corresponding elements and conditions that actually implement the physics (i.e., the variational form).

- Adding the required nodal DOFs similar to what is done with the variables.


Complementary, it is also due mentioning that, historically, the solver also had the purpose of importing the geometry. However, this is a feature that is becoming deprecated in favour of the geometry-based input at the modeler level. Nevertheless, we consider this a technical detail that is out of the scope of this crash course.

Last but important, the solver also creates and holds the strategy (e.g., the Newton-Raphson instance) and linear solver (if required) that actually solve the problem. Again, these are considered more advance features that we prefer to leave out of the scope of the crash course.


The solver can be retrieved from the `AnalysisStage` level by calling the `_GetSolver()` function.

This becomes specially handy when combined with the `GetComputingModelPart()`, namely `_GetSolver().GetComputingModelPart()`.

Note that by doing so, we can access (and play around with or customize) the database of the model part that is actually used for the problem resolution, something that enables an easy customization of the simulation when called in the proper access points of the `AnalysisStage`.

# 4. Processe

Processes and Modelers are the way we provide our users to make customizations to the simulation without having to touch any of the core components described in this section. Essentially, they are pieces of code that will be executed in selected points during the analysis stage.

Processes in particular are the backbone for custom code execution in the analysis stage. As we have seen during the analysis stage they are use in different moments during the simulation execution.

The most important characteristic  of a process is being a class with a fix set of entry points:

- `ExecuteInitialize`: Will be called during the initialize sequence of an `AnalysisStage`, before the initialization of the `Solver`.

- `ExecuteBeforeSolustionLoop`: Will be called during the initialize sequence of an `AnalysisStage`, after the initialization of the `Solver`

- `ExecuteInitializeSolutionStep`: Will be called at the begining of each solution loop, before executing the preconditioners and solvers

- `ExecuteFinalizeSolutionStep`: Will be called at the end of each solution loop, after executing the preconditioners and solvers but before the output stage.

- `ExecuteBeforeOutputStep`: Will be called at the begining of the output stage for every output process active, before printing the results

- `ExecuteAfterOutputStep`: Will be called at the end of the output stage for every output process active, after printing the results

- `ExecuteFinalize`: Will be called during the finalize sequence of an `AnalysisStage`, just before existing the stage.

In order to use a process, you have to add it to one of the lists available in the project parameters with the name of the process that you want to use, and the parameters it expects. There are plenty of processes in any of the examples available through kratos, but just to illustrate let's take a look at one in particular from the ones you downloaded:

```json
"constraints_process_list" : [{
"python_module" : "assign_vector_variable_process",
"kratos_module" : "KratosMultiphysics",
"Parameters"    : {
    "model_part_name" : "Structure.DISPLACEMENT_Ground",
    "variable_name"   : "DISPLACEMENT",
    "constrained"     : [true,true,true],
    "value"           : [0.0,0.0,0.0],
    "interval"        : [0.0,"End"]
}
}],
```

As you can see, here we are seeing the `assign_vector_variable_process` which is found in the module `KratosMultiphysics`. If you were to use a process from another application, the `kratos_module` would change.

The list of `Parameters` are the settings that the process expects, and they may broadly change between processes.

In particular, this process assigns a value to a vector variable, in this case `DISPLACEMENT`. This is typically at operation that you will want to do before solving a step during a given amount of time, indicated in the `interval`.

If we look at the implementation of such process in Kratos, we will see that, among the different entry points that we have commented, this process makes use of `ExecuteInitializeSolutionStep`, with a filter for the time:


```python
def ExecuteInitializeSolutionStep(self):
    for process in self.aux_processes:
        process.ExecuteInitializeSolutionStep()
```
In particular, this process uses a list of subprocesses, which is something that you can also do. The subprocess follows the same logic:


```python
def ExecuteInitializeSolutionStep(self):
    current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

    if self.interval.IsInInterval(current_time):
            self.value = self.table.GetValue(current_time)
            self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)
```

Here we can see a simplified version of its implementation, and as you can see the behavior is the expected.

# 5. Modelers

# 6. Utilities

Finally utilities play a similar role as processes and modeleres, but there is no fixed entry points in which they will be called during a simulation.

The role of utilities is to provide a mechanism to ensure that functions that may be useful for other users are encapuslated and can be used, but is your responsability as a programmer to call them whenever necessary, for example inside a process, or in a function that you created in your custom analysis stage.

Let's see some examples: