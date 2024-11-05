---
title: Analysis stage and Simulation Loop
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

# 2.1 Initialize

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

This may look intimidating at first, but lets go step by step

1) **Modelers**: The first block we encounter is `Modelers` block. We will se what a modeler is in section 4, and in this block we are creating the ones we will need.

```python
# Modelers:
self._CreateModelers()
self._ModelersSetupGeometryModel()
self._ModelersPrepareGeometryModel()
self._ModelersSetupModelPart()
```

- **Solver**: Next is the `Solvers` block. As with `Modelers`, we will present a preview of solvers in section 3, but we can start to see some interesting things that should ring a bell on us: It Imports and prepares the modelpart.

    If youu remember, we have stated several times that the solver was in charge of reading modelparts because it has de list of variables, and is here that the code reads our modelpart.

    Aside from reading, we will also prepare the modelpart making the necessary changes to addapt it to a specific solver, and finally add the correct degrees of freedom.

```python
# Solver
self._GetSolver().ImportModelPart()
self._GetSolver().PrepareModelPart()
self._GetSolver().AddDofs()
```

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

- **Time and Timestep**: Finaly, as the analysis stage is the class in control of the time loop, the last task of the initialization stage is to provide the initial, final and time step of our simulation:

```python
# Get stepping and time settings
self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

# Get the time and the solver
self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
```

# 2.2 Execution

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

Anologously, the `FinalizeSolutionStep` perform the same steps (except for changing the material properties)


# 2.3 Finalize

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

# 3. The Solver

The solver is the responsible of building and solving your system given the information of your geometries, variables and configuration. There are many types of solvers inside Kratos, and it goes beyond the scope of this course to explain all details. 

Explain interface and _GetComputingPart


# 4. Processes and Modeleres

Processes, Utilities and Modelers are the way we privde our users to make small customizations to the simulation without having to touch any of the core components descrived in this section. Essentially they are pieces of code that will be executed in selected points during the analysis stage. 

# 5. Utilities

Finally utilities play a similar role as processes and modeleres, but there is no fixed entry porints in which they will be called during a simulation inside the analysis stage.

The role of utilities is to provide a mechanism to ensure that functions that may be usefull for other users are encapuslated and can be used, but is your responsability as a programmer to call them whenever necessary, for example inside a process, or in a function that you created in your custom analysis stage.

Let's make a couple of examples: