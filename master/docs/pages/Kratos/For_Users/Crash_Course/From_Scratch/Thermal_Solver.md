---
title: Implementing Thermal Solver
keywords: 
tags: [Implementing Thermal Solver]
sidebar: kratos_for_users
summary: 
---

The following tutorial will explain how to **implement a `Solver` from scratch**, in this case applied to the particular case of the *thermal problem*, following the lines of the already presented on the other tutorials. We will skip the most advanced points of the construction of the `Solver` and we will create a **basic working solver** for the sake of academic purposes, so the resulting file will not coincide with the real solver on the [repository](https://github.com/KratosMultiphysics/Kratos/blob/Release-6.0/applications/convection_diffusion_application/python_scripts/convection_diffusion_base_solver.py). 

## Imports 

The solver as any element of *Kratos* requires to import the corresponding libraries and applications. Like out objective on mind is to create a thermal problem, we will import the base `KratosMultiphysics` library as well as the `ConvectionDiffusionApplication`. For importing files from the filesystem we will import `os` *python* library, which will be helpful.

```python
# Python imports
import os

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
``` 

# Constructing the solver

We will create a `Solver`, following was is done in other applications `Solver`, for the sake of consistency, so in first place we define the function `CreateSolver`, which is common among all the solvers and therefore it is necessary to be called that way. This function will use as input a `Model` and configuration `Parameters`:

```python
def CreateSolver(model, custom_settings):
    return ConvectionDiffusionSolver(model, custom_settings)
```

As we see this function call a class called `ConvectionDiffusionSolver`, so we need to define our solver in first place, for that we define the constructor or `__init__` function:

```python
class ConvectionDiffusionSolver(object):
    """The base class for convection-diffusion solvers.
    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.
    """
    def __init__(self, model, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "ThermalModelPart",
            "echo_level": 0,
            "buffer_size": 2,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "computing_model_part_name" : "Thermal",
            "material_import_settings" :{
                "materials_filename": ""
            },
            "clear_storage": false,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "BICGSTABSolver",
                "preconditioner_type": "DiagonalPreconditioner",
                "max_iteration": 5000,
                "tolerance": 1e-9,
                "scaling": false
            },
            "element_replace_settings" : {
                "element_name" : "EulerianConvDiff",
                "condition_name" : "Condition"
            },
            "problem_domain_sub_model_part_list": ["conv_diff_body"],
            "processes_sub_model_part_list": [""]
        }
        """)
    
        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        model_part_name = self.settings["model_part_name"].GetString()
        self.main_model_part = self.model[model_part_name]
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Construction finished")
```

First we define the default `Parameters`, which will overwrite the not define configurations of our `custom_settings`, this is done with the `ValidateAndAssignDefaults` method. We define the `main_model_part` as a part of the solver, using the `self` to define it as an instance attribute.

## Common methods for the solver

The following are common methods defined in all solvers in order to have an interoperability between them to couple different physical problems.

### Add variables

The `AddVariables` is the liable on adding the nodal historical variables to the model part. We will need to define a buffer in order to be able to access this historical variables. We will see how in following sections.

We add the variables corresponding to the standard thermal convection-diffusion problem. We will need to add this variables to the `CONVECTION_DIFFUSION_SETTINGS` (this is something specific of the `ConvectionDiffusionAplication`, and is done in order to solve different types convection diffusion problems).

```python
def AddVariables(self):
    ''' Add nodal solution step variables ( later add to CONVECTION_DIFFUSION_SETTINGS)
    '''
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)
    self.main_model_part.AddNodalSolutionStepVariable(ConvectionDiffusionApplication.PROJECTED_SCALAR1)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONVECTION_VELOCITY)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_FLUX)

    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Variables ADDED")
```

### Add DoFs

In this case we add to the nodes the degrees of freedom corresponding to the problems we want to solve, in our case only the `TEMPERATURE`.

```python
def AddDofs(self):
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.REACTION_FLUX, self.main_model_part)
    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "DOF's ADDED")
```

We the `VariableUtils` class, which allows us to do this operation in parallel, in case we want to use the methods available directly on the `Node` class we can do:

```python
def AddDofs(self):
    for node in self.main_model_part.Nodes:
        node.AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.REACTION_FLUX)          
    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "DOF's ADDED")  
```

### Read external file (`*.mdpa` file)

We read the `*.mdpa` file for that we use the class `ModelPartIO`, which is the class in charge of managing the write/read procedures.

```python
def ReadModelPart(self):
    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Reading model part.")
    problem_path = os.getcwd()
    input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
    # Import model part from mdpa file.
    KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
    KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Finished reading model part.")
```

### Prepare model part

After reading the model part we need to do some operations, as checking the mesh orientation, assigning the physical elements and conditions to the *dummy* `Element`/`Condition` and filling the buffer. 

We mean for *dummy* ` Condition` and `Element` that only contain a geometry and don't solve any physical problem, this is done for the sake of interoperability, which allows us to use the same mesh for the fluid and thermal problem.

First we define the reading of the properties of the problem. We use the  `read_materials_process` defined on the core of *Kratos*.

```python
def _execute_after_reading(self):
    """Import materials. """
    materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
    if (materials_filename != ""):
        import read_materials_process
        # Create a dictionary of model parts.
        # Add constitutive laws and material properties from json file to model parts.
        read_materials_process.ReadMaterialsProcess(self.model, self.settings["material_import_settings"])
        
        # We set the properties that are nodal
        self._assign_nodally_properties()
        
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Materials were successfully imported.")
    else:
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Materials were not imported.")
```

In the case of the convection diffusion problems, the properties are not read from the `Property` but directly from the nodal values, so we create a method in order to copy this properties into the nodes of the mesh:

```python
def _assign_nodally_properties(self):
    # We transfer the values of the con.diff variables to the nodes
    with open(self.settings["material_import_settings"]["materials_filename"].GetString(), 'r') as parameter_file:
        materials = KratosMultiphysics.Parameters(parameter_file.read())
        
    for i in range(materials["properties"].size()):
        model_part = self.main_model_part.GetSubModelPart(materials["properties"][i]["model_part_name"].GetString())
        mat = materials["properties"][i]["Material"]
        
        for key, value in mat["Variables"].items():
            var = KratosMultiphysics.KratosGlobals.GetVariable(key)
            if (self._check_variable_to_set(var)):
                if value.IsDouble():
                    KratosMultiphysics.VariableUtils().SetScalarVar(var, value.GetDouble(), model_part.Nodes)
                elif value.IsVector():
                    KratosMultiphysics.VariableUtils().SetVectorVar(var, value.GetVector(), model_part.Nodes)
                else:
                    raise ValueError("Type of value is not available")
```

We define the `Parameters` necessaries for the replace the `Element`/`Condition`, which depend on the configuration parameters, and the geometry of the problem (2D/3D, triangle, quadrilaterals, etc...). We will use the class `ReplaceElementsAndConditionsProcess` to replace these elements and conditions.

```python
def _get_element_condition_replace_settings(self):
    num_nodes_elements = 0
    if (len(self.main_model_part.Elements) > 0):
        num_nodes_elements = len(self.main_model_part.Elements[1].GetNodes())

    ## Elements
    if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
        if (num_nodes_elements == 3):
            self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D")
        else:
            self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D4N")
    elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
        if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
            self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff3D")
    else:
        raise Exception("DOMAIN_SIZE not set")
    
    ## Conditions
    num_nodes_conditions = 0
    if (len(self.main_model_part.Conditions) > 0):
        num_nodes_conditions = len(self.main_model_part.Conditions[1].GetNodes())
    if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
        self.settings["element_replace_settings"]["condition_name"].SetString("LineCondition2D2N")
    elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
        self.settings["element_replace_settings"]["condition_name"].SetString("SurfaceCondition3D3N")
    else:
        raise Exception("DOMAIN_SIZE not set")

    return self.settings["element_replace_settings"]
```

Fill the buffer is as easy as:

```python
def _set_and_fill_buffer(self):
    """Prepare nodal solution step data containers and time step information. """
    # Set the buffer size for the nodal solution steps data. Existing nodal
    # solution step data may be lost.
    required_buffer_size = self.settings["buffer_size"].GetInt()
    current_buffer_size = self.main_model_part.GetBufferSize()
    buffer_size = max(current_buffer_size, required_buffer_size)
    self.main_model_part.SetBufferSize(buffer_size)
    # Cycle the buffer. This sets all historical nodal solution step data to
    # the current value and initializes the time stepping in the process info.
    delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
    time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
    step =-buffer_size
    time = time - delta_time * buffer_size
    self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
    for i in range(0, buffer_size):
        step = step + 1
        time = time + delta_time
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
        self.main_model_part.CloneTimeStep(time)
```

Finally, now that we have all the tools to prepare the model part we can define the method:

```python
def PrepareModelPartForSolver(self):
        
    # Check and prepare computing model part and import constitutive laws.
    self._execute_after_reading()

    throw_errors = False
    KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors).Execute()
    KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()

    self._set_and_fill_buffer()
    
    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "ModelPart prepared for Solver.")
```

### Initialize method

In order to initialize the solver we need to create the strategy, which will be *Newton-Raphson* in our case
. Additionally we will define the `CONVECTION_DIFFUSION_SETTINGS` which can be used by the elements and conditions from the application.

```python
def Initialize(self):
    """Perform initialization after adding nodal variables and dofs to the main model part. """
    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Initializing ...")

    # We define the convection diffusion settings
    convention_diffusion_settings = KratosMultiphysics.ConvectionDiffusionSettings()
    convention_diffusion_settings.SetDensityVariable(KratosMultiphysics.DENSITY)
    convention_diffusion_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
    convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
    convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.HEAT_FLUX)
    convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.FACE_HEAT_FLUX)
    convention_diffusion_settings.SetProjectionVariable(ConvectionDiffusionApplication.PROJECTED_SCALAR1)
    convention_diffusion_settings.SetConvectionVariable(KratosMultiphysics.CONVECTION_VELOCITY)
    convention_diffusion_settings.SetMeshVelocityVariable(KratosMultiphysics.MESH_VELOCITY)
    convention_diffusion_settings.SetVelocityVariable(KratosMultiphysics.VELOCITY)
    convention_diffusion_settings.SetSpecificHeatVariable(KratosMultiphysics.SPECIFIC_HEAT)
    convention_diffusion_settings.SetReactionVariable(KratosMultiphysics.REACTION_FLUX)
        
    self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convention_diffusion_settings)
        
    # The mechanical solution strategy is created here if it does not already exist.
    computing_model_part = self.main_model_part # We will use the main model part
    #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
    computing_model_part.ProcessInfo[ConvectionDiffusionApplication.THETA] = 0.5
    computing_model_part.ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 1.0
    conv_diff_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
    import linear_solver_factory
    linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
    conv_diff_convergence_criterion = KratosMultiphysics.ResidualCriteria(self.settings["residual_relative_tolerance"].GetDouble(), self.settings["residual_absolute_tolerance"].GetDouble())
    builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
    self.conv_diff_strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,  conv_diff_scheme, linear_solver, conv_diff_convergence_criterion, builder_and_solver, self.settings["max_iteration"].GetInt(), True, False, False)
    self.conv_diff_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
    self.conv_diff_strategy.Initialize()
    self.conv_diff_strategy.Check()
    if self.settings["clear_storage"].GetBool():
        self.conv_diff_strategy.Clear()
    KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Finished initialization.")
```

### Solving methods

The `Solver` as well as the strategy, has the same solving steps as the strategies (*Newton-Rapshon* and so on), this methods are `InitializeSolutionStep`, `Predict`, `SolveSolutionStep` and `FinalizeSolutionStep`, which can be call all together in only one method `Solve`.

```python
def Solve(self):
    if self.settings["clear_storage"].GetBool():
        self.conv_diff_strategy.Clear()
    self.conv_diff_strategy.Solve()

def InitializeSolutionStep(self):
    self.conv_diff_strategy.InitializeSolutionStep()

def Predict(self):
    self.conv_diff_strategy.Predict()

def SolveSolutionStep(self):
    is_converged = self.conv_diff_strategy.SolveSolutionStep()
    return is_converged

def FinalizeSolutionStep(self):
    self.conv_diff_strategy.FinalizeSolutionStep()
```

## Final script

Our thermal solver will be like the following:

```python
# Python imports
import os

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionSolver(main_model_part, custom_settings)

class ConvectionDiffusionSolver(object):
    """The base class for convection-diffusion solvers.
    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.
    """
    def __init__(self, model, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "ThermalModelPart",
            "echo_level": 0,
            "buffer_size": 2,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "computing_model_part_name" : "Thermal",
            "material_import_settings" :{
                "materials_filename": ""
            },
            "clear_storage": false,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "BICGSTABSolver",
                "preconditioner_type": "DiagonalPreconditioner",
                "max_iteration": 5000,
                "tolerance": 1e-9,
                "scaling": false
            },
            "element_replace_settings" : {
                "element_name" : "EulerianConvDiff",
                "condition_name" : "Condition"
            },
            "problem_domain_sub_model_part_list": ["conv_diff_body"],
            "processes_sub_model_part_list": [""]
        }
        """)
    
        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        model_part_name = self.settings["model_part_name"].GetString()
        self.main_model_part = self.model[model_part_name]
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Construction finished")

    def AddVariables(self):
        ''' Add nodal solution step variables ( later add to CONVECTION_DIFFUSION_SETTINGS)
        '''
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(ConvectionDiffusionApplication.PROJECTED_SCALAR1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONVECTION_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_FLUX)

        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            node.AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.REACTION_FLUX)          
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "DOF's ADDED")  

    def ReadModelPart(self):
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Reading model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
        # Import model part from mdpa file.
        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Finished reading model part.")

    def _execute_after_reading(self):
        """Import materials. """
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            import read_materials_process
            # Add constitutive laws and material properties from json file to model parts.
            read_materials_process.ReadMaterialsProcess(self.model, self.settings["material_import_settings"])
            
            # We set the properties that are nodal
            self._assign_nodally_properties()
            
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Materials were successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Materials were not imported.")

    def _assign_nodally_properties(self):
        # We transfer the values of the con.diff variables to the nodes
        with open(self.settings["material_import_settings"]["materials_filename"].GetString(), 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())
            
        for i in range(materials["properties"].size()):
            model_part = self.main_model_part.GetSubModelPart(materials["properties"][i]["model_part_name"].GetString())
            mat = materials["properties"][i]["Material"]
            
            for key, value in mat["Variables"].items():
                var = KratosMultiphysics.KratosGlobals.GetVariable(key)
                if (self._check_variable_to_set(var)):
                    if value.IsDouble():
                        KratosMultiphysics.VariableUtils().SetScalarVar(var, value.GetDouble(), model_part.Nodes)
                    elif value.IsVector():
                        KratosMultiphysics.VariableUtils().SetVectorVar(var, value.GetVector(), model_part.Nodes)
                    else:
                        raise ValueError("Type of value is not available")

    def _get_element_condition_replace_settings(self):
        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            num_nodes_elements = len(self.main_model_part.Elements[1].GetNodes())

        ## Elements
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            if (num_nodes_elements == 3):
                self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D")
            else:
                self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D4N")
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
                self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff3D")
        else:
            raise Exception("DOMAIN_SIZE not set")
        
        ## Conditions
        num_nodes_conditions = 0
        if (len(self.main_model_part.Conditions) > 0):
            num_nodes_conditions = len(self.main_model_part.Conditions[1].GetNodes())
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            self.settings["element_replace_settings"]["condition_name"].SetString("LineCondition2D2N")
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            self.settings["element_replace_settings"]["condition_name"].SetString("SurfaceCondition3D3N")
        else:
            raise Exception("DOMAIN_SIZE not set")

        return self.settings["element_replace_settings"]

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information. """
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        required_buffer_size = self.settings["buffer_size"].GetInt()
        current_buffer_size = self.main_model_part.GetBufferSize()
        buffer_size = max(current_buffer_size, required_buffer_size)
        self.main_model_part.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            self.main_model_part.CloneTimeStep(time)

    def PrepareModelPartForSolver(self):
            
        # Check and prepare computing model part and import constitutive laws.
        self._execute_after_reading()

        throw_errors = False
        KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors).Execute()
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()

        self._set_and_fill_buffer()
        
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Initializing ...")

        # We define the convection diffusion settings
        convention_diffusion_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        convention_diffusion_settings.SetDensityVariable(KratosMultiphysics.DENSITY)
        convention_diffusion_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
        convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
        convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.HEAT_FLUX)
        convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        convention_diffusion_settings.SetProjectionVariable(ConvectionDiffusionApplication.PROJECTED_SCALAR1)
        convention_diffusion_settings.SetConvectionVariable(KratosMultiphysics.CONVECTION_VELOCITY)
        convention_diffusion_settings.SetMeshVelocityVariable(KratosMultiphysics.MESH_VELOCITY)
        convention_diffusion_settings.SetVelocityVariable(KratosMultiphysics.VELOCITY)
        convention_diffusion_settings.SetSpecificHeatVariable(KratosMultiphysics.SPECIFIC_HEAT)
        convention_diffusion_settings.SetReactionVariable(KratosMultiphysics.REACTION_FLUX)
            
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convention_diffusion_settings)
            
        # The mechanical solution strategy is created here if it does not already exist.
        computing_model_part = self.main_model_part # We will use the main model part
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        computing_model_part.ProcessInfo[ConvectionDiffusionApplication.THETA] = 0.5
        computing_model_part.ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 1.0
        conv_diff_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        conv_diff_convergence_criterion = KratosMultiphysics.ResidualCriteria(self.settings["residual_relative_tolerance"].GetDouble(), self.settings["residual_absolute_tolerance"].GetDouble())
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        self.conv_diff_strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,  conv_diff_scheme, linear_solver, conv_diff_convergence_criterion, builder_and_solver, self.settings["max_iteration"].GetInt(), True, False, False)
        self.conv_diff_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        self.conv_diff_strategy.Initialize()
        self.conv_diff_strategy.Check()
        if self.settings["clear_storage"].GetBool():
            self.conv_diff_strategy.Clear()
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Finished initialization.")

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.conv_diff_strategy.Clear()
        self.conv_diff_strategy.Solve()

    def InitializeSolutionStep(self):
        self.conv_diff_strategy.InitializeSolutionStep()

    def Predict(self):
        self.conv_diff_strategy.Predict()

    def SolveSolutionStep(self):
        is_converged = self.conv_diff_strategy.SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.conv_diff_strategy.FinalizeSolutionStep()
```