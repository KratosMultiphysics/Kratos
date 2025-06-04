---
title: Pure diffusion solver derived from main python solver
keywords: 
tags: [Tutorial-Pure-diffusion-solver-derived-from-main-python-solver.md]
sidebar: kratos_for_developers
summary: 
---

TODO: Finish me

# Overview

# Derived solver

```python


# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MyLaplacianApplication")

# Import applications
import KratosMultiphysics.MyLaplacianApplication as Poisson

# Importing the base class
from python_solver import PythonSolver

# Other imports
import os


def CreateSolver(model, custom_settings):
    return PureDiffusionSolver(model, custom_settings)

class PureDiffusionSolver(PythonSolver):
    """The base class for pure diffusion solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_solution_scheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _create_solution_scheme
    _create_convergence_criterion
    _create_linear_solver
    _create_builder_and_solver
    _create_pure_diffusion_solution_strategy

    The pure_diffusion_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions get_pure_diffusion_solution_strategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, model, custom_settings):
        super(PureDiffusionSolver, self).__init__(model, custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "PureDiffusion",
            "domain_size" : -1,
            "echo_level": 0,
            "analysis_type": "linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "computing_model_part_name" : "PureDiffusionModelPart",
            "material_import_settings" :{
                "materials_filename": ""
            },
            "solver_type": "Stationary",
            "time_stepping" : { },
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "block_builder": true,
            "clear_storage": false,
            "move_mesh_flag": false,
            "convergence_criterion": "residual_criterion",
            "solution_relative_tolerance": 1.0e-4,
            "solution_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "AMGCL",
                "smoother_type":"ilu0",
                "krylov_type":"gmres",
                "coarsening_type":"aggregation",
                "max_iteration": 5000,
                "tolerance": 1e-9,
                "scaling": false
            },
            "element_replace_settings" : {
                "element_name" : "MyLaplacianElement",
                "condition_name" : "PointSourceCondition"
            },
            "problem_domain_sub_model_part_list": ["conv_diff_body"],
            "processes_sub_model_part_list": [""],
            "auxiliary_variables_list" : []
        }
        """)

        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings.AddEmptyValue("buffer_size")
        self.settings["buffer_size"].SetInt(self.GetMinimumBufferSize())
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
            self.solver_imports_model_part = False
        else:
            self.main_model_part = KratosMultiphysics.ModelPart(model_part_name) # Model.CreateodelPart()
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.solver_imports_model_part = True

        self.print_on_rank_zero("::[PureDiffusionSolver]:: ", "Construction finished")

    def AddVariables(self):
        ''' Add nodal solution step variables 
        '''
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE);
        self.main_model_part.AddNodalSolutionStepVariable(Poisson.POINT_HEAT_SOURCE);
            
        self.print_on_rank_zero("::[PureDiffusionSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        self.print_warning_on_rank_zero("::[PureDiffusionSolver]:: ", "Please define GetMinimumBufferSize() in your solver")
        return 1

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.REACTION_FLUX,self.main_model_part)
        self.print_on_rank_zero("::[PureDiffusionSolver]:: ", "DOF's ADDED")

    def ImportModelPart(self):
        """This function imports the ModelPart
        """
        if self.solver_imports_model_part:
            self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.is_restarted():
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()

            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors).Execute()
            
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()
        
            self._set_and_fill_buffer()

        # This will be removed once the Model is fully supported! => It wont e necessary anymore
        if not self.model.HasModelPart(self.main_model_part.Name):
            self.model.AddModelPart(self.main_model_part)
            
        if (self.settings["echo_level"].GetInt() > 0):
            self.print_on_rank_zero(self.model)

        KratosMultiphysics.Logger.PrintInfo("::[PureDiffusionSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        self.print_on_rank_zero("::[PureDiffusionSolver]:: ", "Initializing ...")
        # The convection_diffusion solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        pure_diffusion_solution_strategy = self.get_pure_diffusion_solution_strategy()
        pure_diffusion_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        if not self.is_restarted():
            pure_diffusion_solution_strategy.Initialize()
        else:
            # SetInitializePerformedFlag is not a member of SolvingStrategy but
            # is used by ResidualBasedNewtonRaphsonStrategy.
            try:
                pure_diffusion_solution_strategy.SetInitializePerformedFlag(True)
            except AttributeError:
                pass
        self.Check()
        self.print_on_rank_zero("::[PureDiffusionSolver]:: ", "Finished initialization.")

    def GetOutputVariables(self):
        pass

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        pure_diffusion_solution_strategy = self.get_pure_diffusion_solution_strategy()
        pure_diffusion_solution_strategy.Solve()

    def InitializeSolutionStep(self):
        self.get_pure_diffusion_solution_strategy().InitializeSolutionStep()

    def Predict(self):
        self.get_pure_diffusion_solution_strategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_pure_diffusion_solution_strategy().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.get_pure_diffusion_solution_strategy().FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def SetEchoLevel(self, level):
        self.get_pure_diffusion_solution_strategy().SetEchoLevel(level)

    def Clear(self):
        self.get_pure_diffusion_solution_strategy().Clear()

    def Check(self):
        self.get_pure_diffusion_solution_strategy().Check()

    #### Specific internal functions ####

    def get_solution_scheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._create_solution_scheme()
        return self._solution_scheme

    def get_convergence_criterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._create_convergence_criterion()
        return self._convergence_criterion

    def get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def get_builder_and_solver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver

    def get_pure_diffusion_solution_strategy(self):
        if not hasattr(self, '_pure_diffusion_solution_strategy'):
            self._pure_diffusion_solution_strategy = self._create_pure_diffusion_solution_strategy()
        return self._pure_diffusion_solution_strategy

    def import_materials(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]

    #### Private functions ####
    
    def _execute_after_reading(self):
        """Prepare computing model part and import constitutive laws. """
        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        # Assign mesh entities from domain and process sub model parts to the computing model part.
        import check_and_prepare_model_process_convection_diffusion as check_and_prepare_model_process
        check_and_prepare_model_process.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        if not self.model.HasModelPart(self.main_model_part.Name):
            self.model.AddModelPart(self.main_model_part)
        
        # Import constitutive laws.
        materials_imported = self.import_materials()
        if materials_imported:
            self.print_on_rank_zero("::[PureDiffusionSolver]:: ", "Materials were successfully imported.")
        else:
            self.print_on_rank_zero("::[PureDiffusionSolver]:: ", "Materials were not imported.")

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information. """
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        required_buffer_size = self.settings["buffer_size"].GetInt()
        if required_buffer_size < self.GetMinimumBufferSize():
            required_buffer_size = self.GetMinimumBufferSize()
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

    def _get_element_condition_replace_settings(self):
        # Duplicate model part
        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            for elem in self.main_model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break

        ## Elements
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            if (self.settings["element_replace_settings"]["element_name"].GetString() == "MyLaplacianElement"):
                if (num_nodes_elements == 3):
                    self.settings["element_replace_settings"]["element_name"].SetString("MyLaplacianElement")
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            raise Exception("Element not registered")
        else:
            raise Exception("DOMAIN_SIZE not set")
        
        ## Conditions
        num_nodes_conditions = 0
        if (len(self.main_model_part.Conditions) > 0):
            for cond in self.main_model_part.Conditions:
                num_nodes_conditions = len(cond.GetNodes())
                break
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            if (self.settings["element_replace_settings"]["condition_name"].GetString() == "PointSourceCondition"):
                if (num_nodes_conditions == 1):
                    self.settings["element_replace_settings"]["condition_name"].SetString("PointSourceCondition")
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            raise Exception("Condition not registered")
        else:
            raise Exception("DOMAIN_SIZE not set")

        return self.settings["element_replace_settings"]
    
    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("solution_relative_tolerance",self.settings["solution_relative_tolerance"])
        conv_params.AddValue("solution_absolute_tolerance",self.settings["solution_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        return conv_params

    def _create_convergence_criterion(self):
        import base_convergence_criteria_factory as convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.ConvergenceCriteriaFactory(self._get_convergence_criterion_settings())
        return convergence_criterion.convergence_criterion

    def _create_linear_solver(self):
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["block_builder"].GetBool():
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _create_solution_scheme(self):
        """Create the solution scheme for the structural problem.
        """
        solver_type = self.settings["solver_type"].GetString()
        if solver_type == "Stationary":
            convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        else:
            convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme() # Replace BDF and Euler
        return convection_diffusion_scheme
        

    def _create_pure_diffusion_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            pure_diffusion_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                pure_diffusion_solution_strategy = self._create_newton_raphson_strategy()
            else:
                pure_diffusion_solution_strategy = self._create_line_search_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return pure_diffusion_solution_strategy

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        pure_diffusion_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              pure_diffusion_scheme,
                                                              linear_solver,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        pure_diffusion_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        pure_diffusion_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                        pure_diffusion_scheme,
                                        linear_solver,
                                        pure_diffusion_convergence_criterion,
                                        builder_and_solver,
                                        self.settings["max_iteration"].GetInt(),
                                        self.settings["compute_reactions"].GetBool(),
                                        self.settings["reform_dofs_at_each_step"].GetBool(),
                                        self.settings["move_mesh_flag"].GetBool())

```

# Check and prepare model part

```python


# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    """Prepare the computing model part.

    The computing model part is created if it does not exist. Nodes and elements
    from the domain sub model parts are added to the computing model part.
    Conditions are added from the processes sub model parts.
    """
    def __init__(self, main_model_part, Parameters):
        self.main_model_part = main_model_part
        
        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        self.problem_domain_sub_model_part_list = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

    def Execute(self):
        
        problem_domain_sub_model_parts = []
        for i in range(self.problem_domain_sub_model_part_list.size()):
            problem_domain_sub_model_parts.append(self.main_model_part.GetSubModelPart(self.problem_domain_sub_model_part_list[i].GetString()))
        
        processes_parts = []
        for i in range(self.processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        #construct a model part which contains both the skin and the volume
        if (self.main_model_part.HasSubModelPart(self.computing_model_part_name)):
            computational_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        else:
            computational_model_part = self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        computational_model_part.Properties  = self.main_model_part.Properties
        
        for part in problem_domain_sub_model_parts:
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(computational_model_part, part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
            transfer_process.Execute()
        for part in processes_parts:
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(computational_model_part, part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS)
            transfer_process.Execute()

        KratosMultiphysics.Logger.PrintInfo("Computing model part:", computational_model_part)
```
