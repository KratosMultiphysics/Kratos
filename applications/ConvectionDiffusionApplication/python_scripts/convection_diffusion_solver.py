from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
import KratosMultiphysics.base_convergence_criteria_factory as convergence_criteria_factory

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.ConvectionDiffusionApplication.check_and_prepare_model_process_convection_diffusion as check_and_prepare_model_process

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(model, custom_settings):
    return ConvectionDiffusionSolver(model, custom_settings)

class ConvectionDiffusionSolver(PythonSolver):
    """The base class for convection-diffusion solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_solution_scheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _create_solution_scheme
    _create_convergence_criterion
    _create_linear_solver
    _create_builder_and_solver
    _create_convection_diffusion_solution_strategy

    The convection_diffusion_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions get_convection_diffusion_solution_strategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass = True
        super().__init__(model, custom_settings)

        # Convection diffusion variables check
        self._ConvectionDiffusionVariablesCheck(custom_settings)

        model_part_name = self.settings["model_part_name"].GetString()

        # Set default buffer size
        self.min_buffer_size = 1

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
            self.solver_imports_model_part = False
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name) # Model.CreateodelPart()
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.solver_imports_model_part = True

        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "ThermalModelPart",
            "domain_size" : -1,
            "echo_level": 0,
            "analysis_type": "linear",
            "solver_type": "convection_diffusion_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "computing_model_part_name" : "thermal_computing_domain",
            "material_import_settings" :{
                "materials_filename": ""
            },
            "convection_diffusion_variables" : {
                "density_variable"              : "DENSITY",
                "diffusion_variable"            : "CONDUCTIVITY",
                "unknown_variable"              : "TEMPERATURE",
                "volume_source_variable"        : "HEAT_FLUX",
                "surface_source_variable"       : "FACE_HEAT_FLUX",
                "projection_variable"           : "PROJECTED_SCALAR1",
                "convection_variable"           : "CONVECTION_VELOCITY",
                "mesh_velocity_variable"        : "MESH_VELOCITY",
                "transfer_coefficient_variable" : "TRANSFER_COEFFICIENT",
                "velocity_variable"             : "VELOCITY",
                "specific_heat_variable"        : "SPECIFIC_HEAT",
                "reaction_variable"             : "REACTION_FLUX"
            },
            "time_stepping" : {
                "time_step": 1.0
            },
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
                "solver_type": "amgcl",
                "smoother_type":"ilu0",
                "krylov_type":"gmres",
                "coarsening_type":"aggregation",
                "max_iteration": 5000,
                "tolerance": 1e-9,
                "scaling": false
            },
            "element_replace_settings" : {
                "element_name" : "EulerianConvDiff",
                "condition_name" : "ThermalFace"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "auxiliary_variables_list" : []
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def AddVariables(self, target_model_part=None):

        if target_model_part == None:
            target_model_part = self.main_model_part

        ''' Add nodal solution step variables based on provided CONVECTION_DIFFUSION_SETTINGS
        '''
        convention_diffusion_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        density_variable = self.settings["convection_diffusion_variables"]["density_variable"].GetString()
        if (density_variable != ""):
            convention_diffusion_settings.SetDensityVariable(KratosMultiphysics.KratosGlobals.GetVariable(density_variable))
        diffusion_variable = self.settings["convection_diffusion_variables"]["diffusion_variable"].GetString()
        if (diffusion_variable != ""):
            convention_diffusion_settings.SetDiffusionVariable(KratosMultiphysics.KratosGlobals.GetVariable(diffusion_variable))
        unknown_variable = self.settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        if (unknown_variable != ""):
            convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable))
        volume_source_variable = self.settings["convection_diffusion_variables"]["volume_source_variable"].GetString()
        if (volume_source_variable != ""):
            convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable(volume_source_variable))
        surface_source_variable = self.settings["convection_diffusion_variables"]["surface_source_variable"].GetString()
        if (surface_source_variable != ""):
            convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable(surface_source_variable))
        projection_variable = self.settings["convection_diffusion_variables"]["projection_variable"].GetString()
        if (projection_variable != ""):
            convention_diffusion_settings.SetProjectionVariable(KratosMultiphysics.KratosGlobals.GetVariable(projection_variable))
        convection_variable = self.settings["convection_diffusion_variables"]["convection_variable"].GetString()
        if (convection_variable != ""):
            convention_diffusion_settings.SetConvectionVariable(KratosMultiphysics.KratosGlobals.GetVariable(convection_variable))
        mesh_velocity_variable = self.settings["convection_diffusion_variables"]["mesh_velocity_variable"].GetString()
        if (mesh_velocity_variable != ""):
            convention_diffusion_settings.SetMeshVelocityVariable(KratosMultiphysics.KratosGlobals.GetVariable(mesh_velocity_variable))
        transfer_coefficient_variable = self.settings["convection_diffusion_variables"]["transfer_coefficient_variable"].GetString()
        if (transfer_coefficient_variable != ""):
            convention_diffusion_settings.SetTransferCoefficientVariable(KratosMultiphysics.KratosGlobals.GetVariable(transfer_coefficient_variable))
        velocity_variable = self.settings["convection_diffusion_variables"]["velocity_variable"].GetString()
        if (velocity_variable != ""):
            convention_diffusion_settings.SetVelocityVariable(KratosMultiphysics.KratosGlobals.GetVariable(velocity_variable))
        specific_heat_variable = self.settings["convection_diffusion_variables"]["specific_heat_variable"].GetString()
        if (specific_heat_variable != ""):
            convention_diffusion_settings.SetSpecificHeatVariable(KratosMultiphysics.KratosGlobals.GetVariable(specific_heat_variable))
        reaction_variable = self.settings["convection_diffusion_variables"]["reaction_variable"].GetString()
        if (reaction_variable != ""):
            convention_diffusion_settings.SetReactionVariable(KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable))

        target_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convention_diffusion_settings)

        if target_model_part.ProcessInfo.Has(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS):
            if convention_diffusion_settings.IsDefinedDensityVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetDensityVariable())
            if convention_diffusion_settings.IsDefinedDiffusionVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetDiffusionVariable())
            if convention_diffusion_settings.IsDefinedUnknownVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetUnknownVariable())
            if convention_diffusion_settings.IsDefinedVolumeSourceVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetVolumeSourceVariable())
            if convention_diffusion_settings.IsDefinedSurfaceSourceVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetSurfaceSourceVariable())
            if convention_diffusion_settings.IsDefinedProjectionVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetProjectionVariable())
            if convention_diffusion_settings.IsDefinedConvectionVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetConvectionVariable())
            if convention_diffusion_settings.IsDefinedMeshVelocityVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetMeshVelocityVariable())
            if convention_diffusion_settings.IsDefinedTransferCoefficientVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetTransferCoefficientVariable())
            if convention_diffusion_settings.IsDefinedVelocityVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetVelocityVariable())
            if convention_diffusion_settings.IsDefinedSpecificHeatVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetSpecificHeatVariable())
            if convention_diffusion_settings.IsDefinedReactionVariable():
                target_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetReactionVariable())
        else:
            raise Exception("The provided target_model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

        # Adding nodal area variable (some solvers use it. TODO: Ask)
        #target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        # If LaplacianElement is used
        if (self.settings["element_replace_settings"]["element_name"].GetString() == "LaplacianElement"):
            target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def AddDofs(self):
        settings = self.main_model_part.ProcessInfo[KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS]
        if settings.IsDefinedReactionVariable():
            KratosMultiphysics.VariableUtils().AddDof(settings.GetUnknownVariable(), settings.GetReactionVariable(),self.main_model_part)
        else:
            KratosMultiphysics.VariableUtils().AddDof(settings.GetUnknownVariable(), self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "DOF's ADDED")

    def ImportModelPart(self):
        """This function imports the ModelPart"""
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

        if (self.settings["echo_level"].GetInt() > 0):
            KratosMultiphysics.Logger.PrintInfo(self.model)

        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part."""
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Initializing ...")
        # The convection_diffusion solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        convection_diffusion_solution_strategy = self.get_convection_diffusion_solution_strategy()
        convection_diffusion_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        if not self.is_restarted():
            convection_diffusion_solution_strategy.Initialize()
        else:
            # SetInitializePerformedFlag is not a member of SolvingStrategy but
            # is used by ResidualBasedNewtonRaphsonStrategy.
            try:
                convection_diffusion_solution_strategy.SetInitializePerformedFlag(True)
            except AttributeError:
                pass
        self.Check()
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Finished initialization.")

    def GetOutputVariables(self):
        pass

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        convection_diffusion_solution_strategy = self.get_convection_diffusion_solution_strategy()
        convection_diffusion_solution_strategy.Solve()

    def InitializeSolutionStep(self):
        self.get_convection_diffusion_solution_strategy().InitializeSolutionStep()

    def Predict(self):
        self.get_convection_diffusion_solution_strategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_convection_diffusion_solution_strategy().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.get_convection_diffusion_solution_strategy().FinalizeSolutionStep()

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
        self.get_convection_diffusion_solution_strategy().SetEchoLevel(level)

    def Clear(self):
        self.get_convection_diffusion_solution_strategy().Clear()

    def Check(self):
        self.get_convection_diffusion_solution_strategy().Check()

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

    def get_convection_diffusion_solution_strategy(self):
        if not hasattr(self, '_convection_diffusion_solution_strategy'):
            self._convection_diffusion_solution_strategy = self._create_convection_diffusion_solution_strategy()
        return self._convection_diffusion_solution_strategy

    def import_materials(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)

            # We set the properties that are nodal
            self._assign_nodally_properties()
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]

    #### Private functions ####

    def _assign_nodally_properties(self):

        # We transfer the values of the con.diff variables to the nodes
        with open(self.settings["material_import_settings"]["materials_filename"].GetString(), 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        for i in range(materials["properties"].size()):
            model_part = self.model.GetModelPart(materials["properties"][i]["model_part_name"].GetString())
            mat = materials["properties"][i]["Material"]

            for key, value in mat["Variables"].items():
                var = KratosMultiphysics.KratosGlobals.GetVariable(key)
                if (self._check_variable_to_set(var)):
                    if value.IsDouble():
                        KratosMultiphysics.VariableUtils().SetVariable(var, value.GetDouble(), model_part.Nodes)
                    elif value.IsVector():
                        KratosMultiphysics.VariableUtils().SetVariable(var, value.GetVector(), model_part.Nodes)
                    else:
                        raise ValueError("Type of value is not available")

    def _check_variable_to_set(self, var):
        thermal_settings = self.main_model_part.ProcessInfo[KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS]
        if (thermal_settings.IsDefinedDensityVariable()):
            if (thermal_settings.GetDensityVariable() == var):
                return True
        if (thermal_settings.IsDefinedDiffusionVariable()):
            if (thermal_settings.GetDiffusionVariable() == var):
                return True
        if (thermal_settings.IsDefinedVolumeSourceVariable()):
            if (thermal_settings.GetVolumeSourceVariable() == var):
                return True
        if (thermal_settings.IsDefinedSurfaceSourceVariable()):
            if (thermal_settings.GetSurfaceSourceVariable() == var):
                return True
        if (thermal_settings.IsDefinedProjectionVariable()):
            if (thermal_settings.GetProjectionVariable() == var):
                return True
        if (thermal_settings.IsDefinedConvectionVariable()):
            if (thermal_settings.GetConvectionVariable() == var):
                return True
        if (thermal_settings.IsDefinedTransferCoefficientVariable()):
            if (thermal_settings.GetTransferCoefficientVariable() == var):
                return True
        if (thermal_settings.IsDefinedSpecificHeatVariable()):
            if (thermal_settings.GetSpecificHeatVariable() == var):
                return True
        else:
            return False

    def _execute_after_reading(self):
        """Prepare computing model part and import constitutive laws."""
        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        # Assign mesh entities from domain and process sub model parts to the computing model part.
        check_and_prepare_model_process.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Import constitutive laws.
        materials_imported = self.import_materials()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were not imported.")

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information."""
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
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
        for _ in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            self.main_model_part.CloneTimeStep(time)

    def _get_element_condition_replace_settings(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## Elements
        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            for elem in self.main_model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break

        num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)

        element_name = self.settings["element_replace_settings"]["element_name"].GetString()

        if domain_size not in (2,3):
            raise Exception("DOMAIN_SIZE not set")

        if element_name == "EulerianConvDiff":
            if domain_size == 2:
                if num_nodes_elements == 3:
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D")
                else:
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D4N")
            else:
                if num_nodes_elements == 4:
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff3D")
                else:
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff3D8N")
        elif element_name in ("LaplacianElement","AdjointHeatDiffusionElement","QSConvectionDiffusionExplicit","DConvectionDiffusionExplicit"):
            name_string = "{0}{1}D{2}N".format(element_name,domain_size, num_nodes_elements)
            self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        ## Conditions
        num_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().SumAll(len(self.main_model_part.Conditions))

        if num_conditions > 0:
            num_nodes_conditions = 0
            if (len(self.main_model_part.Conditions) > 0):
                for cond in self.main_model_part.Conditions:
                    num_nodes_conditions = len(cond.GetNodes())
                    break

            num_nodes_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)

            condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
            if condition_name in ("FluxCondition","ThermalFace","Condition"):
                name_string = "{0}{1}D{2}N".format(condition_name,domain_size, num_nodes_conditions)
                self.settings["element_replace_settings"]["condition_name"].SetString(name_string)
        else:
            self.settings["element_replace_settings"]["condition_name"].SetString("")

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
        convergence_criterion = convergence_criteria_factory.ConvergenceCriteriaFactory(self._get_convergence_criterion_settings())
        return convergence_criterion.convergence_criterion

    def _create_linear_solver(self):
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["block_builder"].GetBool():
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    @classmethod
    def _create_solution_scheme(self):
        """Create the solution scheme for the convection-diffusion problem."""
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _create_convection_diffusion_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            convection_diffusion_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                convection_diffusion_solution_strategy = self._create_newton_raphson_strategy()
            else:
                convection_diffusion_solution_strategy = self._create_line_search_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return convection_diffusion_solution_strategy

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        convection_diffusion_scheme = self.get_solution_scheme()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              convection_diffusion_scheme,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        convection_diffusion_scheme = self.get_solution_scheme()
        convection_diffusion_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                        convection_diffusion_scheme,
                                        convection_diffusion_convergence_criterion,
                                        builder_and_solver,
                                        self.settings["max_iteration"].GetInt(),
                                        self.settings["compute_reactions"].GetBool(),
                                        self.settings["reform_dofs_at_each_step"].GetBool(),
                                        self.settings["move_mesh_flag"].GetBool())

    def _create_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        convection_diffusion_scheme = self.get_solution_scheme()
        convection_diffusion_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.LineSearchStrategy(computing_model_part,
                            convection_diffusion_scheme,
                            convection_diffusion_convergence_criterion,
                            builder_and_solver,
                            self.settings["max_iteration"].GetInt(),
                            self.settings["compute_reactions"].GetBool(),
                            self.settings["reform_dofs_at_each_step"].GetBool(),
                            self.settings["move_mesh_flag"].GetBool())

    def _ConvectionDiffusionVariablesCheck(self, custom_settings):
        default_settings = self.GetDefaultParameters()
        default_conv_diff_variables = default_settings["convection_diffusion_variables"]
        if not custom_settings.Has("convection_diffusion_variables"):
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "\'convection_diffusion_variables\' not defined, taking default ", default_conv_diff_variables)
        else:
            custom_conv_diff_variables = custom_settings["convection_diffusion_variables"]
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "density_variable", default_conv_diff_variables["density_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "diffusion_variable", default_conv_diff_variables["diffusion_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "unknown_variable", default_conv_diff_variables["unknown_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "volume_source_variable", default_conv_diff_variables["volume_source_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "surface_source_variable", default_conv_diff_variables["surface_source_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "projection_variable", default_conv_diff_variables["projection_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "convection_variable", default_conv_diff_variables["convection_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "mesh_velocity_variable", default_conv_diff_variables["mesh_velocity_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "transfer_coefficient_variable", default_conv_diff_variables["transfer_coefficient_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "velocity_variable", default_conv_diff_variables["velocity_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "specific_heat_variable", default_conv_diff_variables["specific_heat_variable"].GetString())
            self._ConvectionDiffusionSingleVariableCheck(custom_conv_diff_variables, "reaction_variable", default_conv_diff_variables["reaction_variable"].GetString())

    def _ConvectionDiffusionSingleVariableCheck(self, custom_conv_diff_variables, variable_entry, variable_name):
        if not custom_conv_diff_variables.Has(variable_entry):
            custom_conv_diff_variables.AddEmptyValue(variable_entry).SetString(variable_name)
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "\'{0}\' in \'convection_diffusion_variables\' not defined, taking default \'{1}\'.".format(variable_entry, variable_name))
