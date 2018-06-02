from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Other imports
import os


def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionBaseSolver(main_model_part, custom_settings)


class ConvectionDiffusionBaseSolver(object):
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
    _create_mechanical_solution_strategy
    _create_restart_utility

    The mechanical_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions get_mechanical_solution_strategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    settings -- Kratos parameters containing solver settings.
    main_model_part -- the model part used to construct the solver.
    """
    def __init__(self, main_model_part, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "SolverName - please provide a proper one",
            "echo_level": 0,
            "buffer_size": 2,
            "analysis_type": "non_linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "restart_settings" : {
                "load_restart"  : false,
                "save_restart"  : false
            },
            "computing_model_part_name" : "Thermal",
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
                "transfer_coefficient_variable" : "",
                "velocity_variable"             : "VELOCITY",
                "specific_heat_variable"        : "SPECIFIC_HEAT",
                "reaction_variable"             : "REACTION_FLUX"
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
            "processes_sub_model_part_list": [""],
            "auxiliary_variables_list" : []
        }
        """)

        # Adding warnings
        if not custom_settings.Has("convection_diffusion_variables"):
            self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "W-A-R-N-I-N-G: CONVECTION DIFFUSION  VARIABLES NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"])
        else:
            if not custom_settings["convection_diffusion_variables"].Has("density_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "W-A-R-N-I-N-G: DENSITY VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["density_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("diffusion_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "W-A-R-N-I-N-G: DIFUSSION VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["diffusion_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("unknown_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "W-A-R-N-I-N-G: UNKNOWN VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["unknown_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("volume_source_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "W-A-R-N-I-N-G: VOLUME SOURCE VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["volume_source_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("surface_source_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "W-A-R-N-I-N-G: SURFACE SOURCE VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["surface_source_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("projection_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", " W-A-R-N-I-N-G: PROJECTION VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["projection_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("convection_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", " W-A-R-N-I-N-G: CONVECTION VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["convection_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("mesh_velocity_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", " W-A-R-N-I-N-G: MESH VELOCITY VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["mesh_velocity_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("transfer_coefficient_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", " W-A-R-N-I-N-G: TRANSFER COEFFICIENT VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["transfer_coefficient_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("velocity_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", " W-A-R-N-I-N-G: VELOCITY VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["velocity_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("specific_heat_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", " W-A-R-N-I-N-G: SPECIFIC HEAT VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["specific_heat_variable"].GetString())
            if not custom_settings["convection_diffusion_variables"].Has("reaction_variable"):
                self.print_warning_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", " W-A-R-N-I-N-G: REACTION VARIABLE NOT DEFINED, TAKING DEFAULT", default_settings["convection_diffusion_variables"]["reaction_variable"].GetString())
    
        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part
        self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "Construction finished")

        # Set if the analysis is restarted
        if self.settings["restart_settings"].Has("load_restart"):
            load_restart = self.settings["restart_settings"]["load_restart"].GetBool()
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = load_restart
        else:
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def AddVariables(self):
        ''' Add nodal solution step variables based on provided CONVECTION_DIFFUSION_SETTINGS
        '''
        convention_diffusion_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        density_variable = self.settings["convection_diffusion_variables"]["density_variable"].GetString()
        if (density_variable is not ""):
            convention_diffusion_settings.SetDensityVariable(KratosMultiphysics.KratosGlobals.GetVariable(density_variable))
        diffusion_variable = self.settings["convection_diffusion_variables"]["diffusion_variable"].GetString()
        if (diffusion_variable is not ""):
            convention_diffusion_settings.SetDiffusionVariable(KratosMultiphysics.KratosGlobals.GetVariable(diffusion_variable))
        unknown_variable = self.settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        if (unknown_variable is not ""):
            convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable))
        volume_source_variable = self.settings["convection_diffusion_variables"]["volume_source_variable"].GetString()
        if (volume_source_variable is not ""):
            convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable(volume_source_variable))
        surface_source_variable = self.settings["convection_diffusion_variables"]["surface_source_variable"].GetString()
        if (surface_source_variable is not ""):
            convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable(surface_source_variable))
        projection_variable = self.settings["convection_diffusion_variables"]["projection_variable"].GetString()
        if (projection_variable is not ""):
            convention_diffusion_settings.SetProjectionVariable(KratosMultiphysics.KratosGlobals.GetVariable(projection_variable))
        convection_variable = self.settings["convection_diffusion_variables"]["convection_variable"].GetString()
        if (convection_variable is not ""):
            convention_diffusion_settings.SetConvectionVariable(KratosMultiphysics.KratosGlobals.GetVariable(convection_variable))
        mesh_velocity_variable = self.settings["convection_diffusion_variables"]["mesh_velocity_variable"].GetString()
        if (mesh_velocity_variable is not ""):
            convention_diffusion_settings.SetMeshVelocityVariable(KratosMultiphysics.KratosGlobals.GetVariable(mesh_velocity_variable))
        transfer_coefficient_variable = self.settings["convection_diffusion_variables"]["transfer_coefficient_variable"].GetString()
        if (transfer_coefficient_variable is not ""):
            convention_diffusion_settings.SetTransferCoefficientVariable(KratosMultiphysics.KratosGlobals.GetVariable(transfer_coefficient_variable))
        velocity_variable = self.settings["convection_diffusion_variables"]["velocity_variable"].GetString()
        if (velocity_variable is not ""):
            convention_diffusion_settings.SetVelocityVariable(KratosMultiphysics.KratosGlobals.GetVariable(velocity_variable))
        specific_heat_variable = self.settings["convection_diffusion_variables"]["specific_heat_variable"].GetString()
        if (specific_heat_variable is not ""):
            convention_diffusion_settings.SetSpecificHeatVariable(KratosMultiphysics.KratosGlobals.GetVariable(specific_heat_variable))
        reaction_variable = self.settings["convection_diffusion_variables"]["reaction_variable"].GetString()
        if (reaction_variable is not ""):
            convention_diffusion_settings.SetReactionVariable(KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable))
            
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convention_diffusion_settings)
        
        if self.main_model_part.ProcessInfo.Has(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS):
            if convention_diffusion_settings.IsDefinedDensityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetDensityVariable())
            if convention_diffusion_settings.IsDefinedDiffusionVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetDiffusionVariable())
            if convention_diffusion_settings.IsDefinedUnknownVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetUnknownVariable())
            if convention_diffusion_settings.IsDefinedVolumeSourceVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetVolumeSourceVariable())
            if convention_diffusion_settings.IsDefinedSurfaceSourceVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetSurfaceSourceVariable())
            if convention_diffusion_settings.IsDefinedProjectionVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetProjectionVariable())
            if convention_diffusion_settings.IsDefinedConvectionVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetConvectionVariable())
            if convention_diffusion_settings.IsDefinedMeshVelocityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetMeshVelocityVariable())
            if convention_diffusion_settings.IsDefinedTransferCoefficientVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetTransferCoefficientVariable())
            if convention_diffusion_settings.IsDefinedVelocityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetVelocityVariable())
            if convention_diffusion_settings.IsDefinedSpecificHeatVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetSpecificHeatVariable())
            if convention_diffusion_settings.IsDefinedReactionVariable():
                self.main_model_part.AddNodalSolutionStepVariable(convention_diffusion_settings.GetReactionVariable())
        else:
            raise Exception("The provided main_model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")
        
        # Adding nodal area variable (some solvers use it. TODO: Ask)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        # If LaplacianElement is used
        if (self.settings["element_replace_settings"]["element_name"].GetString() == "LaplacianElement"):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        
        self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "Variables ADDED")


    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        # this can safely be called also for restarts, it is internally checked if the dofs exist already
        settings = self.main_model_part.ProcessInfo[KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS]                     
        if settings.IsDefinedReactionVariable():
            KratosMultiphysics.VariableUtils().AddDof(settings.GetUnknownVariable(), settings.GetReactionVariable(),self.main_model_part)
        else:
            KratosMultiphysics.VariableUtils().AddDof(settings.GetUnknownVariable(), self.main_model_part)
        self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "DOF's ADDED")

    def ImportModelPart(self):
        """ Legacy function, use ReadModelPart and PrepareModelPartForSolver instead """
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Importing model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
        if self.is_restarted():
            self.get_restart_utility().LoadRestart()
        elif(self.settings["model_import_settings"]["input_type"].GetString() == "use_input_model_part"): #TODO: change this once a proper way is agreed
            self.PrepareModelPartForSolver() 
            
        else:
            if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
                # Import model part from mdpa file.
                KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
                KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
                KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Finished reading model part from mdpa file.")
                self.PrepareModelPartForSolver()
                    
            else:
                raise Exception("Other model part input options are not yet implemented.")
        KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Finished importing model part.")

    def ReadModelPart(self):
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Reading model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
        if self.is_restarted():
            self.get_restart_utility().LoadRestart()
        elif(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # Import model part from mdpa file.
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "Finished reading model part from mdpa file.")
        else:
            raise Exception("Other model part input options are not yet implemented.")
        KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Finished reading model part.")

    def PrepareModelPartForSolver(self):
            
        if not self.is_restarted():
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()

            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors).Execute()
            
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()
        
            self._set_and_fill_buffer()
        
        self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]::", "ModelPart prepared for Solver.")

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "Initializing ...")
        # The mechanical solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solution_strategy = self.get_mechanical_solution_strategy()
        mechanical_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        if not self.is_restarted():
            mechanical_solution_strategy.Initialize()
        else:
            # SetInitializePerformedFlag is not a member of SolvingStrategy but
            # is used by ResidualBasedNewtonRaphsonStrategy.
            try:
                mechanical_solution_strategy.SetInitializePerformedFlag(True)
            except AttributeError:
                pass
        self.Check()
        self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "Finished initialization.")

    def GetComputingModelPart(self):
        return self.main_model_part #.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    def GetOutputVariables(self):
        pass

    def SaveRestart(self):
        # Check could be integrated in the utility
        # It is here intentionally, this way the utility is only created if it is actually needed!
        if self.settings["restart_settings"].Has("save_restart"):
            if (self.settings["restart_settings"]["save_restart"].GetBool() == True):
                # the check if this step is a restart-output step is done internally
                self.get_restart_utility().SaveRestart()

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solution_strategy = self.get_mechanical_solution_strategy()
        mechanical_solution_strategy.Solve()

    def InitializeSolutionStep(self):
        self.get_mechanical_solution_strategy().InitializeSolutionStep()

    def Predict(self):
        self.get_mechanical_solution_strategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_mechanical_solution_strategy().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.get_mechanical_solution_strategy().FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        return self.delta_time

    def SetDeltaTime(self, dt):
        # This is a TEMPORARY function until the solver can compute dt!
        self.delta_time = dt

    def SetEchoLevel(self, level):
        self.get_mechanical_solution_strategy().SetEchoLevel(level)

    def Clear(self):
        self.get_mechanical_solution_strategy().Clear()

    def Check(self):
        self.get_mechanical_solution_strategy().Check()

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

    def get_mechanical_solution_strategy(self):
        if not hasattr(self, '_mechanical_solution_strategy'):
            self._mechanical_solution_strategy = self._create_mechanical_solution_strategy()
        return self._mechanical_solution_strategy

    def get_restart_utility(self):
        if not hasattr(self, '_restart_utility'):
            self._restart_utility = self._create_restart_utility()
        return self._restart_utility

    def import_materials(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            import read_materials_process
            # Create a dictionary of model parts.
            Model = KratosMultiphysics.Model()
            Model.AddModelPart(self.main_model_part)
            # Add constitutive laws and material properties from json file to model parts.
            read_materials_process.ReadMaterialsProcess(Model, self.settings["material_import_settings"])
            
            # We set the properties that are nodal
            self._assign_nodally_properties()
            
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def validate_and_transfer_matching_settings(self, origin_settings, destination_settings):
        """Transfer matching settings from origin to destination.

        If a name in origin matches a name in destination, then the setting is
        validated against the destination.

        The typical use is for validating and extracting settings in derived classes:

        class A:
            def __init__(self, model_part, a_settings):
                default_a_settings = Parameters('''{
                    ...
                }''')
                a_settings.ValidateAndAssignDefaults(default_a_settings)
        class B(A):
            def __init__(self, model_part, custom_settings):
                b_settings = Parameters('''{
                    ...
                }''') # Here the settings contain default values.
                self.validate_and_transfer_matching_settings(custom_settings, b_settings)
                super().__init__(model_part, custom_settings)
        """
        for name, dest_value in destination_settings.items():
            if origin_settings.Has(name): # Validate and transfer value.
                orig_value = origin_settings[name]
                if dest_value.IsDouble() and orig_value.IsDouble():
                    destination_settings[name].SetDouble(origin_settings[name].GetDouble())
                elif dest_value.IsInt() and orig_value.IsInt():
                    destination_settings[name].SetInt(origin_settings[name].GetInt())
                elif dest_value.IsBool() and orig_value.IsBool():
                    destination_settings[name].SetBool(origin_settings[name].GetBool())
                elif dest_value.IsString() and orig_value.IsString():
                    destination_settings[name].SetString(origin_settings[name].GetString())
                elif dest_value.IsArray() and orig_value.IsArray():
                    if dest_value.size() != orig_value.size():
                        raise Exception('len("' + name + '") != ' + str(dest_value.size()))
                    for i in range(dest_value.size()):
                        if dest_value[i].IsDouble() and orig_value[i].IsDouble():
                            dest_value[i].SetDouble(orig_value[i].GetDouble())
                        elif dest_value[i].IsInt() and orig_value[i].IsInt():
                            dest_value[i].SetInt(orig_value[i].GetInt())
                        elif dest_value[i].IsBool() and orig_value[i].IsBool():
                            dest_value[i].SetBool(orig_value[i].GetBool())
                        elif dest_value[i].IsString() and orig_value[i].IsString():
                            dest_value[i].SetString(orig_value[i].GetString())
                        elif dest_value[i].IsSubParameter() and orig_value[i].IsSubParameter():
                            self.validate_and_transfer_matching_settings(orig_value[i], dest_value[i])
                            if len(orig_value[i].items()) != 0:
                                raise Exception('Json settings not found in default settings: ' + orig_value[i].PrettyPrintJsonString())
                        else:
                            raise Exception('Unsupported parameter type.')
                elif dest_value.IsSubParameter() and orig_value.IsSubParameter():
                    self.validate_and_transfer_matching_settings(orig_value, dest_value)
                    if len(orig_value.items()) != 0:
                        raise Exception('Json settings not found in default settings: ' + orig_value.PrettyPrintJsonString())
                else:
                    raise Exception('Unsupported parameter type.')
                origin_settings.RemoveValue(name)

    def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]

    def print_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))

    def print_warning_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KratosMultiphysics.Logger.PrintWarning(" ".join(map(str,args)))

    #### Private functions ####

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
        """Prepare computing model part and import constitutive laws. """
        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        # Assign mesh entities from domain and process sub model parts to the computing model part.
        import check_and_prepare_model_process_convection_diffusion as check_and_prepare_model_process
        check_and_prepare_model_process.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Import constitutive laws.
        materials_imported = self.import_materials()
        if materials_imported:
            self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "Materials were successfully imported.")
        else:
            self.print_on_rank_zero("::[ConvectionDiffusionBaseSolver]:: ", "Materials were not imported.")

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

    def _get_restart_settings(self):
        restart_settings = self.settings["restart_settings"].Clone()
        restart_settings.AddValue("input_filename", self.settings["model_import_settings"]["input_filename"])
        restart_settings.AddValue("echo_level", self.settings["echo_level"])
        restart_settings.RemoveValue("load_restart")
        restart_settings.RemoveValue("save_restart")

        return restart_settings

    def _get_element_condition_replace_settings(self):
        # Duplicate model part
        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            num_nodes_elements = len(self.main_model_part.Elements[1].GetNodes())

        ## Elements
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
                if (num_nodes_elements == 3):
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D")
                else:
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff2D4N")
            elif (self.settings["element_replace_settings"]["element_name"].GetString() == "LaplacianElement"):
                self.settings["element_replace_settings"]["element_name"].SetString("LaplacianElement2D3N")
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
                if (num_nodes_elements == 4):
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff3D")
                else:
                    self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff3D8N")
            elif (self.settings["element_replace_settings"]["element_name"].GetString() == "LaplacianElement"):
                if (num_nodes_elements == 4):
                    self.settings["element_replace_settings"]["element_name"].SetString("LaplacianElement3D4N")
                elif (num_nodes_elements == 8):
                    self.settings["element_replace_settings"]["element_name"].SetString("LaplacianElement3D8N")
                else:
                    self.settings["element_replace_settings"]["element_name"].SetString("LaplacianElement3D27N")
        else:
            raise Exception("DOMAIN_SIZE not set")
        
        ## Conditions
        num_nodes_conditions = 0
        if (len(self.main_model_part.Conditions) > 0):
            num_nodes_conditions = len(self.main_model_part.Conditions[1].GetNodes())
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            if (self.settings["element_replace_settings"]["condition_name"].GetString() == "FluxCondition"):
                self.settings["element_replace_settings"]["condition_name"].SetString("FluxCondition2D2N")
            elif (self.settings["element_replace_settings"]["condition_name"].GetString() == "ThermalFace"):
                self.settings["element_replace_settings"]["condition_name"].SetString("ThermalFace2D")
            else:
                self.settings["element_replace_settings"]["condition_name"].SetString("LineCondition2D2N")
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            if (self.settings["element_replace_settings"]["condition_name"].GetString() == "FluxCondition"):
                if (num_nodes_conditions == 3):
                    self.settings["element_replace_settings"]["condition_name"].SetString("FluxCondition3D3N")
                else:
                    self.settings["element_replace_settings"]["condition_name"].SetString("FluxCondition3D4N")
            elif (self.settings["element_replace_settings"]["condition_name"].GetString() == "ThermalFace"):
                self.settings["element_replace_settings"]["condition_name"].SetString("ThermalFace3D")
            else:
                if (num_nodes_conditions == 3):
                    self.settings["element_replace_settings"]["condition_name"].SetString("SurfaceCondition3D3N")
                else:
                    self.settings["element_replace_settings"]["condition_name"].SetString("SurfaceCondition3D4N")
        else:
            raise Exception("DOMAIN_SIZE not set")

        #modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        #modeler.GenerateModelPart(self.main_model_part, self.GetComputingModelPart(), self.settings["element_replace_settings"]["element_name"].GetString(), self.settings["element_replace_settings"]["condition_name"].GetString())

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
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

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
        raise Exception("Solution Scheme creation must be implemented in the derived class.")
        

    def _create_mechanical_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            mechanical_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                mechanical_solution_strategy = self._create_newton_raphson_strategy()
            else:
                mechanical_solution_strategy = self._create_line_search_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              mechanical_scheme,
                                                              linear_solver,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                        mechanical_scheme,
                                        linear_solver,
                                        mechanical_convergence_criterion,
                                        builder_and_solver,
                                        self.settings["max_iteration"].GetInt(),
                                        self.settings["compute_reactions"].GetBool(),
                                        self.settings["reform_dofs_at_each_step"].GetBool(),
                                        self.settings["move_mesh_flag"].GetBool())

    def _create_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.LineSearchStrategy(computing_model_part,
                            mechanical_scheme,
                            linear_solver,
                            mechanical_convergence_criterion,
                            builder_and_solver,
                            self.settings["max_iteration"].GetInt(),
                            self.settings["compute_reactions"].GetBool(),
                            self.settings["reform_dofs_at_each_step"].GetBool(),
                            self.settings["move_mesh_flag"].GetBool())

    def _create_restart_utility(self):
        """Create the restart utility. Has to be overridden for MPI/trilinos-solvers"""
        import restart_utility
        rest_utility = restart_utility.RestartUtility(self.main_model_part,
                                                      self._get_restart_settings())
        return rest_utility
