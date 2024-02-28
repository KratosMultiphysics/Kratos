# Importing Python modules
import numpy
import numpy.linalg

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from  KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
have_mesh_moving = CheckIfApplicationsAvailable("MeshMovingApplication")
if have_mesh_moving:
    import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.FluidDynamicsApplication import read_distance_from_file

class EmbeddedFormulation(object):
    """Helper class to define embedded-dependent parameters."""
    def __init__(self, formulation_settings):
        self.element_name = None
        self.condition_name = None
        self.process_info_data = {}
        self.element_has_nodal_properties = False
        self.historical_nodal_properties_variables_list = []
        self.non_historical_nodal_properties_variables_list = []

        if formulation_settings.Has("element_type"):
            element_type = formulation_settings["element_type"].GetString()
            if element_type == "embedded_navier_stokes":
                self._SetUpClassicEmbeddedNavierStokes(formulation_settings)
            elif element_type == "embedded_symbolic_navier_stokes":
                warn_msg  = 'Provided \'element_name\' is \'embedded_symbolic_navier_stokes\'. This has been renamed to \'embedded_weakly_compressible_navier_stokes\'. Use this instead.'
                KratosMultiphysics.Logger.PrintWarning(warn_msg)
                self._SetUpEmbeddedWeaklyCompressibleNavierStokes(formulation_settings)
            elif element_type == "embedded_weakly_compressible_navier_stokes":
                self._SetUpEmbeddedWeaklyCompressibleNavierStokes(formulation_settings)
            elif element_type == "embedded_ausas_navier_stokes":
                self._SetUpClassicEmbeddedAusasNavierStokes(formulation_settings)
            elif element_type == "embedded_symbolic_navier_stokes_discontinuous":
                warn_msg  = 'Provided \'element_name\' is \'embedded_symbolic_navier_stokes_discontinuous\'. This has been renamed to \'embedded_weakly_compressible_navier_stokes_discontinuous\'. Use this instead.'
                KratosMultiphysics.Logger.PrintWarning(warn_msg)
                self._SetUpEmbeddedWeaklyCompressibleNavierStokesDiscontinuous(formulation_settings)
            elif element_type == "embedded_weakly_compressible_navier_stokes_discontinuous":
                self._SetUpEmbeddedWeaklyCompressibleNavierStokesDiscontinuous(formulation_settings)
        else:
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

    def SetProcessInfo(self, model_part):
        for variable,value in self.process_info_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpClassicEmbeddedNavierStokes(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_navier_stokes",
            "is_slip": false,
            "slip_length": 1.0e8,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0,
            "level_set_type": "continuous"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = False

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        if formulation_settings["is_slip"].GetBool():
            self.process_info_data[KratosCFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()

    def _SetUpEmbeddedWeaklyCompressibleNavierStokes(self, formulation_settings):
        #TODO: Remove this after deprecation period is over
        if (formulation_settings["element_type"].GetString() == "embedded_symbolic_navier_stokes"):
            formulation_settings["element_type"].SetString("embedded_weakly_compressible_navier_stokes")

        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_weakly_compressible_navier_stokes",
            "is_slip": false,
            "slip_length": 1.0e8,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0,
            "level_set_type": "continuous"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedWeaklyCompressibleNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.DENSITY]
        self.non_historical_nodal_properties_variables_list = [KratosMultiphysics.SOUND_VELOCITY]

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        if formulation_settings["is_slip"].GetBool():
            self.process_info_data[KratosCFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()

    def _SetUpClassicEmbeddedAusasNavierStokes(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_ausas_navier_stokes",
            "is_slip": true,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0,
            "level_set_type": "discontinuous"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedAusasNavierStokes"
        self.condition_name = "EmbeddedAusasNavierStokesWallCondition"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = False

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()

    def _SetUpEmbeddedWeaklyCompressibleNavierStokesDiscontinuous(self, formulation_settings):
        #TODO: Remove this after deprecation period is over
        if (formulation_settings["element_type"].GetString() == "embedded_symbolic_navier_stokes_discontinuous"):
            formulation_settings["element_type"].SetString("embedded_weakly_compressible_navier_stokes_discontinuous")

        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_weakly_compressible_navier_stokes_discontinuous",
            "is_slip": true,
            "slip_length": 1.0e8,
            "penalty_coefficient": 0.1,
            "dynamic_tau": 1.0,
            "level_set_type": "discontinuous"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedWeaklyCompressibleNavierStokesDiscontinuous"
        self.condition_name = "NavierStokesWallCondition"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.DENSITY]
        self.non_historical_nodal_properties_variables_list = [KratosMultiphysics.SOUND_VELOCITY]

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        if formulation_settings["is_slip"].GetBool():
            self.process_info_data[KratosCFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()
        else:
            self.process_info_data[KratosCFD.SLIP_LENGTH] = 0.0


def CreateSolver(model, custom_settings):
    return NavierStokesEmbeddedMonolithicSolver(model, custom_settings)

class NavierStokesEmbeddedMonolithicSolver(FluidSolver):

    def __GetDistanceModificationDefaultSettings(self, level_set_type):
        if level_set_type == "continuous":
            return self.__GetContinuousDistanceModificationDefaultSettings()
        elif level_set_type == "discontinuous":
            return self.__GetDiscontinuousDistanceModificationDefaultSettings()
        else:
            err_msg = 'Provided level set type is: \'' + level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    @classmethod
    def __GetContinuousDistanceModificationDefaultSettings(cls):
        return KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "distance_threshold": 1e-3,
            "continuous_distance": true,
            "check_at_each_time_step": true,
            "avoid_almost_empty_elements": true,
            "deactivate_full_negative_elements": true
        }''')

    @classmethod
    def __GetDiscontinuousDistanceModificationDefaultSettings(cls):
        return KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "distance_threshold": 1e-3,
            "continuous_distance": false,
            "check_at_each_time_step": true,
            "avoid_almost_empty_elements": false,
            "deactivate_full_negative_elements": false
        }''')

    @classmethod
    def _get_fm_ale_implicit_default_settings(cls):
        return KratosMultiphysics.Parameters("""
        {
            "virtual_model_part_name": "VirtualModelPart",
            "structure_model_part_name": "",
            "linear_solver_settings": {
                "solver_type": "cg",
                "tolerance": 1.0e-8,
                "max_iteration": 1000
            },
            "embedded_nodal_variable_settings": {
                "gradient_penalty_coefficient": 1.0e-3,
                "linear_solver_settings": {
                    "preconditioner_type": "amg",
                    "solver_type": "amgcl",
                    "smoother_type": "ilu0",
                    "krylov_type": "cg",
                    "max_iteration": 1000,
                    "verbosity": 0,
                    "tolerance": 1e-8,
                    "scaling": false,
                    "block_size": 1,
                    "use_block_matrices_if_possible": true
                }
            }
        }
        """)

    @classmethod
    def _get_fm_ale_explicit_default_settings(cls):
        return KratosMultiphysics.Parameters("""
        {
            "virtual_model_part_name": "VirtualModelPart",
            "structure_model_part_name": "",
            "search_radius": 0.0
        }
        """)

    def _get_fm_ale_solver_default_settings(self, mesh_movement):
        if mesh_movement == "implicit":
            return self._get_fm_ale_implicit_default_settings()
        elif mesh_movement == "explicit":
            return self._get_fm_ale_explicit_default_settings()
        else:
            raise Exception("Provided mesh movement \'" + mesh_movement + "\'. Available options are \'implicit\' and \'explicit\'.")

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "embedded_solver_from_defaults",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
            },
            "distance_modification_settings": {
            },
            "maximum_iterations": 7,
            "echo_level": 0,
            "time_order": 2,
            "time_scheme": "bdf2",
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0
            },
            "move_mesh_flag": false,
            "formulation": {
                "element_type": "embedded_element_from_defaults",
                "dynamic_tau": 1.0
            },
            "fm_ale_settings": {
                "fm_ale_step_frequency": 0,
                "mesh_movement": "implicit",
                "embedded_velocity_calculation": "from_fluid_mesh_velocity",
                "rbf_interpolation_search_radius": 0.0,
                "rbf_interpolation_edge_factor": 2.0,
                "fm_ale_solver_settings": {
                }
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesEmbeddedMonolithicSolver, cls).GetDefaultParameters())
        return default_settings

    def ValidateSettings(self):
        """Overriding python_solver ValidateSettings to validate the fm_ale_settings
        """
        super(NavierStokesEmbeddedMonolithicSolver, self).ValidateSettings()
        self.settings["fm_ale_settings"].ValidateAndAssignDefaults(self.GetDefaultParameters()["fm_ale_settings"])
        if self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() > 0:
            mesh_movement = self.settings["fm_ale_settings"]["mesh_movement"].GetString()
            self.settings["fm_ale_settings"]["fm_ale_solver_settings"].ValidateAndAssignDefaults(self._get_fm_ale_solver_default_settings(mesh_movement))

    def __init__(self, model, custom_settings):
        # TODO: DO SOMETHING IN HERE TO REMOVE THE "time_order" FROM THE DEFAULT SETTINGS BUT KEEPING THE BACKWARDS COMPATIBILITY
        super(NavierStokesEmbeddedMonolithicSolver,self).__init__(model,custom_settings)

        self.min_buffer_size = 3
        self.embedded_formulation = EmbeddedFormulation(self.settings["formulation"])
        self.element_name = self.embedded_formulation.element_name
        self.condition_name = self.embedded_formulation.condition_name
        self.level_set_type = self.embedded_formulation.level_set_type
        self.element_integrates_in_time = self.embedded_formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.embedded_formulation.element_has_nodal_properties
        self.historical_nodal_properties_variables_list = self.embedded_formulation.historical_nodal_properties_variables_list
        self.non_historical_nodal_properties_variables_list = self.embedded_formulation.non_historical_nodal_properties_variables_list

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        # If the FM-ALE is required, do a first call to the __GetFmAleVirtualModelPart
        # Note that this will create the virtual model part in the model
        if self._FmAleIsActive():
            self.__GetFmAleVirtualModelPart()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesEmbeddedMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)              # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_PRESSURE)          # Post-process variable (stores the fluid nodes pressure and is set to 0 in the structure ones)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_VELOCITY)          # Post-process variable (stores the fluid nodes velocity and is set to 0 in the structure ones)

        # Adding variables required for the FM-ALE algorithm
        if self._FmAleIsActive():
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_REACTION)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            for variable in self.historical_nodal_properties_variables_list:
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def AddDofs(self):
        # Add formulation DOFs and reactions
        super().AddDofs()

        # Add mesh motion problem DOFs for the FM-ALE algorithm
        if self._FmAleIsActive():
            dofs_and_reactions_to_add = []
            dofs_and_reactions_to_add.append(["MESH_DISPLACEMENT_X", "MESH_REACTION_X"])
            dofs_and_reactions_to_add.append(["MESH_DISPLACEMENT_Y", "MESH_REACTION_Y"])
            dofs_and_reactions_to_add.append(["MESH_DISPLACEMENT_Z", "MESH_REACTION_Z"])
            KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)

            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "FM-ALE DOFs added correctly.")

    def PrepareModelPart(self):
        # Call the base solver PrepareModelPart()
        super(NavierStokesEmbeddedMonolithicSolver, self).PrepareModelPart()

        # Set the extra requirements of the embedded formulation
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Sets the embedded formulation configuration
            self.__SetEmbeddedFormulation()
            ## Setting the nodal distance
            self.__SetDistanceFunction()

    def Initialize(self):
        # If the solver requires an instance of the stabilized embedded_formulation class, set the process info variables
        if hasattr(self, 'embedded_formulation'):
            self.embedded_formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        # Set the distance modification process
        self.GetDistanceModificationProcess().ExecuteInitialize()

        # For the primitive Ausas formulation, set the find nodal neighbours process
        # Recall that the Ausas condition requires the nodal neighbours.
        if (self.settings["formulation"]["element_type"].GetString() == "embedded_ausas_navier_stokes"):
            computing_model_part = self.GetComputingModelPart()
            data_communicator = computing_model_part.GetCommunicator().GetDataCommunicator()
            self.find_nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(
                data_communicator,
                computing_model_part)

        # If required, intialize the FM-ALE utility
        if self._FmAleIsActive():
            self.fm_ale_step = 1
            # Fill the virtual model part geometry. Note that the mesh moving util is created in this first call
            self.__GetFmAleUtility().Initialize(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def AdvanceInTime(self, current_time):
        # Call base solver AdvanceInTime to clone the time step and get the new time
        new_time = super(NavierStokesEmbeddedMonolithicSolver, self).AdvanceInTime(current_time)

        # Save the current step and time in the virtual model part process info
        if self._FmAleIsActive():
            self.__GetFmAleVirtualModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1
            self.__GetFmAleVirtualModelPart().ProcessInfo[KratosMultiphysics.TIME] = new_time

        return new_time

    def InitializeSolutionStep(self):
        # Compute the BDF coefficients
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        # If required, compute the nodal neighbours
        if (self.settings["formulation"]["element_type"].GetString() == "embedded_ausas_navier_stokes"):
            (self.find_nodal_neighbours_process).Execute()

        # Set the virtual mesh values from the background mesh
        self.__SetVirtualMeshValues()

        # Call the base solver InitializeSolutionStep()
        super(NavierStokesEmbeddedMonolithicSolver, self).InitializeSolutionStep()

    def SolveSolutionStep(self):
        # Correct the distance field
        # Note that this is intentionally placed in here (and not in the InitializeSolutionStep() of the solver
        # It has to be done before each call to the Solve() in case an outer non-linear iteration is performed (FSI)
        self.GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Perform the FM-ALE operations
        # Note that this also sets the EMBEDDED_VELOCITY from the MESH_VELOCITY
        self.__DoFmAleOperations()

        # Call the base SolveSolutionStep to solve the embedded CFD problem
        is_converged = super(NavierStokesEmbeddedMonolithicSolver,self).SolveSolutionStep()

        # Undo the FM-ALE virtual mesh movement
        self.__UndoFMALEOperations()

        # Restore the fluid node fixity to its original status
        # Note that this is intentionally placed in here (and not in the FinalizeSolutionStep() of the solver
        # It has to be done after each call to the Solve() and the FM-ALE in case an outer non-linear iteration is performed (FSI)
        self.GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

        return is_converged

    def FinalizeSolutionStep(self):
        # Call the base solver FinalizeSolutionStep()
        super(NavierStokesEmbeddedMonolithicSolver, self).FinalizeSolutionStep()

        # Do the FM-ALE end of step operations
        self.__UpdateFMALEStepCounter()

    #TODO: THIS COULD BE SAFELY REMOVED ONCE WE OLD EMBEDDED ELEMENTS ARE REMOVED
    def _SetPhysicalProperties(self):
        materials_imported = super()._SetPhysicalProperties()

        #TODO: REMOVE THIS ONCE WE REMOVE THE OLD EMBEDDED ELEMENTS
        #TODO: THE SOUND_VELOCITY MUST BE ALWAYS RETRIEVED FROM THE PROPERTIES OR THE NODES AS THE NEW WEAKLY COMPRESSIBLE ELEMENT DO
        for el in self.main_model_part.Elements:
            if el.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
                sound_velocity = el.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
            else:
                sound_velocity = 1.0e+12 # Default sound velocity value
                KratosMultiphysics.Logger.PrintWarning('No \'SOUND_VELOCITY\' value found in Properties {0}. Setting default value {1}'.format(el.Properties.Id, sound_velocity))
            if sound_velocity <= 0.0:
                raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")
        self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = sound_velocity

        return materials_imported

    def _SetNodalProperties(self):
        set_density = KratosMultiphysics.DENSITY in self.historical_nodal_properties_variables_list
        set_sound_velocity = KratosMultiphysics.SOUND_VELOCITY in self.non_historical_nodal_properties_variables_list

        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            # Get DENSITY from properties
            if set_density:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                if rho <= 0.0:
                    raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            # Get SOUND_VELOCITY
            if set_sound_velocity:
                if el.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
                    sound_velocity = el.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
                else:
                    sound_velocity = 1.0e+12 # Default sound velocity value
                    KratosMultiphysics.Logger.PrintWarning('No \'SOUND_VELOCITY\' value found in Properties {0}. Setting default value {1}'.format(el.Properties.Id, sound_velocity))
                if sound_velocity <= 0.0:
                    raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        if set_density:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        if set_sound_velocity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)

    def __SetEmbeddedFormulation(self):
        # Set the SLIP elemental flag
        if (self.settings["formulation"]["is_slip"].GetBool()):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.GetComputingModelPart().Elements)
        else:
            # Set the SLIP elemental flag to false in the entire domain
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, False, self.GetComputingModelPart().Elements)

        # Save the formulation settings in the ProcessInfo
        self.embedded_formulation.SetProcessInfo(self.main_model_part)

    def __SetDistanceFunction(self):
        ## Set the nodal distance function
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            DistanceUtility = read_distance_from_file.DistanceImportUtility(self.main_model_part, self.settings["distance_reading_settings"])
            DistanceUtility.ImportDistance()
        elif (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            KratosMultiphysics.Logger.PrintInfo("Navier Stokes Embedded Solver","Distance function taken from the .mdpa input file.")
            # Recall to swap the distance sign (GiD considers d<0 in the fluid region)
            for node in self.main_model_part.Nodes:
                distance_value = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -distance_value)

    def GetDistanceModificationProcess(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self.__CreateDistanceModificationProcess()
        return self._distance_modification_process

    def __CreateDistanceModificationProcess(self):
        # Set the distance modification settings according to the level set type
        # Note that the distance modification process is applied to the volume model part
        distance_modification_settings = self.settings["distance_modification_settings"]
        distance_modification_settings.ValidateAndAssignDefaults(self.__GetDistanceModificationDefaultSettings(self.level_set_type))
        aux_full_volume_part_name = self.settings["model_part_name"].GetString() + "." + self.settings["volume_model_part_name"].GetString()
        distance_modification_settings["model_part_name"].SetString(aux_full_volume_part_name)
        return KratosCFD.DistanceModificationProcess(self.model, distance_modification_settings)

    def _FmAleIsActive(self):
        return self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() > 0

    def __GetFmAleStructureModelPart(self):
        structure_model_part_name = self.settings["fm_ale_settings"]["fm_ale_solver_settings"]["structure_model_part_name"].GetString()
        if self.model.HasModelPart(structure_model_part_name):
            return self.model.GetModelPart(structure_model_part_name)
        else:
            raise Exception("Structure model part {0} not found in model.\n It is expected to be added in your custom analysis stage file.".format(structure_model_part_name))

    def __GetFmAleVirtualModelPart(self):
        if not hasattr(self, '_virtual_model_part'):
            self._virtual_model_part = self.__CreateFmAleVirtualModelPart()
        return self._virtual_model_part

    def __CreateFmAleVirtualModelPart(self):
        virtual_model_part_name = self.settings["fm_ale_settings"]["fm_ale_solver_settings"]["virtual_model_part_name"].GetString()
        virtual_model_part = self.model.CreateModelPart(virtual_model_part_name)
        return virtual_model_part

    def __GetFmAleUtility(self):
        if not hasattr (self, '_mesh_moving_util'):
            self._mesh_moving_util = self.__CreateFmAleUtility()
        return self._mesh_moving_util

    def __CreateFmAleUtility(self):
        if have_mesh_moving:
            mesh_movement = self.settings["fm_ale_settings"]["mesh_movement"].GetString()
            if (mesh_movement == "implicit"):
                mesh_moving_util = KratosMeshMoving.FixedMeshALEUtilities(
                    self.model,
                    self.settings["fm_ale_settings"]["fm_ale_solver_settings"])
            else:
                raise Exception("FM-ALE mesh_movement set to \'" + mesh_movement + "\'. Available option is \'implicit\'.")

            return mesh_moving_util
        else:
            raise Exception("MeshMovingApplication is required to construct the FM-ALE utility (ExplicitFixedMeshALEUtilities)")

    def __IsFmAleStep(self):
        if self._FmAleIsActive():
            if (self.fm_ale_step == self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt()):
                return True
            else:
                return False
        else:
            return False

    def __UpdateFMALEStepCounter(self):
        if self._FmAleIsActive():
            if (self.__IsFmAleStep()):
                # Reset the FM-ALE steps counter
                self.fm_ale_step = 1
            else:
                # Update the FM-ALE steps counter
                self.fm_ale_step += 1

    def __SetVirtualMeshValues(self):
        if self.__IsFmAleStep():
            # Fill the virtual model part variable values: VELOCITY (n,nn), PRESSURE (n,nn)
            self.__GetFmAleUtility().SetVirtualMeshValuesFromOriginMesh()

    def __DoFmAleOperations(self):
        if self.__IsFmAleStep():
            # Solve the mesh problem
            delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            self.__GetFmAleUtility().ComputeMeshMovement(delta_time)

            # Project the obtained MESH_VELOCITY and historical VELOCITY and PRESSURE values to the origin mesh
            buffer_size = self.main_model_part.GetBufferSize()
            domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

            if (domain_size == 2):
                self.__GetFmAleUtility().ProjectVirtualValues2D(self.main_model_part, buffer_size)
            else:
                self.__GetFmAleUtility().ProjectVirtualValues3D(self.main_model_part, buffer_size)

            # If FM-ALE is performed, calculate the EMBEDDED_VELOCITY
            self.__CalculateEmbeddedVelocity()

    def __UndoFMALEOperations(self):
        if self.__IsFmAleStep():
            # Undo the FM-ALE virtual mesh movement
            self.__GetFmAleUtility().UndoMeshMovement()

    def __CalculateEmbeddedVelocity(self):
        fm_ale_settings = self.settings["fm_ale_settings"]
        embedded_velocity_calculation = fm_ale_settings["embedded_velocity_calculation"].GetString()
        if embedded_velocity_calculation == "from_fluid_mesh_velocity":
            # Use the fluid MESH_VELOCITY as EMBEDDED_VELOCITY
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                KratosMultiphysics.MESH_VELOCITY,
                KratosMultiphysics.EMBEDDED_VELOCITY,
                self.GetComputingModelPart(),
                self.GetComputingModelPart(),
                0)
        elif embedded_velocity_calculation == "from_structure_velocity_rbf_interpolation":
            # If FM-ALE is performed (i.e. moving object), calculate the EMBEDDED_VELOCITY from the immersed object VELOCITY
            # This auxiliary function sets the EMBEDDED_VELOCITY in the non-historical nodal database of the intersected element
            # nodes from a Radial-Basis Function (RBF) interpolation of the structure VELOCITY values. These EMBEDDED_VELOCITY
            # nodal values would be later on interpolated in the cuts to impose the unfitted boundary condition within the element
            KratosCFD.FluidAuxiliaryUtilities.MapVelocityFromSkinToVolumeRBF(
                self.GetComputingModelPart(),
                self.__GetFmAleStructureModelPart(),
                self.__CalculateEmbeddedVelocityMapSearchRadius())
        else:
            err_msg = f"Provided 'embedded_velocity_calculation' vale {embedded_velocity_calculation} is not supported.\n"
            err_msg += "Available options are:\n"
            err_msg += "\t- 'from_fluid_mesh_velocity'\n"
            err_msg += "\t- 'from_structure_velocity_rbf_interpolation'\n"
            raise Exception(err_msg)

    def __CalculateEmbeddedVelocityMapSearchRadius(self):
        # If not computed yet, calculate the search radius
        if not hasattr(self, '__embedded_velocity_map_search_radius'):
            # First check if a user-defined specific value is provided
            # Note that the default value is 0 in order to do this
            fm_ale_settings = self.settings["fm_ale_settings"]
            search_radius = fm_ale_settings["rbf_interpolation_search_radius"].GetDouble()
            if search_radius > 0.0:
                self.__embedded_velocity_map_search_radius = search_radius
            else:
                # If not calculate it from the largest edge distance
                max_edge_factor = fm_ale_settings["rbf_interpolation_edge_factor"].GetDouble()
                max_edge_length = KratosCFD.FluidAuxiliaryUtilities.FindMaximumEdgeLength(self.GetComputingModelPart())
                self.__embedded_velocity_map_search_radius = max_edge_factor * max_edge_length

        return self.__embedded_velocity_map_search_radius
