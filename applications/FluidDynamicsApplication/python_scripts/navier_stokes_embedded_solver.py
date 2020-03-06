from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

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

        if formulation_settings.Has("element_type"):
            element_type = formulation_settings["element_type"].GetString()
            if element_type == "embedded_navier_stokes":
                self._SetUpClassicEmbeddedNavierStokes(formulation_settings)
            elif element_type == "embedded_symbolic_navier_stokes":
                self._SetUpEmbeddedSymbolicNavierStokes(formulation_settings)
            elif element_type == "embedded_ausas_navier_stokes":
                self._SetUpClassicEmbeddedAusasNavierStokes(formulation_settings)
            elif element_type == "embedded_symbolic_navier_stokes_discontinuous":
                self._SetUpEmbeddedSymbolicNavierStokesDiscontinuous(formulation_settings)
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
        self.element_has_nodal_properties = True

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosCFD.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        if formulation_settings["is_slip"].GetBool():
            self.process_info_data[KratosCFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()

    def _SetUpEmbeddedSymbolicNavierStokes(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_symbolic_navier_stokes",
            "is_slip": false,
            "slip_length": 1.0e8,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0,
            "level_set_type": "continuous"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedSymbolicNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        self.element_has_nodal_properties = False

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosCFD.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
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
        self.element_has_nodal_properties = True

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosCFD.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()

    def _SetUpEmbeddedSymbolicNavierStokesDiscontinuous(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_symbolic_navier_stokes_discontinuous",
            "is_slip": true,
            "slip_length": 1.0e8,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0,
            "level_set_type": "discontinuous"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedSymbolicNavierStokesDiscontinuous"
        self.condition_name = "NavierStokesWallCondition"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        self.element_has_nodal_properties = False

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosCFD.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
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
                "gradient_penalty_coefficient": 0.0,
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
    def GetDefaultSettings(cls):
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
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": false,
            "no_skin_parts":[""],
            "time_stepping"                : {
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
                "fm_ale_solver_settings": {
                }
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesEmbeddedMonolithicSolver, cls).GetDefaultSettings())
        return default_settings

    def ValidateSettings(self):
        """Overriding python_solver ValidateSettings to validate the fm_ale_settings
        """
        super(NavierStokesEmbeddedMonolithicSolver, self).ValidateSettings()
        self.settings["fm_ale_settings"].ValidateAndAssignDefaults(self.GetDefaultSettings()["fm_ale_settings"])
        if self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() > 0:
            mesh_movement = self.settings["fm_ale_settings"]["mesh_movement"].GetString()
            self.settings["fm_ale_settings"]["fm_ale_solver_settings"].ValidateAndAssignDefaults(self._get_fm_ale_solver_default_settings(mesh_movement))

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesEmbeddedMonolithicSolver,self).__init__(model,custom_settings)

        self.min_buffer_size = 3
        self.embedded_formulation = EmbeddedFormulation(self.settings["formulation"])
        self.element_name = self.embedded_formulation.element_name
        self.condition_name = self.embedded_formulation.condition_name

        ## Set the formulation level set type
        self.level_set_type = self.embedded_formulation.level_set_type

        ## Set the nodal properties flag
        self.element_has_nodal_properties = self.embedded_formulation.element_has_nodal_properties

        ## Construct the linear solver
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        # If the FM-ALE is required, do a first call to the _get_fm_ale_virtual_model_part
        # Note that this will create the virtual model part in the model
        self.__fm_ale_is_active = self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() > 0
        if self.__fm_ale_is_active:
            self._get_fm_ale_virtual_model_part()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedMonolithicSolver", "Construction of NavierStokesEmbeddedMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
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

        if self.__fm_ale_is_active:
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedMonolithicSolver", "Fluid solver variables added correctly.")

    def PrepareModelPart(self):
        # Call the base solver PrepareModelPart()
        super(NavierStokesEmbeddedMonolithicSolver, self).PrepareModelPart()

        # Set the extra requirements of the embedded formulation
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Sets the embedded formulation configuration
            self._set_embedded_formulation()
            ## Setting the nodal distance
            self._set_distance_function()

    def Initialize(self):
        computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Set the time discretization utility to compute the BDF coefficients
        time_order = self.settings["time_order"].GetInt()
        if time_order == 2:
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
        else:
            raise Exception("Only \"time_order\" equal to 2 is supported. Provided \"time_order\": " + str(time_order))

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],   # Domain size (2,3)
                                                                                        self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1) # DOFs (3,4)

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.

        # Set the distance modification process
        self._GetDistanceModificationProcess().ExecuteInitialize()

        # For the primitive Ausas formulation, set the find nodal neighbours process
        # Recall that the Ausas condition requires the nodal neighbours.
        if (self.settings["formulation"]["element_type"].GetString() == "embedded_ausas_navier_stokes"):
            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            self.find_nodal_neighbours_process = KratosMultiphysics.FindNodalNeighboursProcess(self.GetComputingModelPart())

        # If required, intialize the FM-ALE utility
        if self.__fm_ale_is_active:
            self.fm_ale_step = 1
            # Fill the virtual model part geometry. Note that the mesh moving util is created in this first call
            self._get_mesh_moving_util().Initialize(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedMonolithicSolver", "Solver initialization finished.")

    def AdvanceInTime(self, current_time):
        # Call base solver AdvanceInTime to clone the time step and get the new time
        new_time = super(NavierStokesEmbeddedMonolithicSolver, self).AdvanceInTime(current_time)

        # Save the current step and time in the virtual model part process info
        if self.__fm_ale_is_active:
            self._get_fm_ale_virtual_model_part().ProcessInfo[KratosMultiphysics.STEP] += 1
            self._get_fm_ale_virtual_model_part().ProcessInfo[KratosMultiphysics.TIME] = new_time

        return new_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Compute the BDF coefficients
            (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

            # If required, compute the nodal neighbours
            if (self.settings["formulation"]["element_type"].GetString() == "embedded_ausas_navier_stokes"):
                (self.find_nodal_neighbours_process).Execute()

            # Set the virtual mesh values from the background mesh
            self._set_virtual_mesh_values()

        # Call the base solver InitializeSolutionStep()
        super(NavierStokesEmbeddedMonolithicSolver, self).InitializeSolutionStep()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Correct the distance field
            # Note that this is intentionally placed in here (and not in the InitializeSolutionStep() of the solver
            # It has to be done before each call to the Solve() in case an outer non-linear iteration is performed (FSI)
            self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

            # Perform the FM-ALE operations
            # Note that this also sets the EMBEDDED_VELOCITY from the MESH_VELOCITY
            self._do_fm_ale_operations()

            # Call the base SolveSolutionStep to solve the embedded CFD problem
            is_converged = super(NavierStokesEmbeddedMonolithicSolver,self).SolveSolutionStep()

            # Undo the FM-ALE virtual mesh movement
            self.__UndoFMALEOperations()

            # Restore the fluid node fixity to its original status
            # Note that this is intentionally placed in here (and not in the FinalizeSolutionStep() of the solver
            # It has to be done after each call to the Solve() and the FM-ALE in case an outer non-linear iteration is performed (FSI)
            self._GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        # Call the base solver FinalizeSolutionStep()
        super(NavierStokesEmbeddedMonolithicSolver, self).FinalizeSolutionStep()

        # Do the FM-ALE end of step operations
        if self._TimeBufferIsInitialized():
            self.__UpdateFMALEStepCounter()

    def _SetPhysicalProperties(self):
        # Call the base solver _SetPhysicalProperties()
        materials_imported = super(NavierStokesEmbeddedMonolithicSolver, self)._SetPhysicalProperties()

        # Check if the SOUND_VELOCITY has been defined by the user
        user_defined_sound_velocity = False
        for elem in self.main_model_part.Elements:
            if elem.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
                user_defined_sound_velocity = True
                sound_velocity = elem.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
            break

        # Set the SOUND_VELOCITY value (wave velocity)
        # TODO: Save the SOUND_VELOCITY in the element Properties
        if user_defined_sound_velocity:
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = sound_velocity
        else:
            # If the wave velocity is not defined take a large enough value to consider the fluid as incompressible
            default_sound_velocity = 1e+12
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = default_sound_velocity

        return materials_imported

    def _SetNodalProperties(self):
        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")
        # Transfer the obtained properties to the nodes
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity, self.main_model_part.Nodes)

    def _set_constitutive_law(self):
        ## Construct the constitutive law needed for the embedded element
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian3DLaw()
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian2DLaw()

    def _set_embedded_formulation(self):
        # Set the SLIP elemental flag
        if (self.settings["formulation"]["is_slip"].GetBool()):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.GetComputingModelPart().Elements)
        else:
            # Set the SLIP elemental flag to false in the entire domain
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, False, self.GetComputingModelPart().Elements)

        # Save the formulation settings in the ProcessInfo
        self.embedded_formulation.SetProcessInfo(self.main_model_part)

    def _set_distance_function(self):
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

    def _GetDistanceModificationProcess(self):
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

    def _get_fm_ale_structure_model_part(self):
        structure_model_part_name = self.settings["fm_ale_settings"]["fm_ale_solver_settings"]["structure_model_part_name"].GetString()
        if self.model.HasModelPart(structure_model_part_name):
            return self.model.GetModelPart(structure_model_part_name)
        else:
            raise Exception("Structure model part {0} not found in model.\n It is expected to be added in your custom analysis stage file.".format(structure_model_part_name))

    def _get_fm_ale_virtual_model_part(self):
        if not hasattr(self, '_virtual_model_part'):
            self._virtual_model_part = self._create_fm_ale_virtual_model_part()
        return self._virtual_model_part

    def _create_fm_ale_virtual_model_part(self):
        virtual_model_part_name = self.settings["fm_ale_settings"]["fm_ale_solver_settings"]["virtual_model_part_name"].GetString()
        virtual_model_part = self.model.CreateModelPart(virtual_model_part_name)
        return virtual_model_part

    def _get_mesh_moving_util(self):
        if not hasattr (self, '_mesh_moving_util'):
            self._mesh_moving_util = self._create_mesh_moving_util()
        return self._mesh_moving_util

    def _create_mesh_moving_util(self):
        if have_mesh_moving:
            mesh_movement = self.settings["fm_ale_settings"]["mesh_movement"].GetString()
            if (mesh_movement == "implicit"):
                mesh_moving_util = KratosMeshMoving.FixedMeshALEUtilities(
                    self.model,
                    self.settings["fm_ale_settings"]["fm_ale_solver_settings"])
            elif (mesh_movement == "explicit"):
                mesh_moving_util = KratosMeshMoving.ExplicitFixedMeshALEUtilities(
                    self.model,
                    self.settings["fm_ale_settings"]["fm_ale_solver_settings"])
            else:
                raise Exception("FM-ALE mesh_movement set to \'" + mesh_movement + "\'. Available options are \'implicit\' and \'explicit\'.")

            return mesh_moving_util
        else:
            raise Exception("MeshMovingApplication is required to construct the FM-ALE utility (ExplicitFixedMeshALEUtilities)")

    def _is_fm_ale_step(self):
        if self.__fm_ale_is_active:
            if (self.fm_ale_step == self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt()):
                return True
            else:
                return False
        else:
            return False

    def __UpdateFMALEStepCounter(self):
        if self.__fm_ale_is_active:
            if (self._is_fm_ale_step()):
                # Reset the FM-ALE steps counter
                self.fm_ale_step = 1
            else:
                # Update the FM-ALE steps counter
                self.fm_ale_step += 1

    def _set_virtual_mesh_values(self):
        if self._is_fm_ale_step():
            # Fill the virtual model part variable values: VELOCITY (n,nn), PRESSURE (n,nn)
            self._get_mesh_moving_util().SetVirtualMeshValuesFromOriginMesh()

    def _do_fm_ale_operations(self):
        if self._is_fm_ale_step():
            # Solve the mesh problem
            delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            self._get_mesh_moving_util().ComputeMeshMovement(delta_time)

            # Project the obtained MESH_VELOCITY and historical VELOCITY and PRESSURE values to the origin mesh
            buffer_size = self.main_model_part.GetBufferSize()
            domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

            if (domain_size == 2):
                self._get_mesh_moving_util().ProjectVirtualValues2D(self.main_model_part, buffer_size)
            else:
                self._get_mesh_moving_util().ProjectVirtualValues3D(self.main_model_part, buffer_size)

            # If FM-ALE is performed, use the MESH_VELOCITY as EMBEDDED_VELOCITY
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                KratosMultiphysics.MESH_VELOCITY,
                KratosMultiphysics.EMBEDDED_VELOCITY,
                self.GetComputingModelPart(),
                self.GetComputingModelPart(),
                0)

    def __UndoFMALEOperations(self):
        if self._is_fm_ale_step():
            # Undo the FM-ALE virtual mesh movement
            self._get_mesh_moving_util().UndoMeshMovement()
