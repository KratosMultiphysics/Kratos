from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
try:
    import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
    have_mesh_moving = True
except ImportError:
    have_mesh_moving = False


# Import base class file
from fluid_solver import FluidSolver

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
            "dynamic_tau": 1.0
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedNavierStokes"
        self.condition_name = "NavierStokesWallCondition"

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
            "dynamic_tau": 1.0
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedSymbolicNavierStokes"
        self.condition_name = "NavierStokesWallCondition"

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosCFD.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        if formulation_settings["is_slip"].GetBool():
            self.process_info_data[KratosCFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()

    def _SetUpClassicEmbeddedAusasNavierStokes(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_ausas_navier_stokes",
            "is_slip": true,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedAusasNavierStokes"
        self.condition_name = "EmbeddedAusasNavierStokesWallCondition"

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosCFD.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()

    def _SetUpEmbeddedSymbolicNavierStokesDiscontinuous(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_symbolic_navier_stokes_discontinuous",
            "is_slip": true,
            "slip_length": 1.0e8,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedSymbolicNavierStokesDiscontinuous"
        self.condition_name = "NavierStokesWallCondition"

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosCFD.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        if formulation_settings["is_slip"].GetBool():
            self.process_info_data[KratosCFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()
        else:
            self.process_info_data[KratosCFD.SLIP_LENGTH] = 0.0


def CreateSolver(model, custom_settings):
    return NavierStokesEmbeddedMonolithicSolver(model, custom_settings)

class NavierStokesEmbeddedMonolithicSolver(FluidSolver):

    def _ValidateSettings(self, settings):
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
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
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
                "solver_type"         : "AMGCL"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
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
                "structure_model_part_name": "",
                "search_radius" : 1.0
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):
        super(NavierStokesEmbeddedMonolithicSolver,self).__init__(model,custom_settings)

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        self.min_buffer_size = 3
        self.embedded_formulation = EmbeddedFormulation(self.settings["formulation"])
        self.element_name = self.embedded_formulation.element_name
        self.condition_name = self.embedded_formulation.condition_name

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedMonolithicSolver", "Construction of NavierStokesEmbeddedMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)              # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_PRESSURE)          # Post-process variable (stores the fluid nodes pressure and is set to 0 in the structure ones)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_VELOCITY)          # Post-process variable (stores the fluid nodes velocity and is set to 0 in the structure ones)

        if (self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() != 0):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
            self._get_fm_ale_virtual_model_part().AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self._get_fm_ale_virtual_model_part().AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
            self._get_fm_ale_virtual_model_part().AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
            self._get_fm_ale_virtual_model_part().AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedMonolithicSolver", "Fluid solver variables added correctly.")

    def ImportModelPart(self):
        super(NavierStokesEmbeddedMonolithicSolver, self).ImportModelPart()

    def PrepareModelPart(self):
        super(NavierStokesEmbeddedMonolithicSolver, self).PrepareModelPart()
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Sets DENSITY, DYNAMIC_VISCOSITY and SOUND_VELOCITY
            self._set_physical_properties()
            ## Sets the constitutive law
            self._set_constitutive_law()
            ## Sets the embedded formulation configuration
            self._set_embedded_formulation()
            ## Setting the nodal distance
            self._set_distance_function()

    def Initialize(self):
        computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(computing_model_part,
                                                                            self.settings["time_order"].GetInt())

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

        # For the primitive Ausas formulation, set the find nodal neighbours process
        # Recall that the Ausas condition requires the nodal neighbours.
        if (self.settings["formulation"]["element_type"].GetString() == "embedded_ausas_navier_stokes"):
            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            self.find_nodal_neighbours_process = KratosMultiphysics.FindNodalNeighboursProcess(self.GetComputingModelPart(),
                                                                                               number_of_avg_elems,
                                                                                               number_of_avg_nodes)

        # If required, intialize the FM-ALE utility
        if (self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() != 0):
            self.fm_ale_step = 1
            # Fill the virtual model part geometry. Note that the mesh moving util is created in this first call
            self._get_mesh_moving_util().FillVirtualModelPart(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedMonolithicSolver", "Solver initialization finished.")

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Compute the BDF coefficients
            (self.bdf_process).Execute()

            # If required, compute the nodal neighbours
            if (self.settings["formulation"]["element_type"].GetString() == "embedded_ausas_navier_stokes"):
                (self.find_nodal_neighbours_process).Execute()

            # Perform the FM-ALE operations
            self._do_fm_ale_operations()

            # Fluid solver step initialization
            (self.solver).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Fluid solver finalize solution step
            (self.solver).FinalizeSolutionStep()

            # Do the FM-ALE end of step operations
            self._finalize_fm_ale_step()

    def _set_physical_properties(self):
        ## Set the SOUND_VELOCITY value (wave velocity)
        if self.main_model_part.Properties[1].Has(KratosMultiphysics.SOUND_VELOCITY):
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = self.main_model_part.Properties[1][KratosMultiphysics.SOUND_VELOCITY]
        else:
            # If the wave velocity is not defined take a large enough value to consider the fluid as incompressible
            default_sound_velocity = 1e+12
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = default_sound_velocity

        # Transfer density and (dynamic) viscostity to the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            break

        # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity, self.main_model_part.Nodes)

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
            import read_distance_from_file
            DistanceUtility = read_distance_from_file.DistanceImportUtility(self.main_model_part, self.settings["distance_reading_settings"])
            DistanceUtility.ImportDistance()
        elif (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            KratosMultiphysics.Logger.PrintInfo("Navier Stokes Embedded Solver","Distance function taken from the .mdpa input file.")
            # Recall to swap the distance sign (GiD considers d<0 in the fluid region)
            for node in self.main_model_part.Nodes:
                distance_value = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -distance_value)

    def _get_fm_ale_structure_model_part(self):
        structure_model_part_name = self.settings["fm_ale_settings"]["structure_model_part_name"].GetString()
        if self.model.HasModelPart(structure_model_part_name):
            return self.model.GetModelPart(structure_model_part_name)
        else:
            raise Exception("Structure model part {0} not found in model.\n It is expected to be added in your custom analysis stage file.".format(structure_model_part_name))

    def _get_fm_ale_virtual_model_part(self):
        if not hasattr(self, '_virtual_model_part'):
            self._virtual_model_part = self._create_fm_ale_virtual_model_part()
        return self._virtual_model_part

    def _create_fm_ale_virtual_model_part(self):
        virtual_model_part = self.model.CreateModelPart("VirtualModelPart")
        return virtual_model_part

    def _get_mesh_moving_util(self):
        if not hasattr (self, '_mesh_moving_util'):
            self._mesh_moving_util = self._create_mesh_moving_util()
        return self._mesh_moving_util

    def _create_mesh_moving_util(self):
        if have_mesh_moving:
            mesh_moving_util = KratosMeshMoving.ExplicitMeshMovingUtilities(
                self._get_fm_ale_virtual_model_part(),
                self._get_fm_ale_structure_model_part(),
                self.settings["fm_ale_settings"]["search_radius"].GetDouble())
            return mesh_moving_util
        else:
            raise Exception("MeshMovingApplication is required to construct the FM-ALE utility (ExplicitMeshMovingUtilities)")

    def _is_fm_ale_step(self):
        if (self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() != 0):
            if (self.fm_ale_step == self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt()):
                return True
            else:
                return False

    def _finalize_fm_ale_step(self):
        if (self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() != 0):
            if (self._is_fm_ale_step()):
                # Undo virtual mesh movement
                self._get_mesh_moving_util().UndoMeshMovement()
                # Reset the FM-ALE steps counter
                self.fm_ale_step = 1
            else:
                # Update the FM-ALE steps counter
                self.fm_ale_step += 1

    def _do_fm_ale_operations(self):
        if (self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() != 0):
            # Fill the virtual model part variable values: VELOCITY (n,nn), PRESSURE (n,nn)
            for i_step in range(self.main_model_part.GetBufferSize()):
                KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                    KratosMultiphysics.PRESSURE, self.main_model_part, self._get_fm_ale_virtual_model_part(), i_step)
                KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                    KratosMultiphysics.VELOCITY, self.main_model_part, self._get_fm_ale_virtual_model_part(), i_step)

            # Solve the mesh problem
            delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            self._get_mesh_moving_util().ComputeExplicitMeshMovement(delta_time)

            # Project the obtained MESH_VELOCITY and historical VELOCITY and PRESSURE values to the origin mesh
            buffer_size = self.main_model_part.GetBufferSize()
            domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

            if (domain_size == 2):
                self._get_mesh_moving_util().ProjectVirtualValues2D(self.main_model_part, buffer_size)
            else:
                self._get_mesh_moving_util().ProjectVirtualValues3D(self.main_model_part, buffer_size)
