from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Import base class file
import navier_stokes_embedded_solver

def CreateSolver(model, structure_model_part, custom_settings):
    return NavierStokesEmbeddedFMALEMonolithicSolver(model, custom_settings)

class NavierStokesEmbeddedFMALEMonolithicSolver(navier_stokes_embedded_solver.NavierStokesEmbeddedMonolithicSolver):

    def _ValidateSettings(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "embedded_solver_from_defaults",
            "model_part_name": "",
            "structure_model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "distance_reading_settings": {
                "import_mode": "from_mdpa",
                "distance_file_name": "no_distance_file"
            },
            "maximum_iterations": 7,
            "dynamic_tau": 1.0,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings": {
                "solver_type": "AMGCL"
            },
            "volume_model_part_name": "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step": true,
                "CFL_number": 1,
                "minimum_delta_time": 1e-2,
                "maximum_delta_time": 1.0
            },
            "move_mesh_flag": false,
            "fm_ale_settings": {
                "fm_ale_step_frequency": 1,
                "search_radius" : 1.0
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        if settings["structure_model_part_name"].GetString() == "":
            raise Exception('Please provide the name of the fixed model part as the "structure_model_part_name" (string) parameter!')

        return settings

    def __init__(self, model, custom_settings):
        super(NavierStokesEmbeddedFMALEMonolithicSolver,self).__init__(model,custom_settings)

        self.element_name = "EmbeddedNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.min_buffer_size = 3

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        # Retrieve the structural model part using json input
        structure_model_part_name = self.settings["model_part_name"].GetString()
        if self.model.HasModelPart(structure_model_part_name):
            self.structure_model_part = self.model.GetModelPart(structure_model_part_name)
        else:
            raise Exception("Structural model part {0} not found in model".format(structure_model_part_name))

        ## Set the FM-ALE framework
        self.fm_ale_step_frequency = self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt()
        if (self.fm_ale_step_frequency != 0):
            self.fm_ale_step = 1
            self.virtual_model_part = KratosMultiphysics.ModelPart("VirtualModelPart")
            self.mesh_moving_util = KratosMeshMoving.ExplicitMeshMovingUtilities(
                self.virtual_model_part, self.structure_model_part, self.settings["fm_ale_settings"]["search_radius"].GetDouble())

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedFMALEMonolithicSolver", "Construction of NavierStokesEmbeddedFMALEMonolithicSolver finished.")


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
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_PRESSURE)          # Post-process variable (stores the fluid nodes pressure and is set to 0 in the structure ones)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_VELOCITY)          # Post-process variable (stores the fluid nodes velocity and is set to 0 in the structure ones)

        if (self.fm_ale_step_frequency != 0):
            self.virtual_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
            self.virtual_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedFMALEMonolithicSolver", "Fluid solver variables added correctly.")


    def Initialize(self):
        super(NavierStokesEmbeddedFMALEMonolithicSolver, self).Initialize()

        ## Set the virtual model part geometry
        if (self.fm_ale_step_frequency != 0):
            self.mesh_moving_util.FillVirtualModelPart(self.main_model_part)


    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Compute the BDF coefficients
            (self.bdf_process).Execute()

            # Perform the FM-ALE operations
            if (self._is_fm_ale_step()):
                self._do_fm_ale_operations()
            else:
                self.fm_ale_step += 1

            # Fluid solver step initialization
            (self.solver).InitializeSolutionStep()


    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Fluid solver finalize solution step
            (self.solver).FinalizeSolutionStep()

            # Undo the virtual model part movement
            if (self._is_fm_ale_step()):
                self.mesh_moving_util.UndoMeshMovement()


    def _is_fm_ale_step(self):
        if (self.fm_ale_step_frequency != 0):
            if (self.fm_ale_step == self.fm_ale_step_frequency):
                return True
            else:
                return False

    def _do_fm_ale_operations(self):
        # Fill the virtual model part variable values: VELOCITY (n,nn), PRESSURE (n,nn)
        for i_step in range(self.main_model_part.GetBufferSize()):
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.PRESSURE, self.main_model_part, self.virtual_model_part, i_step)
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.VELOCITY, self.main_model_part, self.virtual_model_part, i_step)

        # Solve the mesh problem
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.mesh_moving_util.ComputeExplicitMeshMovement(delta_time)

        # Project the obtained MESH_VELOCITY and historical VELOCITY and PRESSURE values to the origin mesh
        buffer_size = self.main_model_part.GetBufferSize()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if (domain_size == 2):
            self.mesh_moving_util.ProjectVirtualValues2D(self.main_model_part, buffer_size)
        else:
            self.mesh_moving_util.ProjectVirtualValues3D(self.main_model_part, buffer_size)

