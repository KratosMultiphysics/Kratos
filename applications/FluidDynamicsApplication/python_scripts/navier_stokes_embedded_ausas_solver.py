from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
import navier_stokes_embedded_solver

def CreateSolver(model, custom_settings):
    return NavierStokesEmbeddedAusasMonolithicSolver(model, custom_settings)

class NavierStokesEmbeddedAusasMonolithicSolver(navier_stokes_embedded_solver.NavierStokesEmbeddedMonolithicSolver):

    def _ValidateSettings(self, settings):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "EmbeddedAusas",
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
            "dynamic_tau": 1.0,
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
            "periodic": "periodic",
            "move_mesh_flag": false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):

        super(NavierStokesEmbeddedAusasMonolithicSolver,self).__init__(model,custom_settings)

        self.element_name = "EmbeddedAusasNavierStokes"
        self.condition_name = "EmbeddedAusasNavierStokesWallCondition"
        self.min_buffer_size = 3

        KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedAusasMonolithicSolver", "Construction of NavierStokesEmbeddedAusasMonolithicSolver finished.")


    def Initialize(self):
        # Initialize the solver as in the base embedded solver
        super(NavierStokesEmbeddedAusasMonolithicSolver, self).Initialize()

        # Set the find nodal neighbours process used in the embedded Ausas formulation condition
        if (self.settings["solver_type"].GetString() == "EmbeddedAusas"):
            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            self.find_nodal_neighbours_process = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part,
                                                                                               number_of_avg_elems,
                                                                                               number_of_avg_nodes)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesEmbeddedAusasMonolithicSolver", "Solver initialization finished.")


    def InitializeSolutionStep(self):
        (self.bdf_process).Execute()
        (self.find_nodal_neighbours_process).Execute()
        if self._TimeBufferIsInitialized():
            (self.solver).InitializeSolutionStep()


    def Solve(self):
        (self.bdf_process).Execute()
        (self.find_nodal_neighbours_process).Execute()

        # Note that the first two time steps are dropped to fill the BDF buffer
        if self._TimeBufferIsInitialized():
            (self.solver).Solve()
