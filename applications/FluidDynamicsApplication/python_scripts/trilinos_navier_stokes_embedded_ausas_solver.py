from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
import trilinos_navier_stokes_embedded_solver

def CreateSolver(model, custom_settings):
    return NavierStokesMPIEmbeddedAusasMonolithicSolver(model, custom_settings)

class NavierStokesMPIEmbeddedAusasMonolithicSolver(trilinos_navier_stokes_embedded_solver.NavierStokesMPIEmbeddedMonolithicSolver):

    def _ValidateSettings(self, settings):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "EmbeddedAusas",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_GID_file",
                "distance_file_name"  : "distance_file"
            },
            "maximum_iterations": 7,
            "dynamic_tau": 1.0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"                        : "MultiLevelSolver",
                "max_iteration"                      : 200,
                "tolerance"                          : 1e-6,
                "max_levels"                         : 3,
                "symmetric"                          : false,
                "reform_preconditioner_at_each_step" : true,
                "scaling"                            : true
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "periodic": "periodic",
            "move_mesh_flag": false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):

        super(NavierStokesMPIEmbeddedAusasMonolithicSolver,self).__init__(model,custom_settings)

        self.element_name = "EmbeddedAusasNavierStokes"
        self.condition_name = "EmbeddedAusasNavierStokesWallCondition"
        self.min_buffer_size = 3

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPIEmbeddedAusasMonolithicSolver","Construction of NavierStokesEmbeddedAusasSolver finished.")


    def Initialize(self):
        # Initialize the solver as in the base embedded solver
        super(NavierStokesMPIEmbeddedAusasMonolithicSolver, self).Initialize()

        # Set the find nodal neighbours process used in the embedded Ausas formulation condition
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.find_nodal_neighbours_process = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part,
                                                                                           number_of_avg_elems,
                                                                                           number_of_avg_nodes)
        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPIEmbeddedAusasMonolithicSolver","Monolithic embedded Ausas fluid solver initialization finished.")


    def InitializeSolutionStep(self):
        (self.bdf_process).Execute()
        (self.find_nodal_neighbours_process).Execute()
        (self.solver).InitializeSolutionStep()


    def Solve(self):
        (self.bdf_process).Execute()
        (self.find_nodal_neighbours_process).Execute()

        if self._TimeBufferIsInitialized():
            (self.solver).Solve()
