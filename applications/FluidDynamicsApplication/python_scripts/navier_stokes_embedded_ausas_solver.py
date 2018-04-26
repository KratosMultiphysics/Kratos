from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
import navier_stokes_embedded_solver

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesEmbeddedAusasMonolithicSolver(main_model_part, custom_settings)

class NavierStokesEmbeddedAusasMonolithicSolver(navier_stokes_embedded_solver.NavierStokesEmbeddedMonolithicSolver):

    def __init__(self, main_model_part, custom_settings):

        self.element_name = "EmbeddedAusasNavierStokes"
        self.condition_name = "EmbeddedAusasNavierStokesWallCondition"
        self.min_buffer_size = 3

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "EmbeddedAusas",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
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
                "solver_type"         : "AMGCL_NS_Solver"
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
            "reorder": false,
            "fm_ale_settings": {
                "fm_ale_step_frequency": 0,
                "mesh_solver_settings":{
                    "problem_data":{
                        "parallel_type": "OpenMP"
                    },
                    "solver_settings":{
                        "solver_type": "mesh_solver_structural_similarity"
                    }
                }
            }
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        ## Set the FM-ALE framework
        self.fm_ale_step_frequency = self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt()
        if (self.fm_ale_step_frequency != 0):
            self.fm_ale_step = 1
            self.virtual_model_part = ModelPart("VirtualModelPart")
            import python_solvers_wrapper_mesh_motion
            self.mesh_solver = python_solvers_wrapper_mesh_motion.CreateSolver(self.virtual_model_part,self.settings["mesh_solver_settings"]["solver_settings"])

        if self._IsPrintingRank():
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
        (self.solver).InitializeSolutionStep()


    def Solve(self):
        (self.bdf_process).Execute()
        (self.find_nodal_neighbours_process).Execute()

        # Note that the first two time steps are dropped to fill the BDF buffer
        if (self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] >= 2):
            (self.solver).Solve()
