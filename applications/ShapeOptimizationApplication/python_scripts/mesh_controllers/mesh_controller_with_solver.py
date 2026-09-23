# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================


# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
import time as timer
from KratosMultiphysics.ShapeOptimizationApplication.mesh_controllers.mesh_controller_base import MeshController
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis

try:
    import KratosMultiphysics.MeshingApplication as KMA
    from KratosMultiphysics.MeshingApplication.mmg_process import MmgProcess as automatic_remeshing_process
    if not hasattr(KMA, "MmgProcess2D"):
        automatic_remeshing_process = None
        automatic_remeshing_error_msg = "MeshingApplication is not compiled with '-DINCLUDE_MMG=ON'"
except ImportError as err:
    automatic_remeshing_process = None
    automatic_remeshing_error_msg = str(err)

# # ==============================================================================
class MeshControllerWithSolver(MeshController) :
    # --------------------------------------------------------------------------
    def __init__(self, MeshSolverSettings, model):
        default_settings = KM.Parameters("""
        {
            "apply_mesh_solver" : true,
            "solver_settings" : {
                "domain_size"     : 3,
                "echo_level"      : 0,
                "solver_type"     : "structural_similarity",
                "model_part_name" : "NONE",
                "model_import_settings"              : {
                    "input_type"     : "use_input_model_part"
                },
                "time_stepping" : {
                    "time_step"       : 1.0
                },
                "linear_solver_settings" : {
                    "solver_type" : "amgcl",
                    "smoother_type":"ilu0",
                    "krylov_type": "gmres",
                    "coarsening_type": "aggregation",
                    "max_iteration": 200,
                    "verbosity" : 0,
                    "tolerance": 1e-7
                },
                "compute_reactions"                : false,
                "calculate_mesh_velocity"          : false
            },
            "processes" : {
                "boundary_conditions_process_list" : []
            },
            "use_automatic_remeshing"     : false,
            "skip_first_remeshing_step"   : false,
            "automatic_remeshing_settings": {
                "strategy"        : "optimization",
                "step_frequency"  : 1,
                "automatic_remesh": true,
                "automatic_remesh_parameters": {
                    "automatic_remesh_type": "Ratio",
                    "min_size_ratio": 1.0,
                    "max_size_ratio": 5.0,
                    "refer_type"    : "Mean"
                },
                "echo_level": 0,
                "force_min" : true,
                "force_max" : true
            }
        }""")
        self.MeshSolverSettings = MeshSolverSettings
        self.MeshSolverSettings.ValidateAndAssignDefaults(default_settings)
        self.MeshSolverSettings["solver_settings"].ValidateAndAssignDefaults(default_settings["solver_settings"])
        self.MeshSolverSettings["processes"].ValidateAndAssignDefaults(default_settings["processes"])

        if not self.MeshSolverSettings["solver_settings"].Has("linear_solver_settings"):
            MeshSolverSettings.AddValue("linear_solver_settings", default_settings["solver_settings"]["linear_solver_settings"])
            KM.Logger.PrintInfo("ShapeOpt", "::[MeshControllerWithSolver]:: using default linear solver for mesh motion.")

        if not MeshSolverSettings.Has("problem_data"):
            self.__AddDefaultProblemData(self.MeshSolverSettings)
        else:
            KM.Logger.PrintInfo("ShapeOpt", "::[MeshControllerWithSolver]:: using custom problem data for mesh motion.")

        self.OptimizationModelPart = model[self.MeshSolverSettings["solver_settings"]["model_part_name"].GetString()]

        if self.MeshSolverSettings["processes"]["boundary_conditions_process_list"].size() == 0:
            self.__FixWholeSurface(self.OptimizationModelPart, self.MeshSolverSettings)
            self.has_automatic_boundary_process = True
        else:
            self.has_automatic_boundary_process = False

        self.is_remeshing_used = self.MeshSolverSettings["use_automatic_remeshing"].GetBool()
        self.skip_first_remshing_step = self.MeshSolverSettings["skip_first_remeshing_step"].GetBool()
        if (self.is_remeshing_used):
            automatic_remeshing_process_settings = self.MeshSolverSettings["automatic_remeshing_settings"]
            if (automatic_remeshing_process is None):
                raise RuntimeError("Automatic remeshing requires to import MeshingApplication. Importing failed with following error msg.\n\t" + automatic_remeshing_error_msg)

            self.__CheckAndSetAutomaticMeshRefinementSettings(automatic_remeshing_process_settings)
            self.remeshing_process = automatic_remeshing_process(model, automatic_remeshing_process_settings)

            # remeshing requires to reinitialize the model_part of the mesh solver
            self.MeshSolverSettings["solver_settings"].AddBool("reinitialize_model_part_each_step", True)

            KM.Logger.PrintInfo("ShapeOpt", "Initialized automatic automatic remeshing process")

        self._mesh_moving_analysis = MeshMovingAnalysis(model, self.MeshSolverSettings)

    # --------------------------------------------------------------------------
    def Initialize(self):
        if self.has_automatic_boundary_process:
            KSO.GeometryUtilities(self.OptimizationModelPart).ExtractBoundaryNodes("auto_surface_nodes")

        self._mesh_moving_analysis.Initialize()
        if self.is_remeshing_used:
            self.remeshing_process.ExecuteInitialize()

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, variable):
        KM.Logger.Print("")
        KM.Logger.PrintInfo("ShapeOpt", "Starting to update the mesh...")
        startTime = timer.time()

        time_before_update = self.OptimizationModelPart.ProcessInfo.GetValue(KM.TIME)
        step_before_update = self.OptimizationModelPart.ProcessInfo.GetValue(KM.STEP)
        delta_time_before_update = self.OptimizationModelPart.ProcessInfo.GetValue(KM.DELTA_TIME)

        # Reset step/time iterators such that they match the current iteration after calling RunSolutionLoop (which internally calls CloneTimeStep)
        self.OptimizationModelPart.ProcessInfo.SetValue(KM.STEP, step_before_update-1)
        self.OptimizationModelPart.ProcessInfo.SetValue(KM.TIME, time_before_update-1)
        self.OptimizationModelPart.ProcessInfo.SetValue(KM.DELTA_TIME, 0)

        KM.VariableUtils().CopyVectorVar(variable, KM.MESH_DISPLACEMENT, self.OptimizationModelPart.Nodes)

        if self.has_automatic_boundary_process and self.is_remeshing_used:
            self.OptimizationModelPart.GetSubModelPart("auto_surface_nodes").GetNodes().clear()
            KSO.GeometryUtilities(self.OptimizationModelPart).ExtractBoundaryNodes("auto_surface_nodes")

        if not self._mesh_moving_analysis.time < self._mesh_moving_analysis.end_time:
            self._mesh_moving_analysis.end_time += 1
        self._mesh_moving_analysis.RunSolutionLoop()

        if self.is_remeshing_used and not self.skip_first_remshing_step:
            self.OptimizationModelPart.Set(KM.MODIFIED, False)
            self.remeshing_process.ExecuteInitializeSolutionStep()
            self.remeshing_process.ExecuteFinalizeSolutionStep()

        KSO.MeshControllerUtilities(self.OptimizationModelPart).LogMeshChangeAccordingInputVariable(KM.MESH_DISPLACEMENT)

        self.OptimizationModelPart.ProcessInfo.SetValue(KM.STEP, step_before_update)
        self.OptimizationModelPart.ProcessInfo.SetValue(KM.TIME, time_before_update)
        self.OptimizationModelPart.ProcessInfo.SetValue(KM.DELTA_TIME, delta_time_before_update)

        self.skip_first_remshing_step = False

        KM.Logger.PrintInfo("ShapeOpt", "Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

    # --------------------------------------------------------------------------
    def Finalize(self):
        self._mesh_moving_analysis.Finalize()

        if self.is_remeshing_used and not self.skip_first_remshing_step:
            self.remeshing_process.ExecuteFinalize()

    # --------------------------------------------------------------------------
    @staticmethod
    def __AddDefaultProblemData(mesh_solver_settings):
        problem_data = KM.Parameters("""{
            "echo_level"    : 0,
            "start_time"    : 0.0,
            "end_time"      : 1.0,
            "parallel_type" : "OpenMP"
        }""")

        mesh_solver_settings.AddValue("problem_data", problem_data)

    # --------------------------------------------------------------------------
    @staticmethod
    def __FixWholeSurface(optimization_model_part, mesh_solver_settings):
        optimization_model_part.CreateSubModelPart("auto_surface_nodes")

        auto_process_settings = KM.Parameters(
            """
            {
                "python_module" : "fix_vector_variable_process",
                "kratos_module" : "KratosMultiphysics",
                "help"          : "This process fixes the selected components of a given vector variable without modifying the value of the variable.",
                "process_name"  : "FixVectorVariableProcess",
                "Parameters"    : {
                    "model_part_name"      : \""""+str(optimization_model_part.Name)+""".auto_surface_nodes\",
                    "variable_name"        : "MESH_DISPLACEMENT",
                    "constrained"          : [true,true,true]
                }
            }
            """)

        KM.Logger.PrintInfo("ShapeOpt", "Add automatic process to fix the whole surface to mesh motion solver:")
        mesh_solver_settings["processes"]["boundary_conditions_process_list"].Append(auto_process_settings)

    # --------------------------------------------------------------------------
    def __CheckAndSetAutomaticMeshRefinementSettings(self, parameters):
        if (parameters.Has("interpolate_nodal_values")):
            if (parameters["interpolate_nodal_values"].GetBool()):
                KM.Logger.PrintWarning("ShapeOpt", "Historical value interpolation is not allowed in automatic remeshing. Turning it off.")
            parameters["interpolate_nodal_values"].SetBool(False)
        else:
            parameters.AddBool("interpolate_nodal_values", False)

        if (parameters.Has("interpolate_non_historical")):
            if (parameters["interpolate_non_historical"].GetBool()):
                KM.Logger.PrintWarning("ShapeOpt", "Non-historical value interpolation is not allowed in automatic remeshing. Turning it off.")
            parameters["interpolate_non_historical"].SetBool(False)
        else:
            parameters.AddBool("interpolate_non_historical", False)

        if (parameters.Has("extrapolate_contour_values")):
            if (parameters["extrapolate_contour_values"].GetBool()):
                KM.Logger.PrintWarning("ShapeOpt", "Value extrapolation is not allowed in automatic remeshing. Turning it off.")
            parameters["extrapolate_contour_values"].SetBool(False)
        else:
            parameters.AddBool("extrapolate_contour_values", False)

        if (parameters.Has("model_part_name")):
            if (parameters["model_part_name"].GetString() != self.OptimizationModelPart.Name):
                KM.Logger.PrintWarning("ShapeOpt", "Mismatching model part name provided for automatic remeshing [ " + parameters["model_part_name"].GetString() + " ]. Using the optimization model part [ " + self.OptimizationModelPart.Name + " ].")
            parameters["model_part_name"].SetString(self.OptimizationModelPart.Name)
        else:
            parameters.AddString("model_part_name", self.OptimizationModelPart.Name)

# ==============================================================================
