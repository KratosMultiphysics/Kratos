# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
import time as timer
from KratosMultiphysics.ShapeOptimizationApplication.mesh_controller_base import MeshController
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis

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
                "compute_reactions"         : false,
                "calculate_mesh_velocity"   : false
            },
            "processes" : {
                "boundary_conditions_process_list" : []
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

        self._mesh_moving_analysis = MeshMovingAnalysis(model, self.MeshSolverSettings)

    # --------------------------------------------------------------------------
    def Initialize(self):
        if self.has_automatic_boundary_process:
            KSO.GeometryUtilities(self.OptimizationModelPart).ExtractBoundaryNodes("auto_surface_nodes")

        self._mesh_moving_analysis.Initialize()

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

        if not self._mesh_moving_analysis.time < self._mesh_moving_analysis.end_time:
            self._mesh_moving_analysis.end_time += 1
        self._mesh_moving_analysis.RunSolutionLoop()

        KSO.MeshControllerUtilities(self.OptimizationModelPart).LogMeshChangeAccordingInputVariable(KM.MESH_DISPLACEMENT)

        self.OptimizationModelPart.ProcessInfo.SetValue(KM.STEP, step_before_update)
        self.OptimizationModelPart.ProcessInfo.SetValue(KM.TIME, time_before_update)
        self.OptimizationModelPart.ProcessInfo.SetValue(KM.DELTA_TIME, delta_time_before_update)

        KM.Logger.PrintInfo("ShapeOpt", "Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

    # --------------------------------------------------------------------------
    def Finalize(self):
        self._mesh_moving_analysis.Finalize()

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

# ==============================================================================
