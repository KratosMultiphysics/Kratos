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
from KratosMultiphysics import *
from KratosMultiphysics.MeshMovingApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
import time as timer
from mesh_controller_base import MeshController
from mesh_moving_analysis import MeshMovingAnalysis

# # ==============================================================================
class MeshControllerWithSolver(MeshController) :
    # --------------------------------------------------------------------------
    def __init__(self, MeshSolverSettings, model):
        default_settings = Parameters("""
        {
            "apply_mesh_solver" : true,
            "solver_settings" : {
                "solver_type" : "mesh_solver_structural_similarity",
                "model_part_name"       : "",
                "model_import_settings"              : {
                    "input_type"     : "use_input_model_part",
                    "input_filename"  : "DUMMY_FILENAME"
                },
                "time_stepping" : {
                    "time_step"       : 1.0
                },
                "domain_size"     : 3,
                "linear_solver_settings" : {
                    "solver_type" : "AMGCL",
                    "smoother_type":"ilu0",
                    "krylov_type": "gmres",
                    "coarsening_type": "aggregation",
                    "max_iteration": 200,
                    "verbosity" : 0,
                    "tolerance": 1e-7
                },
                "compute_reactions"         : false,
                "calculate_mesh_velocities" : false
            },
            "boundary_conditions_process_list" : []
        }""")
        self.MeshSolverSettings = MeshSolverSettings
        self.MeshSolverSettings.ValidateAndAssignDefaults(default_settings)

        # add default problem data parameters for mesh motion analysis
        if not MeshSolverSettings.Has("problem_data"):
            problem_data = Parameters("""{
                "echo_level" : 0,
                "start_time" : 0.0,
                "end_time" : 1.0,
                "parallel_type" : "OpenMP"
            }""")

            self.MeshSolverSettings.AddValue("problem_data", problem_data)
        else:
            print("::[MeshControllerWithSolver]::WARNING: using custom problem data for mesh motion.")

        self.OptimizationModelPart = model[self.MeshSolverSettings["solver_settings"]["model_part_name"].GetString()]

        if self.MeshSolverSettings["boundary_conditions_process_list"].size() == 0:
            self.OptimizationModelPart.CreateSubModelPart("auto_surface_nodes")

            auto_process_settings = Parameters(
                """
                {
                    "python_module" : "fix_vector_variable_process",
                    "kratos_module" : "KratosMultiphysics.ShapeOptimizationApplication",
                    "help"          : "This process fixes the selected components of a given vector variable without modifying the value of the variable.",
                    "process_name"  : "FixVectorVariableProcess",
                    "Parameters"    : {
                        "model_part_name"      : \""""+str(self.OptimizationModelPart.Name)+""".auto_surface_nodes",
                        "variable_name"        : "MESH_DISPLACEMENT",
                        "constrained"          : [true,true,true]
                    }
                }
                """)

            print("Add automatic process to fix the whole surface to mesh motion solver:")
            print(auto_process_settings)
            self.MeshSolverSettings["boundary_conditions_process_list"].Append(auto_process_settings)
            self.has_automatic_boundary_process = True
        else:
            self.has_automatic_boundary_process = False

        self._mesh_moving_analysis = MeshMovingAnalysis(model, self.MeshSolverSettings)

    # --------------------------------------------------------------------------
    def Initialize(self):
        if self.has_automatic_boundary_process:
            GeometryUtilities(self.OptimizationModelPart).ExtractBoundaryNodes("auto_surface_nodes")

        self._mesh_moving_analysis.Initialize()

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, variable):
        print("\n> Starting to update the mesh...")
        startTime = timer.time()

        VariableUtils().CopyVectorVar(variable, MESH_DISPLACEMENT, self.OptimizationModelPart.Nodes)

        time_before_mesh_update = self.OptimizationModelPart.ProcessInfo.GetValue(TIME)

        if not self._mesh_moving_analysis.time < self._mesh_moving_analysis.end_time:
            self._mesh_moving_analysis.end_time += 1
        self._mesh_moving_analysis.RunSolutionLoop()

        self.OptimizationModelPart.ProcessInfo.SetValue(TIME, time_before_mesh_update)

        MeshControllerUtilities(self.OptimizationModelPart).LogMeshChangeAccordingInputVariable(MESH_DISPLACEMENT)

        print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

    # --------------------------------------------------------------------------
    def Finalize(self):
        self._mesh_moving_analysis.Finalize()

# ==============================================================================