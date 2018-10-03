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
    def __init__(self, MeshSolverSettings, OptimizationModelPart):
        default_settings = Parameters("""
        {
            "apply_mesh_solver" : true,
            "problem_data" : {
                "echo_level" : 0,
                "time_step" : 1.1,
                "start_time" : 0.0,
                "end_time" : 1.0,
                "parallel_type" : "OpenMP"
            },
            "solver_settings" : {
                "solver_type" : "mesh_solver_structural_similarity",
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
            }
        }""")
        self.MeshSolverSettings = MeshSolverSettings
        self.MeshSolverSettings.ValidateAndAssignDefaults(default_settings)

        self.MeshSolverSettings["problem_data"].AddEmptyValue("domain_size")
        self.MeshSolverSettings["problem_data"]["domain_size"].SetInt(OptimizationModelPart.ProcessInfo[DOMAIN_SIZE])

        self.MeshSolverSettings["problem_data"].AddEmptyValue("model_part_name")
        self.MeshSolverSettings["problem_data"]["model_part_name"].SetString(OptimizationModelPart.Name)

        self.OptimizationModelPart = OptimizationModelPart
<<<<<<< HEAD
        print("------------------------------------------",self.OptimizationModelPart)
        self.mesh_solver = MeshMovingAnalysis(self.MeshSolverSettings, OptimizationModelPart)
=======

        model = Model()
        model.AddModelPart(self.OptimizationModelPart)

        self._mesh_moving_analysis = MeshMovingAnalysis(model, self.MeshSolverSettings)
>>>>>>> master

    # --------------------------------------------------------------------------
    def Initialize(self):
        self._mesh_moving_analysis.Initialize()

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, InputVariable):
        print("\n> Starting to update the mesh...")
        startTime = timer.time()

        VariableUtils().SetToZero_VectorVar(MESH_DISPLACEMENT,self.OptimizationModelPart.Nodes)

        sub_model_part_name = "surface_nodes"
        GeometryUtilities(self.OptimizationModelPart).ExtractBoundaryNodes(sub_model_part_name)
        surface_nodes = self.OptimizationModelPart.GetSubModelPart(sub_model_part_name).Nodes

        VariableUtils().ApplyFixity(MESH_DISPLACEMENT_X, True, surface_nodes)
        VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Y, True, surface_nodes)
        VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Z, True, surface_nodes)
        VariableUtils().CopyVectorVar(SHAPE_UPDATE, MESH_DISPLACEMENT, surface_nodes)

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
