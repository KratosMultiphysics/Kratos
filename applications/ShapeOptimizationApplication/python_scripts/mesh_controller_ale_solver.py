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

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import time as timer

from mesh_controller_base import MeshController

# # ==============================================================================
class MeshControllerUsingALESolver( MeshController) :
    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, MeshMotionSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.MeshMotionSettings = MeshMotionSettings

        mesh_solver_module = __import__(self.MeshMotionSettings["mesh_solver_settings"]["solver_type"].GetString())
        self.mesh_solver = mesh_solver_module.CreateSolver(self.OptimizationModelPart, self.MeshMotionSettings["mesh_solver_settings"])

        self.mesh_solver.AddVariables()

    # --------------------------------------------------------------------------
    def Initialize( self ):
        self.mesh_solver.AddDofs()
        self.mesh_solver.Initialize()
        self.mesh_solver.SetEchoLevel(0)

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable( self, InputVariable ):
        print("\n> Starting to update the mesh...")
        startTime = timer.time()

        VariableUtils().SetToZero_VectorVar(MESH_DISPLACEMENT,self.OptimizationModelPart.Nodes)

        sub_model_part_name = "surface_nodes"
        GeometryUtilities(self.OptimizationModelPart).ExtractSurfaceNodes(sub_model_part_name)
        surface_nodes = self.OptimizationModelPart.GetSubModelPart(sub_model_part_name).Nodes

        VariableUtils().ApplyFixity(MESH_DISPLACEMENT_X, True, surface_nodes)
        VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Y, True, surface_nodes)
        VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Z, True, surface_nodes)
        VariableUtils().CopyVectorVar(SHAPE_UPDATE, MESH_DISPLACEMENT, surface_nodes)

        self.mesh_solver.Solve()

        MeshControllerUtilities( self.OptimizationModelPart ).LogMeshChangeAccordingInputVariable( MESH_DISPLACEMENT )

        print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

# ==============================================================================