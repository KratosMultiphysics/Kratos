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
        self.MeshControllerUtilities = MeshControllerUtilities( self.OptimizationModelPart )

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
        print("\n> Starting to update the mesh")
        startTime = timer.time()

        # Extract surface nodes
        sub_model_part_of_surface = "surface_nodes"     
        GeometryUtilities(self.OptimizationModelPart).ExtractSurfaceNodes(sub_model_part_of_surface)

        # Apply shape update as boundary condition for computation of mesh displacement 
        for node in self.OptimizationModelPart.GetSubModelPart(sub_model_part_of_surface).Nodes:
            node.Fix(MESH_DISPLACEMENT_X)
            node.Fix(MESH_DISPLACEMENT_Y)
            node.Fix(MESH_DISPLACEMENT_Z)              
            disp = Vector(3)
            disp[0] = node.GetSolutionStepValue(SHAPE_UPDATE_X)
            disp[1] = node.GetSolutionStepValue(SHAPE_UPDATE_Y)
            disp[2] = node.GetSolutionStepValue(SHAPE_UPDATE_Z)
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,disp)

        # Solve for mesh-update
        self.mesh_solver.Solve()

        # Update reference mesh (Since shape updates are imposed as incremental quantities)
        self.mesh_solver.get_mesh_motion_solver().UpdateReferenceMesh()

        # Log absolute mesh displacement
        for node in self.OptimizationModelPart.Nodes:
            mesh_change = Vector(3)
            mesh_change[0] = node.GetSolutionStepValue(MESH_CHANGE_X) + node.GetSolutionStepValue(MESH_DISPLACEMENT_X)
            mesh_change[1] = node.GetSolutionStepValue(MESH_CHANGE_Y) + node.GetSolutionStepValue(MESH_DISPLACEMENT_Y)
            mesh_change[2] = node.GetSolutionStepValue(MESH_CHANGE_Z) + node.GetSolutionStepValue(MESH_DISPLACEMENT_Z)
            node.SetSolutionStepValue(MESH_CHANGE,0,mesh_change)     

        print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

        err

    # --------------------------------------------------------------------------    
    def ResetMeshDisplacement( self ):
        self.MeshControllerUtilities.ResetMeshDisplacement()        

# ==============================================================================