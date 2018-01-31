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

# ==============================================================================
def CreateMeshController( OptimizationModelPart, OptimizationSettings ):
    return MeshController( OptimizationModelPart, OptimizationSettings )

# # ==============================================================================

class MeshController:
    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, OptimizationSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.MeshMotionSettings = OptimizationSettings["design_variables"]["mesh_motion"]
        self.MeshControllerUtilities = MeshControllerUtilities( self.OptimizationModelPart )
        # if self.__IsInternalMeshSolverSpecified():
        #     self.MeshMotionSolver.AddVariables() 

    # --------------------------------------------------------------------------
    def Initialize( self ):
        pass

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable( self, InputVariable ):
        self.MeshControllerUtilities.UpdateMeshAccordingInputVariable( InputVariable )  

    # --------------------------------------------------------------------------    
    def SetCurrentMeshAsNewReference( self ):
        self.MeshControllerUtilities.SetCurrentMeshAsNewReference()        

    # --------------------------------------------------------------------------
    def Finalize( self ):
        pass

    # # --------------------------------------------------------------------------
    # def __IsInternalMeshSolverSpecified( self ):
    #     return self.OptimizationSettings["design_variables"]["mesh_motion"]["move_mesh_according_to_shape_update"].GetBool()


    #         mesh_solver_module = __import__( MeshMotionSettings["solver_type"].GetString() )
    # return mesh_solver_module.CreateSolver( OptimizationModelPart, MeshMotionSettings )

# ==============================================================================