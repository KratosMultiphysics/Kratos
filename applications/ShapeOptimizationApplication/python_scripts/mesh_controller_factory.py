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
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import time as timer

# ==============================================================================
def CreateMeshController( OptimizationModelPart, OptimizationSettings ):
    MeshMotionSettings = OptimizationSettings["design_variables"]["mesh_motion"]
    if MeshMotionSettings["apply_ale_mesh_solver"].GetBool():
        from mesh_controller_ale_solver import MeshControllerUsingALESolver
        return MeshControllerUsingALESolver( OptimizationModelPart, MeshMotionSettings)
    else:
        from mesh_controller_basic_updating import MeshControllerBasicUpdating
        return MeshControllerBasicUpdating( OptimizationModelPart )

# # ==============================================================================