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
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
import time as timer
from .mesh_controller_base import MeshController

# ==============================================================================
class MeshControllerBasicUpdating(MeshController):
    # --------------------------------------------------------------------------
    def __init__(self, OptimizationModelPart):
        self.OptimizationModelPart = OptimizationModelPart

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, variable):
        KM.Logger.Print("")
        KM.Logger.PrintInfo("ShapeOpt", "Starting to update the mesh")
        startTime = timer.time()
        KSO.MeshControllerUtilities(self.OptimizationModelPart).UpdateMeshAccordingInputVariable(variable)
        KSO.MeshControllerUtilities(self.OptimizationModelPart).LogMeshChangeAccordingInputVariable(variable)
        KM.Logger.PrintInfo("ShapeOpt", "Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

# ==============================================================================