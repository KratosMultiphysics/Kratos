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
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
import time as timer
from mesh_controller_base import MeshController

# ==============================================================================
class MeshControllerBasicUpdating(MeshController):
    # --------------------------------------------------------------------------
    def __init__(self, OptimizationModelPart):
        self.OptimizationModelPart = OptimizationModelPart

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, variable, design_surface):
        print("\n> Starting to update the mesh")
        startTime = timer.time()
        MeshControllerUtilities(self.OptimizationModelPart).UpdateMeshAccordingInputVariable(variable)
        MeshControllerUtilities(self.OptimizationModelPart).LogMeshChangeAccordingInputVariable(variable)
        print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

# ==============================================================================