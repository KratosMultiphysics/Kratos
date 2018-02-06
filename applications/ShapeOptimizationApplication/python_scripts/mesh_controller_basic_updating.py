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

from mesh_controller_base import MeshController

# ==============================================================================
class MeshControllerBasicUpdating( MeshController ):
    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart ):
        self.OptimizationModelPart = OptimizationModelPart

    # --------------------------------------------------------------------------
    def Initialize( self ):
        pass

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable( self, InputVariable ):
        print("\n> Starting to update the mesh")
        startTime = timer.time()
        MeshControllerUtilities( self.OptimizationModelPart ).UpdateMeshAccordingInputVariable( InputVariable )  
        print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

# ==============================================================================