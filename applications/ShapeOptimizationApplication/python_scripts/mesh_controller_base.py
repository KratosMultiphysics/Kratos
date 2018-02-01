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

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import time as timer

# ==============================================================================
class MeshController:
    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable( self, InputVariable ):
        raise RuntimeError("Mesh controller base class is called. Please check your implementation of the function >> UpdateMeshAccordingInputVariable << .")     

    # --------------------------------------------------------------------------    
    def ResetMeshDisplacement( self ):
        raise RuntimeError("Mesh controller base class is called. Please check your implementation of the function >> ResetMeshDisplacement << .")     

# ==============================================================================