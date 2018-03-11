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

# ==============================================================================
class MeshController:
    # --------------------------------------------------------------------------
    def Initialize( self ):
        raise NotImplementedError("Mesh controller base class is called. Please check your implementation of the function >> Initialize << .")     

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable( self, InputVariable ):
        raise NotImplementedError("Mesh controller base class is called. Please check your implementation of the function >> UpdateMeshAccordingInputVariable << .")     

# ==============================================================================
