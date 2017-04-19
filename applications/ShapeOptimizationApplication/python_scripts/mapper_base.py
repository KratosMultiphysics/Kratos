# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
class Mapper:

    # --------------------------------------------------------------------------
    def map_to_design_space( self, variable ):
        raise RuntimeError("Mapper base class is called. Please check your implementation of the function >> map_to_design_space << .")     

    # --------------------------------------------------------------------------
    def map_to_geometry_space( self, variable ):
        raise RuntimeError("Mapper base class is called. Please check your implementation of the function >> map_to_geometry_space << .")             

# ==============================================================================
