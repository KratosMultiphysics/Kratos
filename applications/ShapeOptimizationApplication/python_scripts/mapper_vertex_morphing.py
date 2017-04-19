# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Import mapper base classes
from mapper_base import Mapper

# ==============================================================================
class MapperVertexMorphing( Mapper) :

    # --------------------------------------------------------------------------
    def __init__( self, designSurface, listOfDampingRegions, optimizationSettings ):
        self.optimizationSettings = optimizationSettings
        self.mapper = VertexMorphingMapper( designSurface, listOfDampingRegions, optimizationSettings )

    # --------------------------------------------------------------------------
    def map_to_design_space( self, constraintIsActive ):
        self.mapper.compute_mapping_matrix()
        self.mapper.map_sensitivities_to_design_space( constraintIsActive )        

    # --------------------------------------------------------------------------
    def map_to_geometry_space( self ):
        self.mapper.map_design_update_to_geometry_space()  

# ==============================================================================
