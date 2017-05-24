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

# ==============================================================================
def CreateMapper( designSurface, optimizationSettings ):
    isMatrixFreeMappingRequired = optimizationSettings["design_variables"]["filter"]["matrix_free_filtering"].GetBool()
    if isMatrixFreeMappingRequired:
        return MapperVertexMorphingMatrixFree( designSurface, optimizationSettings )  
    try:
        isIntegrationMappingRequired = optimizationSettings["design_variables"]["integration"].GetBool()
        if isIntegrationMappingRequired:
            if isIterativeMappingRequired:
                raise ValueError ("Iterative Mapper and Integration cannot be combined yet!")
            else:
                return MapperVertexMorphingIntegration( designSurface, listOfDampingRegions, optimizationSettings )
    except:
        pass

    else:
        return MapperVertexMorphing( designSurface, optimizationSettings )         

# ==============================================================================