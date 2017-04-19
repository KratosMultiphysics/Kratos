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
def CreateMapper( designSurface, listOfDampingRegions, optimizationSettings ):
    isIterativeMappingRequired = optimizationSettings["design_variables"]["iterative_mapping"].GetBool()
    if isIterativeMappingRequired:
        raise NameError("Iterative mapping is not yet supported by the mapper!","")              
    else:
        return MapperVertexMorphing( designSurface, listOfDampingRegions, optimizationSettings )         

# ==============================================================================