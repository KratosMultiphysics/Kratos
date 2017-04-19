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

from mapper_vertex_morphing import MapperVertexMorphing

# ==============================================================================
def CreateMapper( designSurface, listOfDampingRegions, optimizationSettings ):

    design_variables_type = optimizationSettings["design_variables"]["design_variables_type"].GetString()
    
    if design_variables_type == "vertex_morphing":
        return MapperVertexMorphing( designSurface, listOfDampingRegions, optimizationSettings ) 
    else:
        raise NameError("The following design variables type is not supported by the mapper (name may be misspelled): " + design_variables_type)              

# ==============================================================================