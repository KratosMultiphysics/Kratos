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

# ==============================================================================
def CreateMapper( designSurface, optimizationSettings ):

    # default settings string in json format
    default_settings = Parameters("""
    {
        "design_variables_type"     : "vertex_morphing",
        "input_model_part_name"     : "input_model_part_name",
        "design_submodel_part_name" : "design_surface",
        "domain_size"               : 3,
        "filter" : {
            "filter_function_type"       : "linear",
            "filter_radius"              : 1,
            "max_nodes_in_filter_radius" : 10000,
            "matrix_free_filtering"      : false
        },
        "consistent_mapping_to_geometry_space": false,
        "align_sensitivities_with_normal" : true,
        "integration": {
            "integration_method": "node_sum",
            "number_of_gauss_points": 0
        },
        "damping" : {
            "perform_damping" : false,
            "damping_regions" : []
        }
    }""")

    # overwrite the default settings with user-provided parameters
    optimizationSettings["design_variables"].ValidateAndAssignDefaults(default_settings)

    isMatrixFreeMappingRequired = optimizationSettings["design_variables"]["filter"]["matrix_free_filtering"].GetBool()
    if isMatrixFreeMappingRequired:
        return MapperVertexMorphingMatrixFree( designSurface, optimizationSettings )
    else:
        return MapperVertexMorphing( designSurface, optimizationSettings )

# ==============================================================================