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
        "design_variables_type"              : "vertex_morphing",
        "optimization_model_part_name"       : "optimization_model_part_name",
        "design_surface_sub_model_part_name" : "design_surface",
        "domain_size"                  : 3,
        "filter" : {
            "filter_function_type"       : "linear",
            "filter_radius"              : 1.0,
            "max_nodes_in_filter_radius" : 10000,
            "matrix_free_filtering"      : false,
            "integration": {
                "integration_method": "node_sum",
                "number_of_gauss_points": 0
            }
        },
        "consistent_mapping_to_geometry_space": false,
        "damping" : {
            "perform_damping" : false,
            "damping_regions" : []
        }
    }""")

    # overwrite the default settings with user-provided parameters
    optimizationSettings["design_variables"].RecursivelyValidateAndAssignDefaults(default_settings)

    isMatrixFreeMappingRequired = optimizationSettings["design_variables"]["filter"]["matrix_free_filtering"].GetBool()
    integrationMethod = optimizationSettings["design_variables"]["filter"]["integration"]["integration_method"].GetString()
    if isMatrixFreeMappingRequired:
        if optimizationSettings["design_variables"]["consistent_mapping_to_geometry_space"].GetBool():
             raise ValueError ("Matrix free Mapper has no consistent_mapping_to_geometry_space option yet!")
        if integrationMethod != "node_sum":
             raise ValueError ("Matrix free Mapper can only be combined with 'node_sum' integration method!")
        else:
            return MapperVertexMorphingMatrixFree( designSurface, optimizationSettings )
    else:
        if integrationMethod in ["gauss_integration", "area_weighted_sum"]:
            return MapperVertexMorphingImprovedIntegration( designSurface, optimizationSettings )
        elif integrationMethod == "node_sum":
            return MapperVertexMorphing( designSurface, optimizationSettings )
        else:
            raise ValueError ("CreateMapper: integration_method not known!")

# ==============================================================================