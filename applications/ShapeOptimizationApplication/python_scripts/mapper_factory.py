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

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# ==============================================================================
def CreateMapper(origin_model_part, destination_model_part, mapper_settings):
    default_settings = Parameters("""
    {
        "filter_function_type"        : "linear",
        "filter_radius"               : 0.000000000001,
        "max_nodes_in_filter_radius"  : 10000,
        "apply_matrix_free_filtering" : false,
        "apply_consistent_mapping"    : false,
        "apply_improved_integration"  : false,
        "integration_method"          : "gauss_integration",
        "number_of_gauss_points"      : 5
    }""")

    mapper_settings.RecursivelyValidateAndAssignDefaults(default_settings)

    if mapper_settings["apply_matrix_free_filtering"].GetBool():
        if mapper_settings["apply_consistent_mapping"].GetBool():
             raise ValueError ("Matrix free mapper has no option to map consistently yet!")
        if mapper_settings["apply_improved_integration"].GetBool():
             raise ValueError ("Matrix free mapper does not yet allow for an imporoved integration!")
        else:
            return MapperVertexMorphingMatrixFree(origin_model_part, destination_model_part, mapper_settings)
    else:
        if mapper_settings["apply_improved_integration"].GetBool():
            return MapperVertexMorphingImprovedIntegration(origin_model_part, destination_model_part, mapper_settings)
        else:
            return MapperVertexMorphing(origin_model_part, destination_model_part, mapper_settings)

# ==============================================================================