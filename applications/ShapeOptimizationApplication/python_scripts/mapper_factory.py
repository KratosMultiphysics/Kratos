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
def CreateMapper(model_part, mapper_settings):
    default_settings = Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 1.0,
        "max_nodes_in_filter_radius" : 10000,
        "matrix_free_filtering"      : false,
        "integration":
        {
            "integration_method": "node_sum",
            "number_of_gauss_points": 0
        },
        "consistent_mapping_to_geometry_space": false
    }""")

    mapper_settings.RecursivelyValidateAndAssignDefaults(default_settings)

    if mapper_settings["matrix_free_filtering"].GetBool():
        if mapper_settings["consistent_mapping_to_geometry_space"].GetBool():
             raise ValueError ("Matrix free Mapper has no consistent_mapping_to_geometry_space option yet!")
        if mapper_settings["integration"]["integration_method"].GetString() != "node_sum":
             raise ValueError ("Matrix free Mapper can only be combined with 'node_sum' integration method!")
        else:
            return MapperVertexMorphingMatrixFree(model_part, mapper_settings)
    else:
        if mapper_settings["integration"]["integration_method"].GetString() in ["gauss_integration", "area_weighted_sum"]:
            return MapperVertexMorphingImprovedIntegration(model_part, mapper_settings)
        elif mapper_settings["integration"]["integration_method"].GetString() == "node_sum":
            return MapperVertexMorphing(model_part, mapper_settings)
        else:
            raise ValueError ("CreateMapper: integration_method not known!")

# ==============================================================================