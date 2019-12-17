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
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from .mapping import in_plane_vertex_morphing_mapper

# ==============================================================================
def CreateMapper(origin_model_part, destination_model_part, mapper_settings):
    default_settings = KM.Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 0.000000000001,
        "max_nodes_in_filter_radius" : 10000,
        "matrix_free_filtering"      : false,
        "consistent_mapping"         : false,
        "improved_integration"       : false,
        "integration_method"         : "gauss_integration",
        "number_of_gauss_points"     : 5,
        "in_plane_morphing"                   : false,
        "in_plane_morphing_settings"          : {}
    }""")

    mapper_settings.ValidateAndAssignDefaults(default_settings)

    if mapper_settings["in_plane_morphing"].GetBool():
        return in_plane_vertex_morphing_mapper.InPlaneVertexMorphingMapper(origin_model_part, destination_model_part, mapper_settings)
    elif mapper_settings["matrix_free_filtering"].GetBool():
        if mapper_settings["consistent_mapping"].GetBool():
             raise ValueError ("Matrix free mapper has no option to map consistently yet!")
        if mapper_settings["improved_integration"].GetBool():
             raise ValueError ("Matrix free mapper does not yet allow for an improved integration!")
        else:
            return KSO.MapperVertexMorphingMatrixFree(origin_model_part, destination_model_part, mapper_settings)
    else:
        if mapper_settings["improved_integration"].GetBool():
            return KSO.MapperVertexMorphingImprovedIntegration(origin_model_part, destination_model_part, mapper_settings)
        else:
            return KSO.MapperVertexMorphing(origin_model_part, destination_model_part, mapper_settings)

# ==============================================================================