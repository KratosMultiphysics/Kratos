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
        "filter_radius_factor"       : 7.0,
        "max_nodes_in_filter_radius" : 10000,
        "matrix_free_filtering"      : false,
        "consistent_mapping"         : false,
        "improved_integration"       : false,
        "integration_method"         : "gauss_integration",
        "number_of_gauss_points"     : 5,
        "in_plane_morphing"          : false,
        "in_plane_morphing_settings" : {},
        "plane_symmetry"             : false,
        "plane_symmetry_settings"    : {
            "point" : [0.0, 0.0, 0.0],
            "normal": [1.0, 0.0, 0.0]
        },
        "revolution"             : false,
        "revolution_settings"    : {
            "point" : [0.0, 0.0, 0.0],
            "normal": [0.0, 0.0, 1.0]
        }
    }""")

    mapper_vertex_morphing_matrix_free = KSO.MapperVertexMorphingMatrixFree
    mapper_vertex_morphing_improved_integration = KSO.MapperVertexMorphingImprovedIntegration
    mapper_vertex_morphing_symmetric = KSO.MapperVertexMorphingSymmetric
    mapper_vertex_morphing = KSO.MapperVertexMorphing
    if mapper_settings.Has("filter_radius"):
        if mapper_settings["filter_radius"].IsString():
            if mapper_settings["filter_radius"].GetString() == "adaptive":
                mapper_vertex_morphing_matrix_free = KSO.MapperVertexMorphingMatrixFreeAdaptiveRadius
                mapper_vertex_morphing_improved_integration = KSO.MapperVertexMorphingImprovedIntegrationAdaptiveRadius
                mapper_vertex_morphing_symmetric = KSO.MapperVertexMorphingSymmetricAdaptiveRadius
                mapper_vertex_morphing = KSO.MapperVertexMorphingAdaptiveRadius

                if mapper_settings["in_plane_morphing"].GetBool():
                    raise Exception("\"in_plane_morphing\" is not yet supported with \"adaptive\" filter radius.")
                mapper_settings["filter_radius"].SetDouble(-1.0)
            else:
                raise Exception("\"filter_radius\" either should be double value or \"adaptive\".")

    mapper_settings.ValidateAndAssignDefaults(default_settings)

    if mapper_settings["in_plane_morphing"].GetBool():
        return in_plane_vertex_morphing_mapper.InPlaneVertexMorphingMapper(origin_model_part, destination_model_part, mapper_settings)
    elif mapper_settings["matrix_free_filtering"].GetBool():
        if mapper_settings["consistent_mapping"].GetBool():
             raise ValueError ("Matrix free mapper has no option to map consistently yet!")
        if mapper_settings["improved_integration"].GetBool():
             raise ValueError ("Matrix free mapper does not yet allow for an improved integration!")
        else:
            return mapper_vertex_morphing_matrix_free(origin_model_part, destination_model_part, mapper_settings)
    else:
        if mapper_settings["revolution"].GetBool() and mapper_settings["plane_symmetry"].GetBool():
            raise RuntimeError("revolution and plane_symmetry can not be combined!")

        if mapper_settings["improved_integration"].GetBool():
            return mapper_vertex_morphing_improved_integration(origin_model_part, destination_model_part, mapper_settings)
        elif mapper_settings["plane_symmetry"].GetBool():
            return mapper_vertex_morphing_symmetric(origin_model_part, destination_model_part, mapper_settings)
        elif mapper_settings["revolution"].GetBool():
            return mapper_vertex_morphing_symmetric(origin_model_part, destination_model_part, mapper_settings)
        else:
            return mapper_vertex_morphing(origin_model_part, destination_model_part, mapper_settings)

# ==============================================================================