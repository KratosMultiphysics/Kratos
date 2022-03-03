# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory

class ExplicitVertexMorphing():

    def __init__(self, name, model, settings):
        
        self.name = name
        self.model = model
        self.settings = settings
        self.technique_settings = self.settings["technique_settings"]

        default_technique_settings = KM.Parameters("""
        {
            "filter_function_type"       : "linear",
            "filter_radius"              : 0.000000000001,
            "max_nodes_in_filter_radius" : 10000,
            "matrix_free_filtering"      : false,
            "consistent_mapping"         : false,
            "improved_integration"       : false,
            "integration_method"         : "gauss_integration",
            "number_of_gauss_points"     : 5,
            "in_plane_morphing"          : false,
            "in_plane_morphing_settings" : {},
            "sliding_morphing"           : false,
            "sliding_morphing_settings"  : {},
            "plane_symmetry"             : false,
            "plane_symmetry_settings"    : {
                "point" : [0.0, 0.0, 0.0],
                "normal": [0.0, 0.0, 0.0]
            },
            "revolution"                 : false,
            "revolution_settings"        : {
                "point" : [0.0, 0.0, 0.0],
                "normal": [0.0, 0.0, 0.0]
            }
        }""")

        self.settings["technique_settings"].ValidateAndAssignDefaults(default_technique_settings)   

        self.control = KOA.ExplicitVertexMorphing(self.name,self.model,self.settings)   

  
        # if not isinstance(model_parts_names_list, list):
        #     raise RuntimeError("ExplicitVertexMorphing: Requires list of model part names")
        

    def Initialize(self):
        self.control.Initialize()
        # self.ex_vm_mapper = {}
        # for model_part_name in self.model_parts_names:
        #     if not self.model.HasModelPart(model_part_name):
        #         raise RuntimeError("ExplicitVertexMorphing: Model part {} from control {} does not exist in the input model parts".format(model_part_name,self.name))
        #     ex_mapper = mapper_factory.CreateMapper(self.model.GetModelPart(model_part_name), self.model.GetModelPart(model_part_name), self.settings)
        #     ex_mapper.Initialize()
        #     self.ex_vm_mapper[model_part_name] = ex_mapper


