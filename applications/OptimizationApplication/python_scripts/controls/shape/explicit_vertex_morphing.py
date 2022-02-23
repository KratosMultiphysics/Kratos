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
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory

class ExplicitVertexMorphing():

    def __init__(self, name, model, model_parts_names_list, settings):
        
        self.name = name
        self.model = model
        self.model_parts_names = model_parts_names_list
        self.settings = settings

        if not isinstance(model_parts_names_list, list):
            raise RuntimeError("ExplicitVertexMorphing: Requires list of model part names")
        

    def Initialize(self):
        self.ex_vm_mapper = {}
        for model_part_name in self.model_parts_names:
            if not self.model.HasModelPart(model_part_name):
                raise RuntimeError("ExplicitVertexMorphing: Model part {} from control {} does not exist in the input model parts".format(model_part_name,self.name))
            ex_mapper = mapper_factory.CreateMapper(self.model.GetModelPart(model_part_name), self.model.GetModelPart(model_part_name), self.settings)
            ex_mapper.Initialize()
            self.ex_vm_mapper[model_part_name] = ex_mapper


