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

    def __init__(self, model_parts_list, settings):

        self.vm_mapper = {}

        if not isinstance(model_parts_list, list):
            raise RuntimeError("ExplicitVertexMorphing: Requires list of model parts")

        for model_part in model_parts_list:
            self.vm_mapper[model_part.Name] = mapper_factory.CreateMapper(model_part, model_part, settings)

    def Initialize(self):
        for key,value in self.vm_mapper.items():
            value.Initialize()


