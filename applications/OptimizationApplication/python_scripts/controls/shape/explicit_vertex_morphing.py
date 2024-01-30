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
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.OptimizationApplication.controls.shape.shape_control import ShapeControl

class ExplicitVertexMorphing(ShapeControl):

    def __init__(self, name, model, settings):
        self.project_to_normal = False
        self.smooth_surface = False
        self.plane_symmetry = False        
        self.plane_symmetry = False        
        super().__init__(name,model,settings)
        self.technique_settings = self.settings["technique_settings"]

    def Initialize(self):

        super().Initialize()

        self.ex_vm_mapper = {}
        for model_part_name in self.controlling_objects:
            if not self.model.HasModelPart(model_part_name):
                raise RuntimeError("ExplicitVertexMorphing: Model part {} from control {} does not exist in the input model parts".format(model_part_name,self.name))
            ex_mapper = mapper_factory.CreateMapper(self.model.GetModelPart(model_part_name), self.model.GetModelPart(model_part_name), self.technique_settings)
            ex_mapper.Initialize()
            self.ex_vm_mapper[model_part_name] = ex_mapper
    

    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        for mapper in self.ex_vm_mapper.values():
            mapper.InverseMap(derivative_variable_name,mapped_derivative_variable_name)

    def Compute(self):
        for mapper in self.ex_vm_mapper.values():
            mapper.Map(KOA.D_CX,KOA.D_X)   

    def Update(self):
        for model_part_name in self.controlling_objects:
            model_part = self.model.GetModelPart(model_part_name)
            for node in model_part.Nodes:
                shape_update = node.GetSolutionStepValue(KOA.D_X)
                node.X0 += shape_update[0]
                node.X = node.X0
                node.Y0 += shape_update[1]
                node.Y = node.Y0
                node.Z0 += shape_update[2]
                node.Z = node.Z0     

        for mapper in self.ex_vm_mapper.values():
            mapper.Update()  

    def GetControllingObjects(self):      
        return self.controlling_objects  
            


