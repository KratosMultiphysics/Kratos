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

class ExplicitVertexMorphing():

    def __init__(self, name, model, settings):
        
        self.name = name
        self.model = model
        self.settings = settings
        self.technique_settings = self.settings["technique_settings"]
        self.controlling_objects = self.settings["controlling_objects"].GetStringArray() 

        # add vars
        for model_part_name in self.controlling_objects:
            root_model = model_part_name.split(".")[0]
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.SHAPE_CONTROL)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.SHAPE_CONTROL_UPDATE)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.SHAPE_UPDATE)



    def Initialize(self):
        self.ex_vm_mapper = {}
        for model_part_name in self.controlling_objects:
            if not self.model.HasModelPart(model_part_name):
                raise RuntimeError("ExplicitVertexMorphing: Model part {} from control {} does not exist in the input model parts".format(model_part_name,self.name))
            ex_mapper = mapper_factory.CreateMapper(self.model.GetModelPart(model_part_name), self.model.GetModelPart(model_part_name), self.technique_settings)
            ex_mapper.Initialize()
            self.ex_vm_mapper[model_part_name] = ex_mapper

        # initialize control fields
        zero_vec = Vector(3); 
        for model_part_name in self.controlling_objects:
            model_part = self.model.GetModelPart(model_part_name)         
            for node in model_part.Nodes:
                node.SetSolutionStepValue(KOA.SHAPE_CONTROL, zero_vec)        


        dwdw


    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        for mapper in self.ex_vm_mapper.values():
            mapper.InverseMap(derivative_variable_name,mapped_derivative_variable_name)
            


