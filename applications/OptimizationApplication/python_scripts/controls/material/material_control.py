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

class MaterialControl():
    def __init__(self, name, model, settings):
        
        self.name = name
        self.type = "material"
        self.model = model
        self.settings = settings
        self.controlling_objects = self.settings["controlling_objects"].GetStringArray() 

        self.control_variable_name = "CD"
        self.control_update_name = "D_CD"
        self.scalar_fields = ["CD","FD","PD","PE","D_CD","D_PD_D_FD","D_PE_D_FD"]
        self.vector_fields = []
        self.output_names = self.scalar_fields + self.vector_fields

        # add vars
        for model_part_name in self.controlling_objects:
            root_model = model_part_name.split(".")[0]
            for field_name in self.output_names:
                self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(field_name))


    def Initialize(self):
        # initialize zero the fields
        for model_part_name in self.controlling_objects:
            model_part = self.model.GetModelPart(model_part_name)         
            for node in model_part.Nodes:
                for vec_field in self.vector_fields:
                    node.SetSolutionStepValue(KM.KratosGlobals.GetVariable(vec_field), [0.0, 0.0, 0.0])
                for scala_field in self.scalar_fields:
                    node.SetSolutionStepValue(KM.KratosGlobals.GetVariable(scala_field), 0)
 
    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        raise RuntimeError("MaterialControl:MapFirstDerivative: calling base class function") 

    def Compute(self):
        raise RuntimeError("MaterialControl:Compute: calling base class function")

    def Update(self):
        raise RuntimeError("MaterialControl:Update: calling base class function") 

    def GetControllingObjects(self):
        raise RuntimeError("MaterialControl:GetControllingObjects: calling base class function") 

    def GetVariableName(self):
        return self.control_variable_name       

    def GetUpdateName(self):
        return self.control_update_name   

    def GetOutputNames(self):
        return self.output_names 

    def GetType(self):
        return self.type                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
            


