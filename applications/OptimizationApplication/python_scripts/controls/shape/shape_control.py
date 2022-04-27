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

class ShapeControl():
    def __init__(self, name, model, settings):
        
        self.name = name
        self.type = "shape"
        self.model = model
        self.settings = settings
        self.controlling_objects = self.settings["controlling_objects"].GetStringArray() 

        self.control_variable_name = "CX"
        self.control_update_name = "D_CX"
        self.output_names = ["CX","D_CX","D_X","NORMAL","AUXILIARY_FIELD"]

        # add vars
        for model_part_name in self.controlling_objects:
            root_model = model_part_name.split(".")[0]
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.CX)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.D_CX)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.D_X)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KM.NORMAL)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.AUXILIARY_FIELD)


    def Initialize(self):
        # initialize zero the fields
        for model_part_name in self.controlling_objects:
            model_part = self.model.GetModelPart(model_part_name)         
            for node in model_part.Nodes:
                node.SetSolutionStepValue(KOA.CX, [0.0, 0.0, 0.0])  
                node.SetSolutionStepValue(KOA.D_CX, [0.0, 0.0, 0.0])
                node.SetSolutionStepValue(KOA.D_X, [0.0, 0.0, 0.0])
                node.SetSolutionStepValue(KOA.AUXILIARY_FIELD, [0.0, 0.0, 0.0])  
                node.SetSolutionStepValue(KM.NORMAL, [0.0, 0.0, 0.0]) 
 
    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        raise RuntimeError("ShapeControl:MapFirstDerivative: calling base class function") 

    def Compute(self):
        raise RuntimeError("ShapeControl:Compute: calling base class function")

    def Update(self):
        raise RuntimeError("ShapeControl:Update: calling base class function") 

    def GetVariableName(self):
        return self.control_variable_name    

    def GetControllingObjects(self):
        raise RuntimeError("ShapeControl:GetControllingObjects: calling base class function")            

    def GetUpdateName(self):
        return self.control_update_name   

    def GetOutputNames(self):
        return self.output_names 

    def GetType(self):
        return self.type                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
            


