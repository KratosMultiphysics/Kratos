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

class ThicknessControl():
    def __init__(self, name, model, settings):
        
        self.name = name
        self.type = "thickness"
        self.model = model
        self.settings = settings
        self.controlling_objects = self.settings["controlling_objects"].GetStringArray() 

        self.control_variable_name = "CT"
        self.control_update_name = "D_CT"
        self.output_names = ["CT","PT","D_CT","D_PT"]

        # add vars
        for model_part_name in self.controlling_objects:
            root_model = model_part_name.split(".")[0]
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.CT)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.PT)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.D_CT)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.D_PT)


    def Initialize(self):
        # initialize zero the fields
        for model_part_name in self.controlling_objects:
            model_part = self.model.GetModelPart(model_part_name)         
            for node in model_part.Nodes:
                node.SetSolutionStepValue(KOA.CT, 0.0)  
                node.SetSolutionStepValue(KOA.PT, 0.0)
                node.SetSolutionStepValue(KOA.D_CT, 0.0) 
                node.SetSolutionStepValue(KOA.D_PT, 0.0) 
 
    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        raise RuntimeError("ThicknessControl:MapFirstDerivative: calling base class function") 

    def Compute(self):
        raise RuntimeError("ThicknessControl:Compute: calling base class function")

    def Update(self):
        raise RuntimeError("ThicknessControl:Update: calling base class function") 

    def GetVariableName(self):
        return self.control_variable_name       

    def GetUpdateName(self):
        return self.control_update_name   

    def GetOutputNames(self):
        return self.output_names 

    def GetType(self):
        return self.control_update_name                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
            


