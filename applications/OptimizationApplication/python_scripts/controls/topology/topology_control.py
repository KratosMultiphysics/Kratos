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

class TopologyControl():
    def __init__(self, name, model, settings):
        
        self.name = name
        self.type = "topology"
        self.model = model
        self.settings = settings
        self.controlling_objects = self.settings["controlling_objects"].GetStringArray() 

        self.control_variable_name = "CD"
        self.control_update_name = "D_CD"
        self.output_names = ["CD","PD","FD","D_CD"]

        # add vars
        for model_part_name in self.controlling_objects:
            root_model = model_part_name.split(".")[0]
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.CD)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.PD)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KM.DENSITY)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.FD)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.D_CD)


    def Initialize(self):
        # initialize zero the fields
        for model_part_name in self.controlling_objects:
            model_part = self.model.GetModelPart(model_part_name)         
            for node in model_part.Nodes:
                node.SetSolutionStepValue(KOA.CD, 0.0)  
                node.SetSolutionStepValue(KOA.PD, 0.0)
                node.SetSolutionStepValue(KM.DENSITY, 1.0)
                node.SetSolutionStepValue(KOA.FD, 0.0)
                node.SetSolutionStepValue(KOA.D_CD, 0.0)
 
    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        raise RuntimeError("TopologyControl:MapFirstDerivative: calling base class function") 

    def Compute(self):
        raise RuntimeError("TopologyControl:Compute: calling base class function")

    def Update(self):
        raise RuntimeError("TopologyControl:Update: calling base class function") 

    def GetControllingObjects(self):
        raise RuntimeError("TopologyControl:GetControllingObjects: calling base class function") 

    def GetVariableName(self):
        return self.control_variable_name       

    def GetUpdateName(self):
        return self.control_update_name   

    def GetOutputNames(self):
        return self.output_names 

    def GetType(self):
        return self.type                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
            


