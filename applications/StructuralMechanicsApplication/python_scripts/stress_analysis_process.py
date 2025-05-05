import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.process_factory
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import json
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.structural_elements import Panel
from jsonschema import validate, ValidationError
from KratosMultiphysics.StructuralMechanicsApplication.handbook_config_validation import Schema_Validation

def Factory(settings: KratosMultiphysics.Parameters, Model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    
    return StressAnalysisProcess(Model, settings["Parameters"])

class StressAnalysisProcess(KratosMultiphysics.Process):
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters):
        KratosMultiphysics.Process.__init__(self)

        self.model = model
        self.modelpart = model.GetModelPart("Structure")
        self.config_data = settings["Structural_Elements"]
        #Schema_Validation(settings)

        #### initialization
    def ExecuteFinalizeSolutionStep(self):
            
        for i in range(self.config_data.size()):
            structural_element = self.config_data[i]
            match structural_element["type"].GetString():
                case "Panel":
                    sub_model_part = self.modelpart.GetSubModelPart(structural_element["submodelpart"].GetString())
                    panel = Panel.FromKratosParametersObject(sub_model_part=sub_model_part, data=structural_element)
                    print("Panel Name: ", structural_element["submodelpart"].GetString(), 
                      "\n Length: ", panel.length, 
                      "\n Width: ", panel.width, 
                      "\n Ratio: ", panel.aspect_ratio,
                      "\n E: ", panel.E,
                      "\n Nu: ", panel.nu,
                      "\n t: ", panel.thickness,
                      "\n x: ", panel.x_axis_base_vector,
                      "\n y: ", panel.y_axis_base_vector,
                      "\n z: ", panel.z_axis_base_vector,
                      "\n XX: ", panel.xx_panel_stress,
                      "\n YY: ", panel.yy_panel_stress,
                      "\n CL: ", panel.cl)

                case _:
                    print(f"{structural_element['type']} is not a valid structural element type.")