import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
#import KratosMultiphysics.process_factory
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics.StructuralMechanicsApplication.handbook_config_validation import Schema_Validation
from KratosMultiphysics.StructuralMechanicsApplication.structural_elements.panel import Panel


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
        Schema_Validation(settings)

        #### initialization
    def ExecuteFinalizeSolutionStep(self):
            
        for i in range(self.config_data.size()):
            structural_element = self.config_data[i]
            if structural_element["type"].GetString().lower() == "panel":
                sub_model_part = self.modelpart.GetSubModelPart(structural_element["submodelpart"].GetString())
                panel = Panel.FromKratosParametersObject(sub_model_part=sub_model_part, data=structural_element)
                print("Panel Name: ", structural_element["submodelpart"].GetString(), 
                    "\n Length: ", panel.a, 
                    "\n Width: ", panel.b, 
                    "\n Aspect Ratio: ", panel.aspect_ratio,
                    "\n E: ", panel.E,
                    "\n Nu: ", panel.nu,
                    "\n t: ", panel.thickness,
                    "\n x: ", panel.x_axis_base_vector,
                    "\n y: ", panel.y_axis_base_vector,
                    "\n z: ", panel.z_axis_base_vector,
                    "\n XX: ", panel.xx_panel_stress,
                    "\n YY: ", panel.yy_panel_stress,
                    "\n CL: ", panel.cl,
                    "\n BCs: ", panel.boundary_conditions)
                

            elif structural_element["type"].GetString().lower() == "column":
                sub_model_part = self.modelpart.GetSubModelPart(structural_element["submodelpart"].GetString())
                column = 4
                print("Column Structural Element: ", column)

            else:
                continue