import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
#import KratosMultiphysics.process_factory
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics.StructuralMechanicsApplication.handbook_config_validation import Schema_Validation
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.structural_component_factory import CreateStructuralComponent
from KratosMultiphysics.StructuralMechanicsApplication.reserve_factor_response import ReserveFactorResponse

def Factory(settings: KratosMultiphysics.Parameters, Model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    
    return StressAnalysisProcess(Model, settings["Parameters"])

class StressAnalysisProcess(KratosMultiphysics.Process):
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters):
        KratosMultiphysics.Process.__init__(self)

        Schema_Validation(settings)

        self.model = model
        self.modelpart = model.GetModelPart("Structure")
        self.structural_component_definitions = settings["Structural_Elements"]

    def ExecuteBeforeSolutionLoop(self):

        self.structural_components = [CreateStructuralComponent(self.modelpart, component_definition) for component_definition in self.structural_component_definitions.values()]
        for component in self.structural_components:
            component.Initialize()
        
    def ExecuteFinalizeSolutionStep(self):
            
        for component in self.structural_components:
            component.PrepareAnalysis()
            component.RunAnalysis()

            response = ReserveFactorResponse(self.modelpart, component.sub_model_part.Name)
            rf = response.CalculateValue()

            KratosMultiphysics.Logger.PrintInfo(
                "ReserveFactorResponse",
                f"{component.sub_model_part.Name}: RF={rf:.6e}"
            )
