import numpy as np
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.handbook_config_validation import Schema_Validation
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.structural_component_factory import CreateStructuralComponent
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.handbook_method_factory import CreateHandbookMethod
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.analysis_result import AnalysisResult

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return HandbookAnalysisProcess(model, settings["Parameters"])

class HandbookAnalysisProcess(KratosMultiphysics.Process):
    """Runs handbook methods on structural components and stores the reserve factors."""

    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters):
        super().__init__()
        Schema_Validation(settings)
        self.model_part = model.GetModelPart("Structure")
        self.component_definitions = settings["Structural_Components"]
        self.evaluations = []
        self.results = {}

    def ExecuteBeforeSolutionLoop(self):
        for definition in self.component_definitions.values():
            component = CreateStructuralComponent(self.model_part, definition)
            methods = [CreateHandbookMethod(name.GetString()) for name in definition["analysis_methods"].values()]
            for method in methods:
                method.CheckCompatibility(component)
            component.Initialize()
            self.evaluations.append((component, methods))

    def ExecuteFinalizeSolutionStep(self):
        self.results = {}
        for component, methods in self.evaluations:
            component.PrepareAnalysis()
            for method in methods:
                self.results[(component.name, method.name)] = (
                    method.Evaluate(component) if method.IsApplicable(component)
                    else AnalysisResult(method.name, method.category, np.inf))

        print(self.GetReserveFactors())

    def GetReserveFactors(self) -> dict[str, float]:
        """Fixed keys '<component>.<method>'; inf if not applicable, never NaN."""
        return {f"{component_name}.{method_name}": result.value
                for (component_name, method_name), result in self.results.items()}
