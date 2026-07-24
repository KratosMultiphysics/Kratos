import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel import Panel

_COMPONENT_TYPES = {
    "panel": Panel
}

def CreateStructuralComponent(modelpart: KratosMultiphysics.ModelPart, component_definition: KratosMultiphysics.Parameters):

        component_type = component_definition["type"].GetString().lower()

        try:
            component_class = _COMPONENT_TYPES[component_type]
        except KeyError as exc:
            raise ValueError(
                f"Unsupported structural component type: {component_type!r}"
            ) from exc
        
        sub_model_part_name = component_definition["submodelpart"].GetString()

        if not modelpart.HasSubModelPart(sub_model_part_name):
            raise ValueError(f"SubModelPart {sub_model_part_name} was not found in ModelPart {modelpart.Name}")
        
        sub_model_part = modelpart.GetSubModelPart(sub_model_part_name)

        return component_class.FromKratosParametersObject(sub_model_part=sub_model_part, data=component_definition)
