import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA


class ReserveFactorResponse:
    def __init__(self, model_part, sub_model_part_name, variable=SMA.RESPONSE_VALUE):
        self.model_part = model_part
        self.sub_model_part_name = sub_model_part_name
        self.variable = variable

    def CalculateValue(self):
        if not self.model_part.HasSubModelPart(self.sub_model_part_name):
            raise RuntimeError(
                f"SubModelPart '{self.sub_model_part_name}' was not found "
                f"in ModelPart '{self.model_part.Name}'."
            )
        
        sub_model_part = self.model_part.GetSubModelPart(self.sub_model_part_name)

        if not sub_model_part.Has(self.variable):
            raise RuntimeError(
                f"SubModelPart '{self.sub_model_part_name}' has no value for "
                f"'{self.variable.Name()}'."
            )
        
        return sub_model_part.GetValue(self.variable)
