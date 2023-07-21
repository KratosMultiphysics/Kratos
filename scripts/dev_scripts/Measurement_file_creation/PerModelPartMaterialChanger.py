from abc import ABC
from abc import abstractclassmethod

from Measurement_file_creation.MaterialChanger import MaterialChanger

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA


class PerModelPartMaterialChanger(MaterialChanger):

    def __init__(self, model_part_value_dict: "dict[str,dict[str,float]]"):
        self.model_part_value_dict = model_part_value_dict

    def adjust_material_of_model_part(self, model_part: Kratos.ModelPart) -> None:

        if not KratosOA.ModelPartUtils.CheckModelPartStatus(model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(model_part, "element_specific_properties_created")

        for model_part_name, property_value in self.model_part_value_dict.items():
            counter = 0
            for element in model_part.GetSubModelPart(model_part_name).Elements:
                for prop, value in property_value.items():
                    kratos_property = Kratos.KratosGlobals.GetVariable(prop)
                    element.Properties[kratos_property] = value
                    counter += 1

            print(f"PerModelPartMaterialChanger::Adjusted {counter} values in model part {model_part_name}")