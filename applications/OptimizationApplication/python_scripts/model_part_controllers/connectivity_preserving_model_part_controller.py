import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"ConnectivityPreservingModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return ConnectivityPreservingModelPartController(model, parameters["settings"])

class ConnectivityPreservingModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "",
                    "destination_model_part_name": "",
                    "destination_element_name"   : "",
                    "destination_condition_name" : ""
                }
            ]
        }""")

        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        for sub_model_part_settings in self.parameters["transformation_settings"].values():
            sub_model_part_settings.ValidateAndAssignDefaults(default_settings["transformation_settings"].values()[0])

    def ImportModelPart(self) -> None:
        connectivity_preserve_modeller = Kratos.ConnectivityPreserveModeler()
        for transformation_settings in self.parameters["transformation_settings"].values():
            # get source and destination model parts
            source_model_part = self.model[transformation_settings["source_model_part_name"].GetString()]
            destination_model_part_name: str = transformation_settings["destination_model_part_name"].GetString()
            if not self.model.HasModelPart(destination_model_part_name):
                # if the destination model part is not there, then check whether the destination root model part is existing
                destination_root_model_part_name = destination_model_part_name.split(".")[0]
                if not self.model.HasModelPart(destination_root_model_part_name):
                    # destination root model part is also not there. so we can safely assign the source model part solution variables list
                    destination_root_model_part = self.model.CreateModelPart(destination_root_model_part_name)
                else:
                    # if it is existing, then check whether it has a subset of the variables list.
                    destination_root_model_part = self.model[destination_root_model_part_name]
                    if not KratosOA.OptimizationUtils.IsSolutionStepVariablesListASubSet(source_model_part, destination_root_model_part):
                        raise RuntimeError(f"The destination root model part solution step variables at {destination_root_model_part.FullName()} is not a sub set of {source_model_part.FullName()}")
                KratosOA.OptimizationUtils.CopySolutionStepVariablesList(destination_root_model_part, source_model_part)
                destination_model_part = self.model.CreateModelPart(destination_model_part_name)
            else:
                destination_model_part = self.model[destination_model_part_name]

            element_name = transformation_settings["destination_element_name"].GetString()
            condition_name = transformation_settings["destination_condition_name"].GetString()

            if condition_name != "" and element_name != "":
                connectivity_preserve_modeller.GenerateModelPart(source_model_part, destination_model_part, element_name, condition_name)
            elif element_name != "":
                connectivity_preserve_modeller.GenerateModelPart(source_model_part, destination_model_part, element_name)
            elif condition_name != "":
                connectivity_preserve_modeller.GenerateModelPartWithConditions(source_model_part, destination_model_part, condition_name)

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Duplicated {source_model_part.FullName()} in to {destination_model_part.FullName()}.")

    def GetModelPart(self) -> Kratos.ModelPart:
        return None