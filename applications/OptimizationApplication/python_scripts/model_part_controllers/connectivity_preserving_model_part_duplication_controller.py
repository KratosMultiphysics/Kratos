import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"ConnectivityPreservingModelPartDuplicationController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return ConnectivityPreservingModelPartDuplicationController(model, parameters["settings"])

class ConnectivityPreservingModelPartDuplicationController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "source_model_part_name"      : "",
            "destination_model_part_name" : "",
            "destination_element_name"    : "",
            "destination_condition_name"  : ""
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.source_model_part = model[parameters["source_model_part_name"].GetString()]
        destination_model_part_name = parameters["destination_model_part_name"].GetString()

        if not model.HasModelPart(destination_model_part_name):
            self.destination_model_part = model.CreateModelPart(destination_model_part_name)

        self.destination_element_name = parameters["destination_element_name"].GetString()
        self.destination_condition_name = parameters["destination_condition_name"].GetString()

        # set only the domain size just in case if a solver requires this to properly identify the element type.
        # all the other data in model part will be copied when we do the copying later at ImportModelPart
        self.destination_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = self.source_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

    def ImportModelPart(self) -> None:
        connectivity_preserve_modeller = Kratos.ConnectivityPreserveModeler()
        connectivity_preserve_modeller.GenerateModelPart(self.source_model_part, self.destination_model_part, self.destination_element_name, self.destination_condition_name)
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Duplicated {self.source_model_part.FullName()}.")

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.destination_model_part