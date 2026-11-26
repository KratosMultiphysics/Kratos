import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyConstraintVectorDamTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyConstraintVectorDamTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        # The Dam wrapper deliberately ignores "interval": the component
        # processes keep their standard default interval (the table process does
        # not use it either).
        self.components_process_list = []

        for i, suffix in enumerate("XYZ"):
            component_params = KratosMultiphysics.Parameters("{}")
            component_params.AddValue("model_part_name", settings["model_part_name"])
            component_params.AddValue("value", settings["value"][i])
            component_params.AddEmptyValue("variable_name").SetString(variable_name + "_" + suffix)
            if settings["Value_Table"][i].GetInt() == 0:
                component_params.AddValue("constrained", settings["constrained"][i])
                self.components_process_list.append(AssignScalarVariableProcess(Model, component_params))
            else:
                component_params.AddEmptyValue("is_fixed").SetBool(settings["constrained"][i].GetBool())
                component_params.AddValue("table", settings["Value_Table"][i])
                self.components_process_list.append(KratosPoro.ApplyComponentTableProcess(model_part, component_params))

    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
