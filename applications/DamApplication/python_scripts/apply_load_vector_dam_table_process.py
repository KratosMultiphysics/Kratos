import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLoadVectorDamTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyLoadVectorDamTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        # The Dam wrapper deliberately ignores "interval": the component
        # processes keep their standard default interval (the table process does
        # not use it either).
        self.components_process_list = []

        self.factor = settings["modulus"].GetDouble()
        self.direction = [settings["direction"][i].GetDouble() for i in range(3)]
        self.value = [component * self.factor for component in self.direction]

        for suffix, value in zip("XYZ", self.value):
            if abs(value) > 1.0e-15:
                component_params = KratosMultiphysics.Parameters("{}")
                component_params.AddValue("model_part_name", settings["model_part_name"])
                component_params.AddEmptyValue("value").SetDouble(value)
                component_params.AddEmptyValue("variable_name").SetString(variable_name + "_" + suffix)
                if settings["table"].GetInt() == 0:
                    component_params.AddEmptyValue("constrained").SetBool(False)
                    self.components_process_list.append(AssignScalarVariableProcess(Model, component_params))
                else:
                    component_params.AddEmptyValue("is_fixed").SetBool(False)
                    component_params.AddValue("table", settings["table"])
                    self.components_process_list.append(KratosPoro.ApplyComponentTableProcess(model_part, component_params))

    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
