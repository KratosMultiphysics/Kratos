import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

## This process sets the value several scalar variables corresponding to the thermal problem using the AssignScalarVariableProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeThemalParametersScalarValueProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeThemalParametersScalarValueProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.components_process_list = []

        # Parameter name and assigned variable, in the original order. A zero
        # value means "do not create a process for this variable".
        thermal_parameters = (
            ("ThermalDensity", "DENSITY"),
            ("Conductivity", "CONDUCTIVITY"),
            ("SpecificHeat", "SPECIFIC_HEAT"),
        )

        for parameter_name, variable_name in thermal_parameters:
            if settings[parameter_name].GetDouble() != 0:
                component_params = KratosMultiphysics.Parameters("{}")
                component_params.AddValue("model_part_name", settings["model_part_name"])
                component_params.AddValue("value", settings[parameter_name])
                component_params.AddEmptyValue("constrained").SetBool(False)
                component_params.AddEmptyValue("variable_name").SetString(variable_name)
                # Preserve the legacy apply-once behavior of ApplyConstantScalarValueProcess.
                component_params.AddValue("interval", KratosMultiphysics.Parameters("[0.0, 0.0]"))
                self.components_process_list.append(
                    AssignScalarVariableProcess(Model, component_params)
                )


    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
