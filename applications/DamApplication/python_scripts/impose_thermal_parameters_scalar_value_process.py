from KratosMultiphysics import *

## This proces sets the value several scalar variables corresponding to the thermal problem using the ApplyConstantScalarValueProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeThemalParametersScalarValueProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeThemalParametersScalarValueProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []

        if settings["ThermalDensity"].GetDouble() != 0:

            thermal_density = Parameters("{}")
            thermal_density.AddValue("model_part_name", settings["model_part_name"])
            thermal_density.AddEmptyValue("is_fixed").SetBool(True)
            thermal_density.AddValue("value", settings["ThermalDensity"])
            thermal_density.AddEmptyValue("variable_name").SetString("DENSITY")

            self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, thermal_density))

        if settings["Conductivity"].GetDouble() != 0:

            conductivity = Parameters("{}")
            conductivity.AddValue("model_part_name", settings["model_part_name"])
            conductivity.AddEmptyValue("is_fixed").SetBool(True)
            conductivity.AddValue("value", settings["Conductivity"])
            conductivity.AddEmptyValue("variable_name").SetString("CONDUCTIVITY")

            self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, conductivity))

        if settings["SpecificHeat"].GetDouble() != 0:

            specific_heat = Parameters("{}")
            specific_heat.AddValue("model_part_name", settings["model_part_name"])
            specific_heat.AddEmptyValue("is_fixed").SetBool(True)
            specific_heat.AddValue("value", settings["SpecificHeat"])
            specific_heat.AddEmptyValue("variable_name").SetString("SPECIFIC_HEAT")

            self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, specific_heat))


    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

