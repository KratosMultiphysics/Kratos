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

        if settings["ThermalDensity"].GetDouble() != 0:

            thermal_density = KratosMultiphysics.Parameters("{}")
            thermal_density.AddValue("model_part_name", settings["model_part_name"])
            thermal_density.AddValue("value", settings["ThermalDensity"])
            thermal_density.AddEmptyValue("constrained").SetBool(False)
            thermal_density.AddEmptyValue("variable_name").SetString("DENSITY")

            self.components_process_list.append(AssignScalarVariableProcess(Model, thermal_density))

        if settings["Conductivity"].GetDouble() != 0:

            conductivity = KratosMultiphysics.Parameters("{}")
            conductivity.AddValue("model_part_name", settings["model_part_name"])
            conductivity.AddValue("value", settings["Conductivity"])
            conductivity.AddEmptyValue("constrained").SetBool(False)
            conductivity.AddEmptyValue("variable_name").SetString("CONDUCTIVITY")

            self.components_process_list.append(AssignScalarVariableProcess(Model, conductivity))

        if settings["SpecificHeat"].GetDouble() != 0:

            specific_heat = KratosMultiphysics.Parameters("{}")
            specific_heat.AddValue("model_part_name", settings["model_part_name"])
            specific_heat.AddValue("value", settings["SpecificHeat"])
            specific_heat.AddEmptyValue("constrained").SetBool(False)
            specific_heat.AddEmptyValue("variable_name").SetString("SPECIFIC_HEAT")

            self.components_process_list.append(AssignScalarVariableProcess(Model, specific_heat))


    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

