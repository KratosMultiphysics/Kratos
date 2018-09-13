from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeFaceHeatFluxProcess(Model, settings["Parameters"])

class ImposeFaceHeatFluxProcess(Process):

    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        self.components_process_list = []

        ## This process assign and uniform heat flux
        if ("UniformFlux2D" in settings["model_part_name"].GetString()) or ("UniformFlux3D" in settings["model_part_name"].GetString()):
            if settings["table"].GetInt() == 0:
                t_uniform = Parameters("{}")
                t_uniform.AddValue("model_part_name",settings["model_part_name"])
                t_uniform.AddEmptyValue("is_fixed").SetBool(False)
                t_uniform.AddValue("variable_name",settings["variable_name"])
                t_uniform.AddValue("value",settings["value"])
                self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, t_uniform))
            else:
                self.components_process_list.append(DamFixTemperatureConditionProcess(model_part, settings))

        ## This process compute the heat flux according to q = h(t_ambient - t_current).
        ## Setting the extra values to 0.0 it is possible to use the same process.
        if ("TAmbientFlux2D" in settings["model_part_name"].GetString()) or ("TAmbientFlux3D" in settings["model_part_name"].GetString()):
            t_ambient = Parameters("{}")
            t_ambient.AddValue("model_part_name",settings["model_part_name"])
            t_ambient.AddValue("variable_name",settings["variable_name"])
            t_ambient.AddValue("h_0",settings["h_0"])
            t_ambient.AddValue("ambient_temperature",settings["ambient_temperature"])
            t_ambient.AddValue("table_ambient_temperature",settings["table_ambient_temperature"])
            t_ambient.AddEmptyValue("emisivity").SetDouble(0.0)
            t_ambient.AddEmptyValue("delta_R").SetDouble(0.0)
            t_ambient.AddEmptyValue("absorption_index").SetDouble(0.0)
            t_ambient.AddEmptyValue("total_insolation").SetDouble(0.0)
            self.components_process_list.append(DamTSolAirHeatFluxProcess(model_part, t_ambient))

        ## This process compute the heat flux according to q = h(t_sol_air - t_current)
        if ("TSolAirFluxCondition2D" in settings["model_part_name"].GetString()) or ("TSolAirFluxCondition3D" in settings["model_part_name"].GetString()):
            self.components_process_list.append(DamTSolAirHeatFluxProcess(model_part, settings))

    def ExecuteInitialize(self):
        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
