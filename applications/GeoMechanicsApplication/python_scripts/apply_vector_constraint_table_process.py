import KratosMultiphysics as Core
import KratosMultiphysics.GeoMechanicsApplication as Geo


def Factory(settings, Model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyVectorConstraintTableProcess(Model, settings["Parameters"])


## All the python processes should be derived from "python_process"

class ApplyVectorConstraintTableProcess(Core.Process):
    def __init__(self, Model, settings):
        super().__init__()
        self.model_part = Model[settings["model_part_name"].GetString()]
        self._CreateParametersForComponents(settings)
        self._AddProcesses()


    def _CreateParametersForComponents(self, settings):
        self.x_params = self._CreateParametersForComponent(settings, "X")
        self.y_params = self._CreateParametersForComponent(settings, "Y")
        self.z_params = self._CreateParametersForComponent(settings, "Z")


    def _AddProcesses(self):
        self.components_process_list = []
        self._AddProcessUsing(self.x_params)
        self._AddProcessUsing(self.y_params)
        self._AddProcessUsing(self.z_params)


    def _CreateParametersForComponent(self, settings, component):
        index = self.ComponentToIndex(component)
        if not settings["active"][index].GetBool():
            return None

        parameters = Core.Parameters("{}")
        parameters.AddValue("model_part_name", settings["model_part_name"])
        if settings.Has("is_fixed"):
            parameters.AddValue("is_fixed", settings["is_fixed"][index])
        parameters.AddValue("value", settings["value"][index])
        variable_name = settings["variable_name"].GetString()
        parameters.AddEmptyValue("variable_name").SetString(f"{variable_name}_{component}")
        if settings["table"][index].GetInt() != 0:
            parameters.AddValue("table", settings["table"][index])

        return parameters


    def _AddProcessUsing(self, parameters):
        if not parameters:
            return

        if parameters.Has("table"):
            process = Geo.ApplyComponentTableProcess(self.model_part, parameters)
        else:
            process = Core.ApplyConstantScalarValueProcess(self.model_part, parameters)

        self.components_process_list.append(process)


    @classmethod
    def ComponentToIndex(cls, component):
        if component == "X":
            return 0
        if component == "Y":
            return 1
        if component == "Z":
            return 2

        raise LookupError(f"Unknown component '{component}'")


    def ExecuteInitialize(self):
        for component in self.components_process_list:
            component.ExecuteInitialize()


    def ExecuteInitializeSolutionStep(self):
        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
