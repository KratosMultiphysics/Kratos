import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM


def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(Model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    return AssignBodyAccelerationToMaterialPointProcess(Model, settings["Parameters"])


class AssignBodyAccelerationToMaterialPointProcess(KratosMultiphysics.Process):
    """Assigns a body acceleration field directly to material point elements."""

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "MP_VOLUME_ACCELERATION",
                "mesh_id"              : 0,
                "value"                : [],
                "modulus"              : 1.0,
                "component"            : [0.0, 0.0, 0.0],
                "constrained"          : true,
                "interval"             : [0.0, "End"],
                "local_axes"           : {},
                "set_initial_mp_acceleration" : false
            }
            """)

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = self.GetDefaultParameters()
        self._AllowStringOrNumberValues(settings, default_settings)

        self.interval = KratosMultiphysics.IntervalUtility(settings)
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.set_initial_mp_acceleration = settings["set_initial_mp_acceleration"].GetBool()

        if self.variable != KratosMPM.MP_VOLUME_ACCELERATION:
            msg = 'Parameter "variable_name" only accepts "MP_VOLUME_ACCELERATION".'
            raise Exception(msg)

        if settings["value"].size() != 0:
            if settings["value"].size() != 3:
                raise Exception('Parameter "value" must be an array of length 3.')
            self.value_functions = self._CreateComponentFunctions(settings["value"], settings["local_axes"])
        else:
            self.modulus_function = self._CreateFunction(settings["modulus"], settings["local_axes"])
            self.component_functions = self._CreateComponentFunctions(settings["component"], settings["local_axes"])
            self.value_functions = None

    @staticmethod
    def _AllowStringOrNumberValues(settings, default_settings):
        if settings.Has("value") and settings["value"].size() != 0:
            default_settings["value"].SetStringArray(["0.0", "0.0", "0.0"])
        if settings.Has("modulus") and settings["modulus"].IsString():
            default_settings["modulus"].SetString("0.0")
        if settings.Has("component"):
            default_settings["component"].SetStringArray(["0.0", "0.0", "0.0"])

    @staticmethod
    def _CreateFunction(value_parameter, local_axes):
        if value_parameter.IsNumber():
            function_expression = str(value_parameter.GetDouble())
        elif value_parameter.IsString():
            function_expression = value_parameter.GetString()
        else:
            raise Exception("Function values must be numbers or strings.")

        return KratosMultiphysics.GenericFunctionUtility(function_expression, local_axes)

    def _CreateComponentFunctions(self, value_parameters, local_axes):
        return [self._CreateFunction(value_parameters[i], local_axes) for i in range(3)]

    def ExecuteBeforeSolutionLoop(self):
        self._AssignBodyAcceleration(self.set_initial_mp_acceleration)

    def ExecuteInitializeSolutionStep(self):
        self._AssignBodyAcceleration(False)

    def _AssignBodyAcceleration(self, set_mp_acceleration):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if not self.interval.IsInInterval(current_time):
            return

        process_info = self.model_part.ProcessInfo
        for element in self.model_part.Elements:
            mp_coord = element.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, process_info)[0]
            value = self._EvaluateValue(mp_coord, current_time)

            element.SetValuesOnIntegrationPoints(self.variable, [value], process_info)
            if set_mp_acceleration:
                element.SetValuesOnIntegrationPoints(KratosMPM.MP_ACCELERATION, [value], process_info)

    def _EvaluateValue(self, mp_coord, current_time):
        value = KratosMultiphysics.Vector(3)

        if self.value_functions is not None:
            for i in range(3):
                value[i] = self._CallFunction(self.value_functions[i], mp_coord, current_time)
        else:
            modulus = self._CallFunction(self.modulus_function, mp_coord, current_time)
            for i in range(3):
                value[i] = modulus * self._CallFunction(self.component_functions[i], mp_coord, current_time)

        return value

    @staticmethod
    def _CallFunction(function, mp_coord, current_time):
        return function.CallFunction(mp_coord[0], mp_coord[1], mp_coord[2], current_time, 0.0, 0.0, 0.0)
