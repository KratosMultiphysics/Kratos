import KratosMultiphysics as Kratos


def Factory(settings, model):
    if (not isinstance(model, Kratos.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return InitializePreviousSolutionStepValuesProcess(model, settings["Parameters"])


class InitializePreviousSolutionStepValuesProcess(Kratos.Process):
    def __init__(self, model, settings):
        """This process sets nodal values of previous time step values

        Args:
            Model (Kratos.Model): [description]
            settings (Kratos.Parameters): [description]

        Raises:
            Exception: [description]
        """

        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""
        {
                "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
                "variable_name"  : "PLEASE_SPECIFY_VARIABLE_NAME",
                "value"          : 0.0,
                "interval"       : [0.0, 1e+30],
                "local_axes"     : {},
                "echo_level"     : 0
        }""")

        self.interval_utility = Kratos.IntervalUtility(settings)
        if (settings.Has("value") and settings["value"].IsString()):
            default_parameters["value"].SetString("0.0")
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.echo_level = settings["echo_level"].GetInt()
        self.variable = Kratos.KratosGlobals.GetVariable(
            settings["variable_name"].GetString())

        if not self.model_part.HasNodalSolutionStepVariable(self.variable):
            raise Exception("Solution step variable {:s} not found in solution step variables list in {:d}.".format(
                self.variable.Name(), self.model_part.FullName()))

        if (self.model_part.GetBufferSize() < 2):
            raise Exception(
                "Buffer size should be greater than 2. [ Buffer size = {:d}, model_part_name = {:s} ].".format(
                    self.model_part.GetBufferSize(), self.model_part.FullName()))

        self.value_is_numeric = False
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        else:
            self.function_string = settings["value"].GetString()
            self.aux_function = Kratos.PythonGenericFunctionUtility(
                self.function_string, settings["local_axes"])

            if self.aux_function.DependsOnSpace():
                self.cpp_apply_function_utility = Kratos.ApplyFunctionToNodesUtility(
                    self.model_part.Nodes, self.aux_function)

        self.variable_utils = Kratos.VariableUtils()

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Created previous solution step initialization process for {:s} in {:s}".format(
                self.variable.Name(), self.model_part.FullName()))

    def Execute(self):
        old_time = self.model_part.ProcessInfo.GetPreviousSolutionStepInfo()[
            Kratos.TIME]

        if self.interval_utility.IsInInterval(old_time):
            if self.value_is_numeric:
                self.variable_utils.SetVariable(
                    self.variable, self.value, self.model_part.Nodes, 1)
            else:
                if not self.aux_function.DependsOnSpace():  # depends on time only
                    self.value = self.aux_function.CallFunction(
                        0.0, 0.0, 0.0, old_time, 0.0, 0.0, 0.0)
                    self.variable_utils.SetVariable(
                        self.variable, self.value, self.model_part.Nodes, 1)
                # most general case - space varying function (possibly also time varying)
                else:
                    # needs a way with the new master to set previous time steps values
                    self.cpp_apply_function_utility.ApplyFunction(
                        self.variable, old_time, 1)

            if self.echo_level > 1:
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Initialized previous solution step {:s} in {:s}".format(
                    self.variable.Name(), self.model_part.FullName()))
