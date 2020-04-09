# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return BaseBenchmarkProcess(model, settings["Parameters"])

class BaseBenchmarkProcess(KM.Process):

    def __init__(self, model, settings ):
        super(BaseBenchmarkProcess, self).__init__()

        default_settings = KM.Parameters("""
            {
                "model_part_name"      : "model_part",
                "variables_list"       : [],
                "exact_variables_list" : [],
                "error_variables_list" : [],
                "benchmark_settings"   : {}
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

        self.variables = GenerateVariableListFromInput(settings["variables_list"])
        self.exact_variables = GenerateVariableListFromInput(settings["exact_variables_list"])
        self.error_variables = GenerateVariableListFromInput(settings["error_variables_list"])
        self.benchmark_settings = settings["benchmark_settings"]

    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KM.TIME]

        for node in self.model_part.Nodes:
            for (variable, exact_variable, error_variable) in zip(self.variables, self.exact_variables, self.error_variables):
                if variable == SW.HEIGHT:
                    exact_value = self.Height(node, time)
                elif variable == KM.VELOCITY:
                    exact_value = self.Velocity(node, time)
                elif variable == KM.MOMENTUM:
                    exact_value = self.Momentum(node, time)

                fem_value = node.GetSolutionStepValue(variable)

                node.SetValue(exact_variable, exact_value)
                node.SetValue(error_variable, fem_value - exact_value)

    def Check(self):
        if len(self.variables) != len(self.exact_variables):
            raise Exception("The input variables list does not match the input exact variables list")

        if len(self.variables) != len(self.error_variables):
            raise Exception("The input variables list does not match the input error variables list")

        for (var, exact, error) in zip(self.variables, self.exact_variables, self.error_variables):
            if KM.KratosGlobals.GetVariableType(var.Name()) != KM.KratosGlobals.GetVariableType(exact.Name()):
                msg = var.Name() + " variable type does not match the " + exact.Name() + " variable type"
                raise Exception(msg)

            if KM.KratosGlobals.GetVariableType(var.Name()) != KM.KratosGlobals.GetVariableType(error.Name()):
                msg = var.Name() + " variable type does not match the " + error.Name() + " variable type"
                raise Exception(msg)

    def Height(self, coordinates, time):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")

    def Velocity(self, coordinates, time):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")

    def Momentum(self, coordinates, time):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")
