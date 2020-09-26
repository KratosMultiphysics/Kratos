# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return BaseBenchmarkProcess(model, settings["Parameters"])

class BaseBenchmarkProcess(KM.Process):
    """The base class for the benchmarks."""

    def __init__(self, model, settings ):
        """The constructor of the BaseBenchmarkProcess.

        It is intended to be called from the constructor of deriving classes.
        """
        super().__init__()

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
        default_settings["benchmark_settings"] = self._GetBenchmarkDefaultSettings()
        settings.RecursivelyValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

        self.variables = GenerateVariableListFromInput(settings["variables_list"])
        self.exact_variables = GenerateVariableListFromInput(settings["exact_variables_list"])
        self.error_variables = GenerateVariableListFromInput(settings["error_variables_list"])
        self.benchmark_settings = settings["benchmark_settings"]

    def ExecuteInitialize(self):
        """This method sets the topography and the initial conditions"""
        time = self.model_part.ProcessInfo[KM.TIME]
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(SW.TOPOGRAPHY, self._Topography(node))
            node.SetSolutionStepValue(SW.HEIGHT, self._Height(node, time))
            node.SetSolutionStepValue(KM.VELOCITY, self._Velocity(node, time))
            node.SetSolutionStepValue(KM.MOMENTUM, self._Momentum(node, time))
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.model_part)

    def ExecuteFinalizeSolutionStep(self):
        """This method computes the exact values of the benchmark and computes the error of the simulation"""
        time = self.model_part.ProcessInfo[KM.TIME]
        for node in self.model_part.Nodes:
            for (variable, exact_variable, error_variable) in zip(self.variables, self.exact_variables, self.error_variables):
                if variable == SW.HEIGHT:
                    exact_value = self._Height(node, time)
                elif variable == KM._VELOCITY:
                    exact_value = self._Velocity(node, time)
                elif variable == KM.MOMENTUM:
                    exact_value = self._Momentum(node, time)

                fem_value = node.GetSolutionStepValue(variable)

                node.SetValue(exact_variable, exact_value)
                node.SetValue(error_variable, fem_value - exact_value)

    def Check(self):
        """This method checks if the input values have physical sense."""

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

    @classmethod
    def _GetBenchmarkDefaultSettings(cls):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark settings")

    def _Topography(self, coordinates):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")

    def _Height(self, coordinates, time):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")

    def _Velocity(self, coordinates, time):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")

    def _Momentum(self, coordinates, time):
        return [self._Height(coordinates, time)*v for v in self._Velocity(coordinates, time)]
