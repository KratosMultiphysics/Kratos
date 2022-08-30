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

        self.model = model
        self.settings = settings

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
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)

        self.model_part = self.model[self.settings["model_part_name"].GetString()]
        self.variables = GenerateVariableListFromInput(self.settings["variables_list"])
        self.exact_variables = GenerateVariableListFromInput(self.settings["exact_variables_list"])
        self.error_variables = GenerateVariableListFromInput(self.settings["error_variables_list"])


    def ExecuteInitialize(self):
        """Set the topography and the initial conditions."""
        KM.Timer.Start("Benchmark/Initial state")
        time = self.model_part.ProcessInfo[KM.TIME]
        for node in self.model_part.Nodes:
            if self._Topography(node):
                node.SetSolutionStepValue(SW.TOPOGRAPHY, self._Topography(node))
            node.SetSolutionStepValue(SW.HEIGHT, self._Height(node, time))
            node.SetSolutionStepValue(KM.VELOCITY, self._Velocity(node, time))
            node.SetSolutionStepValue(KM.MOMENTUM, self._Momentum(node, time))
            node.SetSolutionStepValue(SW.FREE_SURFACE_ELEVATION, self._FreeSurfaceElevation(node, time))
        KM.Timer.Stop("Benchmark/Initial state")


    def ExecuteBeforeOutputStep(self):
        """Compute the exact values of the benchmark and the error of the simulation."""
        KM.Timer.Start("Benchmark/Exact values")
        time = self.model_part.ProcessInfo[KM.TIME]
        for (variable, exact_variable, error_variable) in zip(self.variables, self.exact_variables, self.error_variables):

            if variable == SW.HEIGHT:
                exact_value_function = self._Height
            elif variable == KM.VELOCITY:
                exact_value_function = self._Velocity
            elif variable == KM.MOMENTUM:
                exact_value_function = self._Momentum
            elif variable == SW.FREE_SURFACE_ELEVATION:
                exact_value_function = self._FreeSurfaceElevation

            for node in self.model_part.Nodes:
                exact_value = exact_value_function(node, time)
                fem_value = node.GetSolutionStepValue(variable)

                node.SetValue(exact_variable, exact_value)
                node.SetValue(error_variable, fem_value - exact_value)
        KM.Timer.Stop("Benchmark/Exact values")


    def Check(self):
        """Check if the input values have physical sense."""

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
        return 0

    def _Height(self, coordinates, time):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")

    def _Velocity(self, coordinates, time):
        raise Exception("Calling the base class of the benchmark. Please, implement the custom benchmark")

    def _Momentum(self, coordinates, time):
        return [self._Height(coordinates, time)*v for v in self._Velocity(coordinates, time)]

    def _FreeSurfaceElevation(self, coordinates, time):
        return self._Topography(coordinates) + self._Height(coordinates, time)
