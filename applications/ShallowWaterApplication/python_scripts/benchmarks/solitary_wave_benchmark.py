# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess
from KratosMultiphysics.ShallowWaterApplication.utilities.wave_factory import SolitaryWaveFactory

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SolitaryWaveBenchmark(model, settings["Parameters"])

class SolitaryWaveBenchmark(BaseBenchmarkProcess):
    """Solitary wave benchmark.

    Propagation of a solitary wave along the x axis.
    """

    def __init__(self, model, settings ):
        """Constructor of the benchmark.

        The base class validates the settings and set the model_part, variables and benchmark_settings
        """

        super().__init__(model, settings)

        benchmark_settings = settings["benchmark_settings"]
        self.boundary_model_part = model.GetModelPart(benchmark_settings["boundary_model_part_name"].GetString())
        process_info = self.boundary_model_part.ProcessInfo
        self.depth = benchmark_settings["depth"].GetDouble()
        self.wave = SolitaryWaveFactory(self.depth, benchmark_settings, process_info)
        self.x_shift = benchmark_settings["x_shift"].GetDouble()
        self.t_shift = benchmark_settings["t_shift"].GetDouble()


    @classmethod
    def _GetBenchmarkDefaultSettings(cls):
        return KM.Parameters("""
            {
                "boundary_model_part_name" : "",
                "wave_theory"              : "boussinesq",
                "depth"                    : 1.0,
                "amplitude"                : 1.0,
                "x_shift"                  : 0.0,
                "t_shift"                  : 0.0
            }
            """
            )


    def _Topography(self, coordinates):
        return -self.depth


    def _FreeSurfaceElevation(self, coordinates, time):
        x = coordinates.X
        return self.wave.eta(x - self.x_shift, time - self.t_shift)


    def _Height(self, coordinates, time):
        return self._FreeSurfaceElevation(coordinates, time) - self._Topography(coordinates)


    def _Velocity(self, coordinates, time):
        x = coordinates.X
        u_x = self.wave.u(x - self.x_shift, time - self.t_shift)
        return [u_x, 0.0, 0.0]


    def ExecuteInitializeSolutionStep(self):
        time = self.boundary_model_part.ProcessInfo[KM.TIME]
        for node in self.boundary_model_part.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, self._Velocity(node, time))
        KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, True, self.boundary_model_part.Nodes)


    def ExecuteFinalizeSolutionStep(self):
        KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, False, self.boundary_model_part.Nodes)
