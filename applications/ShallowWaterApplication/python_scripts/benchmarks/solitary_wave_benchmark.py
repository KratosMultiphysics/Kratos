import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
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

    def ExecuteInitialize(self):
        # Construction of the wave settings
        benchmark_settings = self.settings["benchmark_settings"]
        self.boundary_model_part = self.model.GetModelPart(benchmark_settings["boundary_model_part_name"].GetString())
        self.wave = SolitaryWaveFactory(self.boundary_model_part, benchmark_settings["wave_specifications"])
        self.x_shift = benchmark_settings["wave_specifications"]["x_shift"].GetDouble()
        self.t_shift = benchmark_settings["wave_specifications"]["t_shift"].GetDouble()
        self.get_depth_from_model_part = benchmark_settings["wave_specifications"]["get_depth_from_model_part"].GetBool()

        # Here the base class sets the topography and initial conditions
        super().ExecuteInitialize()


    @classmethod
    def _GetBenchmarkDefaultSettings(cls):
        return KM.Parameters("""{
            "boundary_model_part_name"  : "",
            "wave_specifications"       : {
                "wave_theory"               : "boussinesq",
                "amplitude"                 : 1.0,
                "depth"                     : 1.0,
                "get_depth_from_model_part" : false,
                "x_shift"                   : 0.0,
                "t_shift"                   : 0.0
            }
        }""")


    def _Topography(self, coordinates):
        if self.get_depth_from_model_part:
            return 0
        else:
            return -self.wave.depth


    def _FreeSurfaceElevation(self, coordinates, time):
        x = coordinates.X
        return self.wave.eta(x - self.x_shift, time - self.t_shift)


    def _Height(self, coordinates, time):
        return self._FreeSurfaceElevation(coordinates, time) + self.wave.depth


    def _Velocity(self, coordinates, time):
        x = coordinates.X
        u_x = self.wave.u(x - self.x_shift, time - self.t_shift)
        return [u_x, 0.0, 0.0]


    def ExecuteBeforeSolutionLoop(self):
        SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.model_part)
        SW.ShallowWaterUtilities().SetMinimumValue(self.model_part, SW.HEIGHT, 0.0)
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.model_part)


    def ExecuteInitializeSolutionStep(self):
        time = self.boundary_model_part.ProcessInfo[KM.TIME]
        for node in self.boundary_model_part.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, self._Velocity(node, time))
            node.SetSolutionStepValue(SW.HEIGHT, self._Height(node, time))
        KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, True, self.boundary_model_part.Nodes)
        KM.VariableUtils().ApplyFixity(SW.HEIGHT, True, self.boundary_model_part.Nodes)


    def ExecuteFinalizeSolutionStep(self):
        KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, True, self.boundary_model_part.Nodes)
        KM.VariableUtils().ApplyFixity(SW.HEIGHT, False, self.boundary_model_part.Nodes)
