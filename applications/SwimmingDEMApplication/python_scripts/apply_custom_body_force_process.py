# Importing the Kratos Library
import KratosMultiphysics

# other imports
from multiple_points_output_process import MultiplePointsOutputProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomBodyForceProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "BODY_FORCE",
                "benchmark_name"       : "vortex",
                "benchmark_parameters" : {},
                "print_point_output"   : false,
                "output_parameters"    : {}
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())

        benchmark_module = __import__(settings["benchmark_name"].GetString(), fromlist=[None])
        self.benchmark = benchmark_module.CreateManufacturedSolution(settings["benchmark_parameters"])

        self.print_output = settings["print_point_output"].GetBool()
        if self.print_output:
            self.output_process = MultiplePointsOutputProcess(model, settings["output_parameters"])

    def ExecuteInitialize(self):
        if self.print_output:
            self.output_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        '''
        This is the place to set the initial conditions given by the benchmark
        '''
        pass

    def ExecuteInitializeSolutionStep(self):
        self._SetBodyForce()

    def ExecuteBeforeOutputStep(self):
        self._ComputeVelocityBenchmark()

    def ExecuteFinalize(self):
        if self.print_output:
            self._ComputeVelocityError()
            self._CopyVelocityAsNonHistorical()
            self.output_process.ExecuteFinalizeSolutionStep()
            self.output_process.ExecuteFinalize

    def _SetBodyForce(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            value = self.benchmark.BodyForce(node.X, node.Y, node.Z, current_time)
            node.SetSolutionStepValue(self.variable, value)

    def _ComputeVelocityError(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            fem_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            exact_vel = self.benchmark.Velocity(node.X, node.Y, node.Z, current_time)
            error = abs((fem_vel[0]**2 + fem_vel[1]**2 + fem_vel[2]**2)**0.5 - (exact_vel[0]**2 + exact_vel[1]**2 + exact_vel[2]**2)**0.5)
            node.SetValue(KratosMultiphysics.NODAL_ERROR, error)

    def _CopyVelocityAsNonHistorical(self):
        for node in self.model_part.Nodes:
            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetValue(KratosMultiphysics.VELOCITY, vel)

    def _ComputeVelocityBenchmark(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            vel = self.benchmark.Velocity(node.X, node.Y, node.Z, current_time)
            node.SetValue(KratosMultiphysics.Y, vel)
