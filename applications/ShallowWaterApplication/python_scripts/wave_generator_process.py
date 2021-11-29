import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from math import pi, sqrt

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveGeneratorProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class WaveGeneratorProcess(KM.Process):

    __formulation = {
        "primitive_variables" : SW.Formulation.PrimitiveVariables,
        "conserved_variables" : SW.Formulation.ConservativeVariables
    }

    def __init__(self, model, settings ):
        KM.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "model_part_name"   : "model_part",
            "formulation"       : "primitive_variables",
            "interval"          : [0.0, 1e30],
            "wave_length"       : 10.0,
            "wave_height"       : 1.0,
            "smooth_time"       : 10.0
        }
        """)

        # Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.interval = KM.IntervalUtility(settings)
        self.formulation = self.__formulation[settings["formulation"].GetString()]
        self.smooth_time = settings["smooth_time"].GetDouble()
        self.fix_dofs = True

        # Wave parameters
        wave_height = settings["wave_height"].GetDouble()
        wave_length = settings["wave_length"].GetDouble()
        self.depth = self._CalculateMeanDepth()

        self.wave_amplitude = 0.5 * wave_height
        self.wavenumber = 2 * pi / wave_length
        self.frequency = sqrt(self._DispersionRelation(self.wavenumber))


    def ExecuteInitialize(self):
        wave_velocity = self.frequency / self.wavenumber
        velocity_amplitude = self.wave_amplitude * wave_velocity / self.depth

        # Creation of the parameters for the c++ process
        velocity_parameters = KM.Parameters("""{}""")
        velocity_parameters.AddEmptyValue("amplitude").SetDouble(velocity_amplitude)
        velocity_parameters.AddEmptyValue("period").SetDouble(self.wave_period)
        velocity_parameters.AddEmptyValue("phase_shift").SetDouble(self.wave_period / 4)
        velocity_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)

        self.velocity_process = SW.ApplySinusoidalFunctionToVector(self.model_part, KM.VELOCITY, velocity_parameters)

        KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KM.DOMAIN_SIZE])
        SW.ShallowWaterUtilities().NormalizeVector(self.model_part, KM.NORMAL)


    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()


    def ExecuteInitializeSolutionStep(self):
        if self._IsInInterval():
            self.velocity_process.ExecuteInitializeSolutionStep()
            if self.formulation == SW.Formulation.ConservativeVariables:
                SW.ShallowWaterUtilities().ComputeMomentum(self.model_part)
            if self.fix_dofs:
                if self.formulation == SW.Formulation.PrimitiveVariables:
                    KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, True, self.model_part.Nodes)
                    KM.VariableUtils().ApplyFixity(KM.VELOCITY_Y, True, self.model_part.Nodes)
                if self.formulation == SW.Formulation.ConservativeVariables:
                    KM.VariableUtils().ApplyFixity(KM.MOMENTUM_X, True, self.model_part.Nodes)
                    KM.VariableUtils().ApplyFixity(KM.MOMENTUM_Y, True, self.model_part.Nodes)


    def _IsInInterval(self):
        """ Returns if we are inside the time interval or not """
        current_time = self.model_part.ProcessInfo[KM.TIME]
        return self.interval.IsInInterval(current_time)


    def _CalculateMeanDepth(self):
        sum_depths = -KM.VariableUtils().SumHistoricalVariable(SW.TOPOGRAPHY, self.model_part, 0)
        mean_depth = sum_depths / self.model_part.NumberOfNodes()
        return mean_depth


    def _DispersionRelation(self, wavenumber):
        g = self.model_part.ProcessInfo()[KM.GRAVITY_Z]
        kh = wavenumber * self.depth
        beta = -0.531
        alpha = 0.5 * beta**2 + beta
        return g * kh * wavenumber * (1 -(alpha + 1/3) * kh**2) / (1 -alpha * kh**2)
