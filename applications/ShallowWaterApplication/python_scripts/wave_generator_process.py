import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.ShallowWaterApplication.utilities import wave_theory_utilities


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

    __wave_theory = {
        "boussinesq"     : wave_theory_utilities.BoussinesqTheory,
        "linear_theory"  : wave_theory_utilities.LinearTheory,
        "shallow_theory" : wave_theory_utilities.ShallowTheory
    }

    def GetDefaultParameters(self):
        """The default settings depend on the specifications of the wave.
        The user settings can be in the project parameters or in the process info."""

        default_parameters = KM.Parameters("""
        {
            "model_part_name"   : "model_part",
            "formulation"       : "primitive_variables",
            "wave_theory"       : "boussinesq",
            "interval"          : [0.0, 1e30]
        }
        """)
        if self.settings.Has("wavelength"):
            default_parameters.SetDouble("wavelength", 0.0)

        if self.settings.Has("period"):
            default_parameters.SetDouble("period", 0.0)
        
        if self.settings.Has("amplitude"):
            default_parameters.SetDouble("amplitude", 0.0)

        return default_parameters


    def __init__(self, model, settings ):
        """"""
        KM.Process.__init__(self)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # Get the custom settings
        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.interval = KM.IntervalUtility(self.settings)
        self.formulation = self.__formulation[self.settings["formulation"].GetString()]
        self.fix_dofs = True


    def ExecuteInitialize(self):
        """Initialize the wave theory and the periodic function."""

        # Setup the wave theory
        wave_theory_class = self.__wave_theory[self.settings["wave_theory"].GetString()]
        depth = self._CalculateMeanDepth()
        gravity = self.model_part.ProcessInfo[KM.GRAVITY_Z]
        wave_theory = wave_theory_class(depth, gravity)
        wave_theory.SetWaveSpecifications(self.settings, self.model_part.ProcessInfo)

        # Creation of the parameters for the c++ process
        velocity_parameters = KM.Parameters("""{}""")
        velocity_parameters.AddEmptyValue("amplitude").SetDouble(wave_theory.horizontal_velocity)
        velocity_parameters.AddEmptyValue("period").SetDouble(wave_theory.period)
        velocity_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)
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
        """Returns if we are inside the time interval or not."""
        current_time = self.model_part.ProcessInfo[KM.TIME]
        return self.interval.IsInInterval(current_time)


    def _CalculateMeanDepth(self):
        sum_depths = -KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.TOPOGRAPHY, self.model_part, 0)
        mean_depth = sum_depths / self.model_part.NumberOfNodes()
        return mean_depth
