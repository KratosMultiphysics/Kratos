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
        # Json input
        "primitive_variables" : SW.Formulation.PrimitiveVariables,
        "conserved_variables" : SW.Formulation.ConservativeVariables
    }

    __variables = {
        # Json input
        "free_surface" : SW.Variables.FreeSurfaceVariable,
        "velocity" : SW.Variables.VelocityVariable,
        "free_surface_and_velocity" : SW.Variables.FreeSurfaceAndVelocity
    }

    def __init__(self, model, settings ):
        KM.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "model_part_name"   : "model_part",
            "formulation"       : "primitive_variables",
            "variables"         : "velocity",
            "interval"          : [0.0, 1e30],
            "wave_period"       : 10.0,
            "wave_height"       : 1.0,
            "smooth_time"       : 10.0
        }
        """)

        # Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.interval = KM.IntervalUtility(settings)
        self.formulation = self.__formulation[settings["formulation"].GetString()]
        self.variables = self.__variables[settings["variables"].GetString()]
        self.smooth_time = settings["smooth_time"].GetDouble()
        self.fix_dofs = True

        # Wave parameters
        self.wave_height = settings["wave_height"].GetDouble()
        self.wave_period = settings["wave_period"].GetDouble()

    def ExecuteInitialize(self):
        wave_amplitude = 0.5 * self.wave_height
        depth = -self.model_part.Nodes.__iter__().__next__().GetSolutionStepValue(SW.TOPOGRAPHY)
        gravity = self.model_part.ProcessInfo[KM.GRAVITY_Z]
        wave_velocity = sqrt(depth * gravity)
        velocity_amplitude = wave_amplitude * wave_velocity / depth

        # Creation of the parameters for the c++ process
        free_surface_parameters = KM.Parameters("""{}""")
        free_surface_parameters.AddEmptyValue("amplitude").SetDouble(wave_amplitude)
        free_surface_parameters.AddEmptyValue("period").SetDouble(self.wave_period)
        free_surface_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)
        free_surface_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)
        free_surface_parameters.AddEmptyValue("smooth_time").SetDouble(self.smooth_time)

        velocity_parameters = KM.Parameters("""{}""")
        velocity_parameters.AddEmptyValue("amplitude").SetDouble(velocity_amplitude)
        velocity_parameters.AddEmptyValue("period").SetDouble(self.wave_period)
        velocity_parameters.AddEmptyValue("phase_shift").SetDouble(self.wave_period / 4)
        velocity_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)
        free_surface_parameters.AddEmptyValue("smooth_time").SetDouble(self.smooth_time)

        if self.variables == SW.Variables.VelocityVariable:
            velocity_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)

        self.free_surface_process = SW.ApplySinusoidalFunctionToScalar(self.model_part, SW.FREE_SURFACE_ELEVATION, free_surface_parameters)
        self.velocity_process = SW.ApplySinusoidalFunctionToVector(self.model_part, KM.VELOCITY, velocity_parameters)

        KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KM.DOMAIN_SIZE])
        SW.ShallowWaterUtilities().NormalizeVector(self.model_part, KM.NORMAL)

    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        if self._IsInInterval():
            set_free_surface = self.variables == SW.Variables.FreeSurfaceVariable or self.variables == SW.Variables.FreeSurfaceAndVelocity
            set_velocity = self.variables == SW.Variables.VelocityVariable or self.variables == SW.Variables.FreeSurfaceAndVelocity

            # Set the free surface if needed
            if set_free_surface:
                self.free_surface_process.ExecuteInitializeSolutionStep()
                SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.model_part)
                if self.fix_dofs:
                    KM.VariableUtils().ApplyFixity(SW.HEIGHT, True, self.model_part.Nodes)
            
            # Set the velocity if needed
            if set_velocity:
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
