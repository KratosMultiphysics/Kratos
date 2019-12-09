import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveGeneratorProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class WaveGeneratorProcess(KM.Process):

    __formulation = {
        # Json input
        "reduced_variables" : SW.Formulation.REDUCED_VARIABLES,
        "conserved_variables" : SW.Formulation.CONSERVED_VARIABLES
    }

    __variables = {
        # Json input
        "free_surface" : SW.Variables.FREE_SURFACE_VARIABLE,
        "velocity" : SW.Variables.VELOCITY_VARIABLE,
        "free_surface_and_velocity" : SW.Variables.FREE_SURFACE_AND_VELOCITY
    }

    def __init__(self, model, settings ):
        KM.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "model_part_name"   : "model_part",
            "formulation"       : "reduced_variables",
            "variables"         : "free_surface",
            "interval"          : [0.0, 1e30],
            "wave_length"       : 10.0,
            "wave_period"       : 10.0,
            "wave_height"       : 1.0
        }
        """)

        # Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.interval = KM.IntervalUtility(settings)
        self.formulation = self.__formulation[settings["formulation"].GetString()]
        self.variables = self.__variables[settings["variables"].GetString()]
        self.fix_dofs = True

        # Definition of pi number
        import math

        # Wave parameters
        wave_height = settings["wave_height"].GetDouble()
        wave_period = settings["wave_period"].GetDouble()
        wave_length = settings["wave_length"].GetDouble()

        # Creation of the parameters for the c++ process
        free_surface_parameters = KM.Parameters("""{}""")
        free_surface_parameters.AddEmptyValue("amplitude").SetDouble(0.5 * wave_height)
        free_surface_parameters.AddEmptyValue("period").SetDouble(wave_period)
        free_surface_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)
        free_surface_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)

        velocity_parameters = KM.Parameters("""{}""")
        velocity_parameters.AddEmptyValue("amplitude").SetDouble(math.pi * wave_height / wave_period)
        velocity_parameters.AddEmptyValue("period").SetDouble(wave_period)
        velocity_parameters.AddEmptyValue("phase_shift").SetDouble(wave_period / 4)
        velocity_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)

        if self.variables == SW.Variables.VELOCITY_VARIABLE:
            velocity_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)

        self.free_surface_process = SW.ApplySinusoidalFunctionToScalar(self.model_part, SW.FREE_SURFACE_ELEVATION, free_surface_parameters)
        self.velocity_process = SW.ApplySinusoidalFunctionToVector(self.model_part, KM.VELOCITY, velocity_parameters)

        KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KM.DOMAIN_SIZE])
        SW.ShallowWaterUtilities().NormalizeVector(self.model_part, KM.NORMAL)

    def ExecuteInitializeSolutionStep(self):
        if self._IsInInterval():
            set_free_surface = self.variables == SW.Variables.FREE_SURFACE_VARIABLE or self.variables == SW.Variables.FREE_SURFACE_AND_VELOCITY
            set_velocity = self.variables == SW.Variables.VELOCITY_VARIABLE or self.variables == SW.Variables.FREE_SURFACE_AND_VELOCITY

            # Set the free surface if needed
            if set_free_surface:
                self.free_surface_process.ExecuteInitializeSolutionStep()
            
            # Set the velocity if needed
            if set_velocity:
                self.velocity_process.ExecuteInitializeSolutionStep()
            
            # Compute the free surface
            SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.model_part)

            # Compute the momentum if needed
            if self.formulation == SW.Formulation.CONSERVED_VARIABLES and set_velocity:
                SW.ShallowWaterUtilities().ComputeMomentum(self.model_part)

            # Fix the free surface if needed
            if self.fix_dofs and self.variables == set_free_surface:
                KM.VariableUtils().ApplyFixity(SW.HEIGHT, True, self.model_part.Nodes)
            
            # Fix the velocity or the momentum if needed
            if self.fix_dofs and self.variables == set_velocity:
                if self.formulation == SW.Formulation.REDUCED_VARIABLES:
                    KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, True, self.model_part.Nodes)
                    KM.VariableUtils().ApplyFixity(KM.VELOCITY_Y, True, self.model_part.Nodes)
                if self.formulation == SW.Formulation.CONSERVED_VARIABLES:
                    KM.VariableUtils().ApplyFixity(KM.MOMENTUM_X, True, self.model_part.Nodes)
                    KM.VariableUtils().ApplyFixity(KM.MOMENTUM_Y, True, self.model_part.Nodes)

    def _IsInInterval(self):
        """ Returns if we are inside the time interval or not """
        current_time = self.model_part.ProcessInfo[KM.TIME]
        return self.interval.IsInInterval(current_time)
