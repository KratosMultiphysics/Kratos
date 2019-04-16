import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveGeneratorProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class WaveGeneratorProcess(KratosMultiphysics.Process):

    __formulation = {
        # Json input
        "reduced_variables" : Shallow.Formulation.REDUCED_VARIABLES,
        "conserved_variables" : Shallow.Formulation.CONSERVED_VARIABLES
    }

    __variables = {
        # Json input
        "free_surface" : Shallow.Variables.FREE_SURFACE_VARIABLE,
        "velocity" : Shallow.Variables.VELOCITY_VARIABLE,
        "free_surface_and_velocity" : Shallow.Variables.FREE_SURFACE_AND_VELOCITY
    }

    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
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
        self.interval = KratosMultiphysics.IntervalUtility(settings)
        self.formulation = self.__formulation[settings["formulation"].GetString()]
        self.variables = self.__variables[settings["variables"].GetString()]

        # Definition of pi number
        import math

        # Wave parameters
        wave_height = settings["wave_height"].GetDouble()
        wave_period = settings["wave_period"].GetDouble()
        wave_length = settings["wave_length"].GetDouble()

        # Creation of the parameters for the c++ process
        free_surface_parameters = KratosMultiphysics.Parameters("""{}""")
        free_surface_parameters.AddEmptyValue("amplitude").SetDouble(0.5 * wave_height)
        free_surface_parameters.AddEmptyValue("period").SetDouble(wave_period)
        free_surface_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)
        free_surface_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)
        
        velocity_parameters = KratosMultiphysics.Parameters("""{}""")
        velocity_parameters.AddEmptyValue("amplitude").SetDouble(math.pi * wave_height / wave_period)
        velocity_parameters.AddEmptyValue("period").SetDouble(wave_period)
        velocity_parameters.AddEmptyValue("phase_shift").SetDouble(wave_period / 4)
        velocity_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)

        if self.variables == Shallow.Variables.VELOCITY_VARIABLE:
            velocity_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)

        self.free_surface_process = Shallow.ApplySinusoidalFunctionToScalar(self.model_part, Shallow.FREE_SURFACE_ELEVATION, free_surface_parameters)
        self.velocity_process = Shallow.ApplySinusoidalFunctionToScalar(self.model_part, Shallow.FREE_SURFACE_ELEVATION, velocity_parameters)
        self.variables_utility = Shallow.ShallowWaterVariablesUtility(self.model_part)

    def ExecuteInitializeSolutionStep(self):
        if self._IsInInterval():
            # Set the free surface if needed
            if self.variables == Shallow.Variables.FREE_SURFACE_VARIABLE or self.variables == Shallow.Variables.FREE_SURFACE_AND_VELOCITY:
                self.free_surface_process.ExecuteInitializeSolutionStep()
            
            # Set the velocity if needed
            if self.variables == Shallow.Variables.VELOCITY_VARIABLE or self.variables == Shallow.Variables.FREE_SURFACE_AND_VELOCITY:
                self.velocity_process.ExecuteInitializeSolutionStep()
            
            # Compute the free surface
            self.variables_utility.ComputeHeightFromFreeSurface()

            # Compute the momentum if needed
            if self.formulation == Shallow.Formulation.CONSERVED_VARIABLES:
                self.variables_utility.ComputeMomentum()

            # Fix the free surface if needed
            if self.variables == Shallow.Variables.FREE_SURFACE_VARIABLE or self.variables == Shallow.Variables.FREE_SURFACE_AND_VELOCITY:
                KratosMultiphysics.VariableUtils().ApplyFixity(Shallow.HEIGHT, True, self.model_part.Nodes)
            
            # Fix the velocity or the momentum if needed
            if self.variables == Shallow.Variables.VELOCITY_VARIABLE or self.variables == Shallow.Variables.FREE_SURFACE_AND_VELOCITY:
                if self.formulation == Shallow.Formulation.REDUCED_VARIABLES:
                    KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY, True, self.model_part.Nodes)
                if self.formulation == Shallow.Formulation.CONSERVED_VARIABLEs:
                    KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MOMENTUM, True, self.model_part.Nodes)


    def _IsInInterval(self):
        """ Returns if we are inside the time interval or not """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        return self.interval.IsInInterval(current_time)
