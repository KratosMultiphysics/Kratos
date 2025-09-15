import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.ShallowWaterApplication.utilities.wave_theory_utilities as wave_theory
import KratosMultiphysics.ShallowWaterApplication.utilities.solitary_wave_utilities as solitary_wave


def WaveTheoryFactory(model_part, settings):
    default_settings = KM.Parameters("""{
        "wave_theory"               : "boussinesq",
        "period"                    : 0.0,
        "wavelength"                : 0.0,
        "amplitude"                 : 0.0,
        "depth"                     : 0.0,
        "get_depth_from_model_part" : true,
        "x_shift"                   : 0.0,
        "t_shift"                   : 0.0
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    wave_modules = {
        "boussinesq"      : wave_theory.BoussinesqTheory,
        "linear_theory"   : wave_theory.LinearTheory,
        "shallow_theory"  : wave_theory.ShallowTheory,
    }
    wave_module = wave_modules[settings["wave_theory"].GetString()]
    gravity = model_part.ProcessInfo[KM.GRAVITY_Z]
    depth = _GetDepth(settings, model_part)
    period = settings["period"].GetDouble()
    wavelength = settings["wavelength"].GetDouble()
    amplitude = settings["amplitude"].GetDouble()

    return wave_module(depth, gravity, period=period, wavelength=wavelength, amplitude=amplitude)


def SolitaryWaveFactory(model_part, settings):
    default_settings = KM.Parameters("""{
        "wave_theory"               : "boussinesq",
        "amplitude"                 : 0.0,
        "depth"                     : 0.0,
        "get_depth_from_model_part" : true,
        "x_shift"                   : 0.0,
        "t_shift"                   : 0.0
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    wave_modules = {
        "goring"      : solitary_wave.GoringSolution,
        "rayleigh"    : solitary_wave.RayleighSolution,
        "boussinesq"  : solitary_wave.BoussinesqSolution,
    }
    wave_module = wave_modules[settings["wave_theory"].GetString()]
    gravity = model_part.ProcessInfo[KM.GRAVITY_Z]
    depth = _GetDepth(settings, model_part)
    amplitude = settings["amplitude"].GetDouble()

    return wave_module(depth, gravity, amplitude=amplitude)


def _GetDepth(settings, model_part):
    if settings["get_depth_from_model_part"].GetBool():
        sum_depths = -KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.TOPOGRAPHY, model_part, 0)
        mean_depth = sum_depths / model_part.NumberOfNodes()
        return mean_depth
    else:
        return settings["depth"].GetDouble()
