import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.ShallowWaterApplication.utilities.wave_theory_utilities as wave_theory
import KratosMultiphysics.ShallowWaterApplication.utilities.solitary_wave_utilities as solitary_wave
import numbers


def WaveTheoryFactory(depth, settings, process_info):
    if not settings.Has("wave_theory"):
        raise Exception("Missing 'wave_theory' field in the project parameters.")

    wave_modules = {
        "boussinesq"      : wave_theory.BoussinesqTheory,
        "linear_theory"   : wave_theory.LinearTheory,
        "shallow_theory"  : wave_theory.ShallowTheory,
    }
    wave_module = wave_modules[settings["wave_theory"].GetString()]
    gravity = process_info.GetValue(KM.GRAVITY_Z)
    depth = _CheckDepth(depth)
    period = _CheckAndGetIfAvailable(settings, process_info, "period", SW.PERIOD)
    wavelength = _CheckAndGetIfAvailable(settings, process_info, "wavelength", SW.WAVELENGTH)
    amplitude = _CheckAndGetIfAvailable(settings, process_info, "amplitude", SW.AMPLITUDE)

    return wave_module(depth, gravity, period=period, wavelength=wavelength, amplitude=amplitude)


def SolitaryWaveFactory(depth, settings, process_info):
    if not settings.Has("wave_theory"):
        raise Exception("Missing 'wave_theory' field in the project parameters.")

    wave_modules = {
        "goring"      : solitary_wave.GoringSolution,
        "rayleigh"    : solitary_wave.RayleighSolution,
        "boussinesq"  : solitary_wave.BoussinesqSolution,
    }
    wave_module = wave_modules[settings["wave_theory"].GetString()]
    gravity = process_info.GetValue(KM.GRAVITY_Z)
    depth = _CheckDepth(depth)
    amplitude = _CheckAndGetIfAvailable(settings, process_info, "amplitude", SW.AMPLITUDE)

    return wave_module(depth, gravity, amplitude=amplitude)


def _CheckAndGetIfAvailable(parameters, process_info, name, variable):
    if parameters.Has(name):
        return parameters[name].GetDouble()
    elif process_info.Has(variable):
        return process_info.GetValue(variable)
    else:
        return 0


def _CheckDepth(depth):
    if isinstance(depth, KM.ModelPart):
        sum_depths = -KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.TOPOGRAPHY, depth, 0)
        mean_depth = sum_depths / depth.NumberOfNodes()
        return mean_depth
    elif isinstance(depth, numbers.Number):
        return depth
    else:
        raise Exception("'depth' argument must be a model part or a numeric value.")
