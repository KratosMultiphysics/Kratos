from numpy import pi, sqrt, tanh
from scipy.optimize import root
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW


class WaveTheory:
    '''Base class for waves calculations.'''

    def __init__(self, depth, gravity=9.81, *, period=0, wavelength=0, amplitude=0):
        self.gravity = gravity
        self.depth = depth
        if wavelength > 0 and amplitude > 0:
            raise Exception('WaveTheory. Specify only the wavelength or the period.')
        self.SetPeriod(period)
        self.SetWavelength(wavelength)
        self.SetAmplitude(amplitude)

    def SetPeriod(self, period):
        if period > 0:
            self.period = period
            self.wavelength = self._CalculateWavelength(period)

    def SetWavelength(self, wavelength):
        if wavelength > 0:
            self.wavelength = wavelength
            self.period = self._CalculatePeriod(wavelength)

    def SetAmplitude(self, amplitude):
        if amplitude > 0:
            self.amplitude = amplitude

    def SetWaveSpecifications(self, parameters = KM.Parameters(), process_info = KM.ProcessInfo()):
        _CheckAndSetWaveSpecifications(self, parameters, process_info)

    @property
    def horizontal_velocity(self):
        return self._HorizontalVelocity(self.amplitude, self.frequency, self.wavenumber)

    @property
    def phase_speed(self):
        return self._PhaseSpeed(self.wavenumber)

    @property
    def wavenumber(self):
        return 2 * pi / self.wavelength

    @property
    def frequency(self):
        return 2 * pi / self.period

    def _DispersionRelation(self, wavenumber):
        KM.Logger.PrintWarning('WaveTheory base class: it is not possible to calculate the disperison relation.')
        return 0

    def _HorizontalVelocity(self, amplitude, frequency, wavenumber):
        KM.Logger.PrintWarning('WaveTheory base class: it is not possible to calculate the horizontal velocity.')
        return 0

    def _PhaseSpeed(self, wavenumber):
        KM.Logger.PrintWarning('WaveTheory base class: it is not possible to calculate the phase speed.')
        return 0

    def _CalculateFrequency(self, wavenumber):
        return sqrt(self._DispersionRelation(wavenumber))

    def _CalculatePeriod(self, wavelength):
        wavenumber = 2 * pi / wavelength
        frequency = self._CalculateFrequency(wavenumber)
        return 2 * pi / frequency

    def _CalculateWavenumber(self, frequency):
        k0 = frequency / sqrt(self.gravity * self.depth)
        func = lambda k: self._CalculateFrequency(k) - frequency
        result = root(func, k0)
        if len(result.x) == 1:
            return result.x[0]
        else:
            return result.x

    def _CalculateWavelength(self, period):
        frequency = 2*pi / period
        wavenumber = self._CalculateWavenumber(frequency)
        return 2 * pi / wavenumber


class BoussinesqTheory(WaveTheory):
    '''Boussinesq theory for dispersive waves.'''

    beta = -0.531
    alpha = 0.5 * beta**2 + beta

    def _DispersionRelation(self, wavenumber):
        g = self.gravity
        kh = wavenumber * self.depth
        return g * kh * wavenumber * (1 -(self.alpha + 1/3) * kh**2) / (1 -self.alpha * kh**2)

    def _HorizontalVelocity(self, amplitude, frequency, wavenumber):
        kh = wavenumber * self.depth
        return frequency * amplitude / kh / (1 -(self.alpha + 1/3) * kh**2)

    def _PhaseSpeed(self, wavenumber):
        kh = wavenumber * self.depth
        gh = self.gravity * self.depth
        return sqrt(gh * (1 -(self.alpha + 1/3) * kh**2) / (1 -self.alpha * kh**2))


class LinearTheory(WaveTheory):
    '''Linear theory for intermediate water.'''

    def _DispersionRelation(self, wavenumber):
        g = self.gravity
        kh = wavenumber * self.depth
        return g * wavenumber * tanh(kh)

    def _HorizontalVelocity(self, amplitude, frequency, wavenumber):
        kh = wavenumber * self.depth
        return frequency * amplitude / kh # Note: 1/h * int{cosh(k*(z+h)) dz}_{-h}^{0} = sinh(kh) / kh

    def _PhaseSpeed(self, wavenumber):
        kh = wavenumber * self.depth
        gh = self.gravity * self.depth
        return sqrt(gh * tanh(kh) / kh)


class ShallowTheory(WaveTheory):
    '''Linear theory for shallow water.'''

    def _DispersionRelation(self, wavenumber):
        g = self.gravity
        kh = wavenumber * self.depth
        return g * wavenumber * kh

    def _HorizontalVelocity(self, amplitude, frequency, wavenumber):
        kh = wavenumber * self.depth
        return frequency * amplitude / kh

    def _PhaseSpeed(self, wavenumber):
        gh = self.gravity * self.depth
        return sqrt(gh)


def _CheckAndSetWaveSpecifications(wave_theory, parameters, process_info):
    # Check and get the wave specifications if provided
    period = _CheckAndGetIfAvailable(parameters, process_info, "period", SW.PERIOD)
    wavelength = _CheckAndGetIfAvailable(parameters, process_info, "wavelength", SW.WAVELENGTH)
    amplitude = _CheckAndGetIfAvailable(parameters, process_info, "amplitude", SW.AMPLITUDE)

    # Check if the wave specification is unique
    if period and wavelength:
        raise Exception("WaveGeneratorProcess. Provide the period or the wavelength. Both parameters are incompatible.")
    if not period and not wavelength:
        raise Exception("WaveGeneratorProcess. Please, specify the wavelength or the period in the project paramenters or hte process info.")
    if not amplitude:
        raise Exception("WaveGeneratorProcess. Please, specify the amplitude in the project parameters or the process info.")

    # Apply the user settings
    wave_theory.SetWavelength(wavelength)
    wave_theory.SetPeriod(period)
    wave_theory.SetAmplitude(amplitude)


def _CheckAndGetIfAvailable(parameters, process_info, name, variable):
    if parameters.Has(name):
        return parameters[name].GetDouble()
    elif process_info.Has(variable):
        return process_info.GetValue(variable)
    else:
        return 0
