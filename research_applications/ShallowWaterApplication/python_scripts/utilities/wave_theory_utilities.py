from math import pi, sqrt, tanh
from scipy.optimize import root


class WaveTheory:
    '''Base class for waves calculations.'''

    def __init__(self, depth, gravity=9.81, *, period=0, wavelength=0, amplitude=0):
        self.gravity = gravity
        self.depth = depth
        if self.depth <= 0:
            raise Exception('WaveTheory. The water depth must be greather than 0.')
        if self.gravity <= 0:
            raise Exception('WaveTheory. The gravity must be greather than 0.')
        if wavelength > 0 and period > 0:
            raise Exception('WaveTheory. Specify only the wavelength or the period.')
        if period > 0:
            self.SetPeriod(period)
        if wavelength > 0:
            self.SetWavelength(wavelength)
        if amplitude > 0:
            self.SetAmplitude(amplitude)

    def SetPeriod(self, period):
        self.period = period
        self.wavelength = self._CalculateWavelength(period)

    def SetWavelength(self, wavelength):
        self.wavelength = wavelength
        self.period = self._CalculatePeriod(wavelength)

    def SetAmplitude(self, amplitude):
        self.amplitude = amplitude

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
        raise Exception('WaveTheory base class: it is not possible to calculate the disperison relation.')

    def _HorizontalVelocity(self, amplitude, frequency, wavenumber):
        raise Exception('WaveTheory base class: it is not possible to calculate the horizontal velocity.')

    def _PhaseSpeed(self, wavenumber):
        raise Exception('WaveTheory base class: it is not possible to calculate the phase speed.')

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
        return sqrt(max(0, gh * (1 -(self.alpha + 1/3) * kh**2) / (1 -self.alpha * kh**2)))


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
