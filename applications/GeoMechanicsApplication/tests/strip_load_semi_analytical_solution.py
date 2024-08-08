from typing import Callable

import math
import cmath


class StripLoad:
    """
    Class for solving the strip load problem semi-analytically

    A line load is applied to the surface of an elastic half-space dynamically. The solution is based on the following
    publication: An Introduction to Soil Dynamics , Verruijt A., 2009, Delft University of Technology, Chapter 12.2

    Attributes:
        - young (float): Young's modulus [Pa]
        - poisson (float): Poisson's ratio [-]
        - density (float): density [kg/m^3]
        - load (float): load value [N/m^2]
        - integral_stepsize (float): dimensionless time step size for solving integrals (default = 0.001) [-]
        - shear_modulus (float): shear modulus [Pa]
        - p_wave_modulus (float): P-wave modulus [Pa]
        - cs (float): shear wave velocity [m/s]
        - cp (float): compression wave velocity [m/s]
        - eta (float): ratio of shear wave to compression wave velocity [-]
        - epsilon (float): epsilon to avoid division by zero [-]

    """

    def __init__(self,
                 youngs_modulus: float,
                 poisson_ratio: float,
                 density: float,
                 load: float,
                 integral_stepsize: float = 0.001):
        """
        Constructor for the strip load class

        Args:
            - youngs_modulus (float): Young's modulus [Pa]
            - poisson_ratio (float): Poisson's ratio [-]
            - density (float): density [kg/m^3]
            - load (float): load value [N/m^2]
            - integral_stepsize (float): dimensionless time step size for solving integrals (default = 0.001) [-]
        """

        self.young = youngs_modulus
        self.poisson = poisson_ratio
        self.density = density
        self.load = load

        # dimensionless time step size for solving integrals
        self.integral_stepsize = integral_stepsize

        # calculate derived elastic properties
        self.shear_modulus = self.young / (2 * (1 + self.poisson))
        self.p_wave_modulus = self.young * (1 - self.poisson) / ((1 + self.poisson) * (1 - 2 * self.poisson))
        self.cs = math.sqrt(self.shear_modulus / self.density)  # shear wave velocity
        self.cp = math.sqrt(self.p_wave_modulus / self.density)  # compression wave velocity

        self.eta = self.cs / self.cp  # ratio of shear wave to compression wave velocity

        # epsilon to avoid division by zero
        self.epsilon = self.integral_stepsize**2

    @staticmethod
    def __calculate_radius(x_norm: float, z_norm: float) -> float:
        """
        Calculate the radius

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]

        Returns:
            - float: dimensionless radius [-]

        """

        return math.sqrt(x_norm**2 + z_norm**2)

    def __simpsons_rule(self, lower_limit: float, upper_limit: float, function: Callable, x_norm: float, z_norm: float,
                        radius: float, k0: float) -> float:
        """
        Calculate the integral of a function using Simpson's rule

        Args:
            - lower_limit (float): lower limit of the integral [-]
            - upper_limit (float): upper limit of the integral [-]
            - function (Callable): function to integrate [-]
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - radius (float): dimensionless radius [-]
            - k0 (float): initial imaginary factor for the integral [-]

        Returns:
            - float: integral of the function [-]

        """

        integral = 0
        # set dimensionless integration time parameter to lower limit
        kappa = lower_limit

        # solve integral with Simpson's rule
        h_1 = 0
        while kappa ** 2 < upper_limit ** 2:
            kappa += self.integral_stepsize
            h_2 = function(x_norm, z_norm, kappa, radius, k0)

            kappa += self.integral_stepsize
            h_3 = function(x_norm, z_norm, kappa, radius, k0)

            integral += (h_1 + 4 * h_2 + h_3) * self.integral_stepsize / 3
            h_1 = h_3

        return integral

    def __k1(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Imaginary part of the imaginary factor for first integral

        Args:
            - x_norm (float): dimensionless x coordinate normalized with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalized with the line load length [-]
            - kappa (float): dimensionless integration time parameter [-]

        Returns:
            - float: imaginary factor for the first integral [-]

        """

        radius = self.__calculate_radius(x_norm, z_norm)
        kappa_r = kappa**2 - self.eta**2 * radius**2

        if kappa_r <= 0:
            k1 = 0
        else:

            # dimensionless complex variables following from Laplace and Fourier transforms
            a = complex(kappa * math.sqrt(kappa_r), self.eta**2 * x_norm * z_norm)
            beta_p = complex(kappa * x_norm / radius**2, (z_norm / radius**2) * math.sqrt(kappa_r))
            b1 = 1 - 2 * beta_p**2
            gp = cmath.sqrt(self.eta**2 - beta_p**2)
            gs = cmath.sqrt(1 - beta_p**2)

            numerator = a * b1**2
            denominator = b1**2 + 4 * beta_p**2 * gp * gs

            # only get the imaginary part
            k1 = (numerator / denominator).imag / math.pi

        return k1

    def __h1(self, x_norm: float, z_norm: float, kappa: float, radius: float, k1_0: float) -> float:
        """
        Integrand of the first integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless integration time parameter [-]
            - radius (float): dimensionless radius [-]
            - k1_0 (float): initial imaginary factor for the first integral [-]

        Returns:
            - float: integrand of the first integral [-]

        """
        pz = kappa ** 2 - self.eta ** 2 * z_norm ** 2
        pr = kappa ** 2 - self.eta ** 2 * radius ** 2

        return (self.__k1(x_norm, z_norm, kappa) - k1_0) / (pz * math.sqrt(pr))

    def __f1(self, x_norm: float, z_norm: float, tau: float) -> float:
        """
        Calculates the first integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - tau (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: first integral solution [-]

        """

        # avoid division by zero
        if abs(x_norm) < self.epsilon:
            x_norm = self.epsilon if x_norm >= 0 else -self.epsilon

        radius = self.__calculate_radius(x_norm, z_norm)

        # dimensionless time parameter
        lower_limit_time = self.eta * radius + self.epsilon

        k1_ini = self.__k1(x_norm, z_norm, lower_limit_time)

        abs_x_norm = abs(x_norm)

        # solve first integral with simpson's rule
        if tau**2 <= lower_limit_time**2:
            first_integral = 0
        else:
            first_integral = ((k1_ini / (self.eta**2 * z_norm * abs_x_norm)) * math.atan(
                (z_norm * math.sqrt(tau**2 - self.eta**2 * radius**2)) / (tau * abs_x_norm)))

            # solve first integral with Simpson's rule
            first_integral += self.__simpsons_rule(lower_limit_time, tau, self.__h1, x_norm, z_norm, radius, k1_ini)

        if (tau > self.eta * z_norm) and (x_norm < 0):
            first_integral += 1

        return first_integral

    def __k2(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Real part of the imaginary factor for the second integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless integration time parameter [-]

        Returns:
            - float: factor for the second integral [-]

        """

        radius = self.__calculate_radius(x_norm, z_norm)
        kappa_r = kappa**2 - radius**2

        if kappa_r <= 0:
            k2 = 0
        else:
            # dimensionless complex variables following from Laplace and Fourier transforms
            beta_s = complex(kappa * x_norm / radius**2, (z_norm / radius**2) * math.sqrt(kappa_r))
            b1 = 1 - 2 * beta_s**2
            gp = cmath.sqrt(self.eta**2 - beta_s**2)
            gs = cmath.sqrt(1 - beta_s**2)

            numerator = 4 * beta_s * (1 - beta_s**2) * gp
            denominator = b1**2 + 4 * beta_s**2 * gp * gs

            k2 = (numerator / denominator).real / math.pi

        return k2

    def __h2(self, x_norm: float, z_norm: float, kappa: float, radius: float, k2_0: float) -> float:
        """
        Integrand of the second integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless integration time parameter [-]
            - radius (float): dimensionless radius [-]
            - k2_0 (float): initial imaginary factor for the second integral [-]

        Returns:
            - float: integrand of the second integral [-]

        """

        pr = kappa**2 - radius**2
        return (self.__k2(x_norm, z_norm, kappa) - k2_0) / math.sqrt(pr)

    def __f2(self, x_norm: float, z_norm: float, tau: float) -> float:
        """
        Calculates the second integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - tau (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: second integral solution [-]

        """

        if abs(x_norm) < self.epsilon:
            x_norm = self.epsilon if x_norm >= 0 else -self.epsilon

        radius = self.__calculate_radius(x_norm, z_norm)

        tau_r = tau / radius

        lower_limit_time = radius + self.epsilon

        k2 = self.__k2(x_norm, z_norm, lower_limit_time)

        if tau**2 <= lower_limit_time**2:
            second_integral = 0
        else:
            second_integral = k2 * math.log(tau_r + math.sqrt(tau_r**2 - 1))

            # solve second integral with simpson's rule
            second_integral += self.__simpsons_rule(lower_limit_time, tau, self.__h2, x_norm, z_norm, radius, k2)

        return second_integral

    def __h3(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Integrand of the third integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless integration time parameter [-]

        Returns:
            - float: integrand of the third integral [-]

        """

        radius = self.__calculate_radius(x_norm, z_norm)

        tau_q = self.eta * abs(x_norm) + z_norm * math.sqrt(1 - self.eta**2)

        if (kappa**2 >= radius**2) or (kappa <= tau_q) or (x_norm**2 <= self.eta**2 * radius**2):
            h3 = 0
        else:
            beta_q = (x_norm * kappa - z_norm * math.sqrt(radius**2 - kappa**2)) / radius**2

            b2 = (1 - 2 * beta_q**2)**2
            numerator = 4 * beta_q * (1 - beta_q**2) * b2 * math.sqrt(beta_q**2 - self.eta**2)

            denominator = b2**2 + 16 * beta_q**4 * (1 - beta_q**2) * (beta_q**2 - self.eta**2)

            h3 = -numerator / (math.pi * denominator * math.sqrt(radius**2 - kappa**2))

        return h3

    def __f3(self, x_norm: float, z_norm: float, tau: float) -> float:
        """
        Calculates the third integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - tau (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: third integral solution [-]

        """

        radius = self.__calculate_radius(x_norm, z_norm)
        tau_s = radius
        tau_q = self.eta * abs(x_norm) + z_norm * math.sqrt(1 - self.eta**2)

        if (tau <= tau_q) or (x_norm**2 <= self.eta**2 * radius**2):
            third_integral = 0
        else:

            # tau_s is the minimum of tau_s and tau
            tau_s = min(tau, tau_s)

            # determine integral stepsize as a function of tau_s and tau_q
            n_steps = 1000
            stepsize = (tau_s - tau_q) / n_steps

            third_integral = 0
            lower_limit_time = tau_q + stepsize / 2

            # set dimensionless integration time parameter to lower limit
            kappa = lower_limit_time

            # solve third integral
            while kappa < tau_s:
                third_integral += self.__h3(x_norm, z_norm, kappa) * stepsize
                kappa += stepsize

        return third_integral

    def calculate_normalised_vertical_stress(self, x_norm: float, z_norm: float, tau: float):
        """
        Calculate the normalised vertical stress

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - tau (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: normalised vertical stress [-]

        """
        normalised_vertical_stress = self.__f1(x_norm + 1, z_norm, tau) - self.__f1(x_norm - 1, z_norm, tau) + \
            self.__f2(x_norm + 1, z_norm, tau) - self.__f2(x_norm - 1, z_norm, tau) + \
            self.__f3(x_norm + 1, z_norm, tau) - self.__f3(x_norm - 1, z_norm, tau)
        return normalised_vertical_stress

    def calculate_vertical_stress(self, x: float, z: float, t: float, load_length: float, load_value: float) -> float:
        """
        Calculate the vertical stress

        Args:
            - x (float): x coordinate [m]
            - z (float): depth from load [m]
            - t (float): time [s]
            - load_length (float): length of the line load [m]
            - load_value (float): value of the line load [N/m^2]

        Returns:
            - float: vertical stress [Pa]

        """
        x_norm = x / load_length
        z_norm = z / load_length
        tau = self.cs * t / load_length

        normalised_vertical_stress = self.calculate_normalised_vertical_stress(x_norm, z_norm, tau)
        vertical_stress = normalised_vertical_stress * load_value
        return vertical_stress
