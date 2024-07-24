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
        - load (float): load value [N/m]
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
            - load (float): load value [N/m]
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

    def __k1(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Imaginary part of the imaginary factor for first integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: imaginary factor for the first integral [-]
        """

        radius = math.sqrt(x_norm**2 + z_norm**2)
        kappa_r = kappa**2 - self.eta**2 * radius**2

        if kappa_r <= 0:
            k1 = 0
        else:

            # dimensionless complex variables following from laplace and fourier transforms
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

    def __f1(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Calculates the first integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: first integral solution [-]

        """

        # avoid division by zero
        if abs(x_norm) < self.epsilon:
            x_norm = self.epsilon if x_norm >= 0 else -self.epsilon

        radius = math.sqrt(x_norm**2 + z_norm**2)

        # dimensionless time parameter
        time_parameter = self.eta * radius + self.epsilon

        k1_ini = self.__k1(x_norm, z_norm, time_parameter)

        x_new = max(math.sqrt(radius**2 - z_norm**2), self.epsilon)

        # solve first integral with simpson's rule
        if kappa**2 <= time_parameter**2:
            first_integral = 0
        else:
            first_integral = ((k1_ini / (self.eta**2 * z_norm * x_new)) * math.atan(
                (z_norm * math.sqrt(kappa**2 - self.eta**2 * radius**2)) / (kappa * x_new)))
            h1_1 = 0

            # solve first integral with simpson's rule
            while time_parameter**2 < kappa**2:
                time_parameter += self.integral_stepsize

                pz = time_parameter**2 - self.eta**2 * z_norm**2
                pr = time_parameter**2 - self.eta**2 * radius**2

                h1_2 = (self.__k1(x_norm, z_norm, time_parameter) - k1_ini) / (pz * math.sqrt(pr))

                time_parameter += self.integral_stepsize

                pz = time_parameter**2 - self.eta**2 * z_norm**2
                pr = time_parameter**2 - self.eta**2 * radius**2

                h1_3 = (self.__k1(x_norm, z_norm, time_parameter) - k1_ini) / (pz * math.sqrt(pr))

                first_integral += (h1_1 + 4 * h1_2 + h1_3) * self.integral_stepsize / 3
                h1_1 = h1_3

        if (kappa > self.eta * z_norm) and (x_norm < 0):
            first_integral += 1

        return first_integral

    def __k2(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Real part of the imaginary factor for the second integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]
        """

        radius = math.sqrt(x_norm**2 + z_norm**2)
        kappa_r = kappa**2 - radius**2

        if kappa_r <= 0:
            k2 = 0
        else:
            # dimensionless complex variables following from laplace and fourier transforms
            beta_s = complex(kappa * x_norm / radius**2, (z_norm / radius**2) * math.sqrt(kappa_r))
            b1 = 1 - 2 * beta_s**2
            gp = cmath.sqrt(self.eta**2 - beta_s**2)
            gs = cmath.sqrt(1 - beta_s**2)

            numerator = 4 * beta_s * (1 - beta_s**2) * gp
            denominator = b1**2 + 4 * beta_s**2 * gp * gs

            k2 = (numerator / denominator).real / math.pi

        return k2

    def __h2(self, x_norm: float, z_norm: float, time_parameter: float, radius: float, k2_0: float) -> float:
        """
        Imaginary factor for the second integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - time_parameter (float): dimensionless integration time parameter [-]
            - radius (float): dimensionless radius [-]
            - k2_0 (float): initial imaginary factor for the second integral [-]
        """

        pr = time_parameter**2 - radius**2
        return (self.__k2(x_norm, z_norm, time_parameter) - k2_0) / math.sqrt(pr)

    def __f2(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Calculates the second integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: second integral solution [-]
        """

        if abs(x_norm) < self.epsilon:
            x_norm = self.epsilon if x_norm >= 0 else -self.epsilon

        radius = math.sqrt(x_norm**2 + z_norm**2)

        kappa_r = kappa / radius

        time_parameter = radius + self.epsilon

        k2 = self.__k2(x_norm, z_norm, time_parameter)

        if kappa**2 <= time_parameter**2:
            second_integral = 0
        else:
            second_integral = k2 * math.log(kappa_r + math.sqrt(kappa_r**2 - 1))
            h2_1 = 0

            # solve second integral with simpson's rule
            while time_parameter**2 < kappa**2:

                time_parameter += self.integral_stepsize
                h2_2 = self.__h2(x_norm, z_norm, time_parameter, radius, k2)

                time_parameter += self.integral_stepsize
                h2_3 = self.__h2(x_norm, z_norm, time_parameter, radius, k2)

                second_integral += (h2_1 + 4 * h2_2 + h2_3) * self.integral_stepsize / 3
                h2_1 = h2_3

        return second_integral

    def __k3(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Factor for the third integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: imaginary factor for the third integral [-]

        """

        radius = math.sqrt(x_norm**2 + z_norm**2)

        kappa_q = self.eta * abs(x_norm) + z_norm * math.sqrt(1 - self.eta**2)

        if (kappa**2 >= radius**2) or (kappa <= kappa_q) or (x_norm**2 <= self.eta**2 * radius**2):
            k3 = 0
        else:
            beta_q = x_norm * kappa / radius**2 - (z_norm / radius**2) * math.sqrt(radius**2 - kappa**2)

            b2 = (1 - 2 * beta_q**2)**2
            numerator = 4 * beta_q * (1 - beta_q**2) * b2 * math.sqrt(beta_q**2 - self.eta**2)

            denominator = b2**2 + 16 * beta_q**4 * (1 - beta_q**2) * (beta_q**2 - self.eta**2)

            k3 = -numerator / (math.pi * denominator * math.sqrt(radius**2 - kappa**2))

        return k3

    def __f3(self, x_norm: float, z_norm: float, kappa: float) -> float:
        """
        Calculates the third integral

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: third integral solution [-]

        """

        radius = math.sqrt(x_norm**2 + z_norm**2)
        kappa_s = radius
        kappa_q = self.eta * abs(x_norm) + z_norm * math.sqrt(1 - self.eta**2)

        if (kappa <= kappa_q) or (x_norm**2 <= self.eta**2 * radius**2):
            third_integral = 0
        else:

            # kappa_s is the minimum of kappa_s and kappa
            kappa_s = min(kappa, kappa_s)

            # determine integral stepsize as a function of kappa_s and kappa_q
            n_steps = 1000
            stepsize = (kappa_s - kappa_q) / n_steps

            third_integral = 0
            time_parameter = kappa_q + stepsize / 2

            # solve third integral
            while time_parameter < kappa_s:
                third_integral += self.__k3(x_norm, z_norm, time_parameter) * stepsize
                time_parameter += stepsize

        return third_integral

    def calculate_normalised_vertical_stress(self, x_norm: float, z_norm: float, kappa: float):
        """
        Calculate the normalised vertical stress

        Args:
            - x_norm (float): dimensionless x coordinate normalised with the line load length [-]
            - z_norm (float): dimensionless z coordinate normalised with the line load length [-]
            - kappa (float): dimensionless time parameter, normalised with line load length and shear wave velocity [-]

        Returns:
            - float: normalised vertical stress [-]

        """
        normalised_vertical_stress = self.__f1(x_norm + 1, z_norm, kappa) - self.__f1(x_norm - 1, z_norm, kappa) + \
            self.__f2(x_norm + 1, z_norm, kappa) - self.__f2(x_norm - 1, z_norm, kappa) + \
            self.__f3(x_norm + 1, z_norm, kappa) - self.__f3(x_norm - 1, z_norm, kappa)
        return normalised_vertical_stress

    def calculate_vertical_stress(self, x: float, z: float, t: float, load_length: float, load_value: float) -> float:
        """
        Calculate the vertical stress

        Args:
            - x (float): x coordinate [m]
            - z (float): depth from load [m]
            - t (float): time [s]
            - load_length (float): length of the line load [m]
            - load_value (float): value of the line load [N/m]

        Returns:
            - float: vertical stress [Pa]

        """
        x_norm = x / load_length
        z_norm = z / load_length
        kappa = self.cs * t / load_length

        normalised_vertical_stress = self.calculate_normalised_vertical_stress(x_norm, z_norm, kappa)
        vertical_stress = normalised_vertical_stress * load_value * load_length
        return vertical_stress


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    porosity = 0.0
    density_solid = 1020

    line_load_length = 1
    load_value = -1000

    strip_load = StripLoad(2.55e3 * 36, 0.25, (1 - porosity) * density_solid, load_value)

    cs = strip_load.cs

    # dimensionless time parameter

    # end_time = 10 * line_load_length / cs
    end_time = 1

    # dimensionless x coordinate
    start_x = 0
    end_x = 10

    n_steps = 100

    ts = [end_time * i / n_steps for i in range(n_steps)]
    kappas = [cs * t / line_load_length for t in ts]
    kappa = cs * end_time / line_load_length

    xs = [start_x + (end_x - start_x) * i / n_steps for i in range(n_steps)]

    # z coordinate
    x = 5 * line_load_length
    z = 1 * line_load_length

    all_sigma_zz = []

    for x in xs:

        xi = x / line_load_length
        zeta = z / line_load_length

        sigma_zz_normalised = strip_load.calculate_normalised_vertical_stress(xi, zeta, kappa)
        sigma_zz = sigma_zz_normalised * abs(strip_load.load)

        all_sigma_zz.append(sigma_zz)

    plt.plot(xs, all_sigma_zz)

    plt.show()
