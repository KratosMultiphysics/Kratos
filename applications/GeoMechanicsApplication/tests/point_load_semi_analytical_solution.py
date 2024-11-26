import math

class PointLoad:
    """
    Class for solving the point load problem semi-analytically

    A point load is applied to the surface of an elastic half-space dynamically. The solution is based on the following
    publication: An Introduction to Soil Dynamics , Verruijt A., 2009, Delft University of Technology, Chapter 13.2

    Attributes:
        - youngs_modulus (float): Young's modulus [Pa]
        - poissons_ratio (float): Poisson's ratio [-]
        - density (float): density [kg/m^3]
        - load (float): load value [N]
        - shear_modulus (float): shear modulus [Pa]
        - cs (float): shear wave velocity [m/s]

    """

    def __init__(self,
                 youngs_modulus: float,
                 poissons_ratio: float,
                 density: float,
                 load: float,
                 integral_stepsize: float = 0.001):
        """
        Constructor for the point load class

        Args:
            - youngs_modulus (float): Young's modulus [Pa]
            - poissons_ratio (float): Poisson's ratio [-]
            - density (float): density [kg/m^3]
            - load (float): load value [N]
        """

        self.youngs_modulus = youngs_modulus
        self.poissons_ratio = poissons_ratio
        self.density = density
        self.load = load

        # dimensionless time step size for solving integrals
        self.integral_stepsize = integral_stepsize

        # calculate derived elastic properties
        self.shear_modulus = self.youngs_modulus / (2 * (1 + self.poissons_ratio))
        self.cs = math.sqrt(self.shear_modulus / self.density)  # shear wave velocity


    def __wpekeris(self, nu: float, t: float) -> float:
        fac = 1.0 / (2.0 * math.pi * math.pi)
        eps = 1.0e-12
        eps1 = 0.001
        nn = (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu))
        n = math.sqrt(nn)
        e = 1.0e-12
        f = 1.0
        b = (1.0 - nu) / 8.0

        if nu > 0.1:
            while f > e:
                a = b
                b = (1.0 - nu) / (8.0 * (1.0 + a) * (nu + a))
                f = math.fabs(b - a)
        else:
            while f > e:
                a = b
                b = math.sqrt((1.0 - nu) / (8.0 * (1.0 + a) * (1.0 + nu / a)))
                f = math.fabs(b - a)

        tr = math.sqrt(1.0 + b)
        if t <= n:
            g = 0.0
        elif t <= 1.0:
            tt = t * t
            xa = 0.0
            xb = math.pi / 2.0
            k = int(1000.0 * (xb - xa))
            dx = (xb - xa) / k
            fa = 0.0
            g = 0.0
            for j in range(k+1):
                s = math.sin(xa + j * dx)
                ss = s * s
                yy = nn + (tt - nn) * ss
                a = (1.0 - 2.0 * yy) * (1.0 - 2.0 * yy) * (tt - nn) * ss
                b = (1.0 + 8.0 * yy * (-1.0 + yy * (3.0 - 2.0 * nn - 2.0 * yy * (1.0 - nn))))
                f = a / b
                if j > 0: 
                    g -= fac * (f + fa) * dx
                fa = f
        elif t > tr:
            ta = t
            if t < (tr - eps):
                ta = tr + eps1
            tt = ta * ta
            xr = math.asin(math.sqrt((tr * tr - nn) / (tt - nn)))
            xa = 0.0
            xb = xr - eps1
            k = int(1000.0 * (xb - xa))
            dx = (xb - xa) / k
            fa = 0.0
            g = 0.0
            for j in range(k+1):
                s = math.sin(xa + j * dx)
                ss = s * s
                yy = nn + (tt - nn) * ss
                a = (1.0 - 2.0 * yy) * (1.0 - 2.0 * yy) * (tt - nn) * ss
                b = (1.0 + 8.0 * yy * (-1.0 + yy * (3.0 - 2.0 * nn - 2.0 * yy * (1.0 - nn))))
                f = a / b
                if j > 0:
                    g -= fac * (f + fa) * dx
                fa = f

            xa = xr + eps1
            xb = math.pi / 2.0
            k = int(1000.0 * (xb - xa))
            dx = (xb - xa) / k
            for j in range(k+1):
                s = math.sin(xa + j * dx)
                ss = s * s
                yy = nn + (tt - nn) * ss
                a = (1.0 - 2.0 * yy) * (1.0 - 2.0 * yy) * (tt - nn) * ss
                b = (1.0 + 8.0 * yy * (-1.0 + yy * (3.0 - 2.0 * nn - 2.0 * yy * (1.0 - nn))))
                f = a / b
                if j > 0: 
                    g -= fac * (f + fa) * dx
                fa = f

            xr = math.asin(math.sqrt((tr * tr - 1.0) / (tt - 1.0)))
            xa = 0.0
            xb = xr - eps1
            k = int(1000.0 * (xb - xa))
            dx = (xb - xa) / k
            for j in range(k+1):
                s = math.sin(xa + j * dx)
                ss = s * s
                yy = 1.0 + (tt - 1.0) * ss
                a = 4.0 * yy * (tt - 1.0) * (yy - nn) * ss
                b = (1.0 + 8.0 * yy * (-1.0 + yy * (3.0 - 2.0 * nn - 2.0 * yy * (1.0 - nn))))
                f = a / b
                if j > 0:
                    g -= fac * (f + fa) * dx
                fa = f

            xa = xr + eps1
            xb = math.pi / 2.0
            k = int(1000.0 * (xb - xa))
            dx = (xb - xa) / k
            for j in range(k+1):
                s = math.sin(xa + j * dx)
                ss = s * s
                yy = 1.0 + (tt - 1.0) * ss
                a = 4.0 * yy * (tt - 1.0) * (yy - nn) * ss
                b = (1.0 + 8.0 * yy * (-1.0 + yy * (3.0 - 2.0 * nn - 2.0 * yy * (1.0 - nn))))
                f = a / b
                if j > 0:
                    g -= fac * (f + fa) * dx
                fa = f
        else:
            ta = t
            if t > (tr - eps):
                ta = tr - eps1
            tt = ta * ta
            xa = 0.0
            xb = math.pi / 2.0
            k = int(1000.0 * (xb - xa))
            dx = (xb - xa) / k
            g = 0.0
            for j in range(k+1):
                s = math.sin(xa + j * dx)
                ss = s * s
                yy = nn + (tt - nn) * ss
                a = (1.0 - 2.0 * yy) * (1.0 - 2.0 * yy) * (tt - nn) * ss
                b = (1.0 + 8.0 * yy * (-1.0 + yy * (3.0 - 2.0 * nn - 2.0 * yy * (1.0 - nn))))
                f = a / b
                if j > 0:
                    g -= fac * (f + fa) * dx
                fa = f

            for j in range(k+1):
                s = math.sin(xa + j * dx)
                ss = s * s
                yy = 1.0 + (tt - 1.0) * ss
                a = 4.0 * yy * (tt - 1.0) * (yy - nn) * ss
                b = (1.0 + 8.0 * yy * (-1.0 + yy * (3.0 - 2.0 * nn - 2.0 * yy * (1.0 - nn))))
                f = a / b
                if j > 0:
                    g -= fac * (f + fa) * dx
                fa = f

        return g


    def calculate_vertical_displacement(self, radius: float, t: float, load_value: float) -> float:
        """
        Calculate the vertical displacement

        Args:
            - radius (float): horizontal distance to the point load [m]
            - t (float): time [s]
            - load_value (float): value of the point load [N]

        Returns:
            - float: vertical displacement [m]
        """

        return self.__wpekeris(self.poissons_ratio, self.cs * t / radius)
