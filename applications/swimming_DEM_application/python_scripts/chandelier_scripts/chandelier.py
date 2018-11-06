import math
import cmath
import mpmath
import matplotlib.pyplot as plt
import numpy as np
import chandelier_parameters as pp

class AnalyticSimulator:
    def __init__(self, pp):
        self.pp = pp
        self.CalculateNonDimensionalVars()
        self.CalculateABC()
        self.CalculateRoots()
        self.CalculateAs()

    def CalculateNonDimensionalVars(self):
        pp = self.pp
        rho_f = float(pp.rho_f)
        rho_p = float(pp.rho_p)
        a = float(pp.a)
        nu = float(pp.nu)
        R = float(pp.R)
        g = float(pp.g)
        omega = float(pp.omega)
        x0 = float(pp.x0)
        y0 = float(pp.y0)
        z0 = float(pp.z0)
        u0 = float(pp.u0)
        v0 = float(pp.v0)
        volume = 4.0 / 3.0 * math.pi * a ** 3

        gamma = rho_p / rho_f
        S = a ** 2 * omega / (9.0 * nu) # pseudo Stokes' Number
        Fr = (R * omega) ** 2 / (g * R) # pseudo Froude's Number
        Re_T = 4.0 / 9.0 * (gamma - 1.0) * g * a ** 3 / nu ** 2 # terminal velocity Reynold's number
        nDw0 = 2.0 * (1.0 - gamma) * S / Fr # Non-dimensional terminal sedimentation velocity

        Z0 = (x0 + 1j * y0) / R
        U0 = 1j * Z0
        #U0 = Z0(u0 + 1j * v0) / (R * omega)
        NDz0 = z0 / R
        self.basset_dimensional_coeff = rho_f * volume * R * omega ** 2
        self.gamma = gamma
        self.S = S
        self.Fr = Fr
        self.Re_T = Re_T
        self.U0 = U0
        self.Z0 = Z0
        self.NDz0 = NDz0
        self.NDw0 = nDw0

    def CalculateABC(self):
        include_lift = self.pp.include_lift
        include_history_force = self.pp.include_history_force

        Cs = self.pp.Cs
        gamma = self.gamma
        S = self.S

        if include_lift:
            A = 1.0/ S - 1.j * Cs / (math.pi * math.sqrt(2 * S))
            B = 3.0 - Cs / (math.pi * math.sqrt(2.0 * S)) - 1.j / S
        else:
            A = 1.0/ S
            B = 3.0 - 1.j / S

        C = - 3. / math.sqrt(math.pi * S)

        A /= (2.0 * gamma + 1)
        B /= (2.0 * gamma + 1)
        C /= (2.0 * gamma + 1)

        if not include_history_force:
            C = 0.0

        self.A = A
        self.B = B
        self.C = C
        self.C1 = - C * (gamma + 0.5)

    def CalculateRoots(self):
        A = self.A
        B = self.B
        C = self.C

        from numpy.polynomial import Polynomial as P
        x0 = B
        x1 = 1.j * C * math.sqrt(math.pi)
        x2 = A
        x3 = - C * math.sqrt(math.pi)
        x4 = 1.0
        p = P([x0, x1, x2, x3, x4])

        self.roots = p.roots()

    def GetProductOfRootDifferences(self, roots, i):
        product = 1.0

        for j in range(4):

            if not (j == i):
                product *= (roots[i] - roots[j])

        return product

    def CalculateAs(self):
        X = self.roots
        U0 = self.U0
        Z0 = self.Z0
        B = self.B
        C = self.C

        A_i_list = 4 * [None]

        for i in range(4):
            A_i_list[i] = (U0 * (X[i] ** 2 - C * math.sqrt(math.pi) * X[i]) - B * Z0) / self.GetProductOfRootDifferences(X, i)

        self.As = A_i_list

    def CalculatePosition(self, NDcoors, t, NDvel = None):
        A = self.As
        X = self.roots
        NDz0 = self.NDz0
        NDw0 = self.NDw0

        NDxy = 0.
        suma = sum([A[i] * X[i] for i in range(4)])

        for i in range(4):
            NDxy += A[i] / X[i] * cmath.exp(X[i] ** 2 * t) * mpmath.erfc(- X[i] * math.sqrt(t))

        NDcoors[0] = float(NDxy.real)
        NDcoors[1] = float(NDxy.imag)
        NDcoors[2] = float(NDz0 + t * NDw0)

        if NDvel is None:
            NDxy = 0.0
            for i in range(4):
                NDxy += A[i] * X[i] * cmath.exp(X[i] ** 2 * t) * mpmath.erfc(- X[i] * math.sqrt(t))
            NDvel[0] = float(NDxy.real)
            NDvel[1] = float(NDxy.imag)
            NDvel[2] = float(NDw0)

    def CalculateTrajectory(self):
        pp = self.pp
        n = pp.n_t_steps + 1 # number of time instants (including t = final_time)

        self.times = [pp.final_time * i / pp.n_t_steps for i in range(n)]
        self.NDtimes = [t * pp.omega for t in self.times]
        self.NDx = [None] * n
        self.NDy = [None] * n
        self.NDz = [None] * n
        NDcoors  = [None] * 3

        for k in range(n):
            self.CalculatePosition(NDcoors, self.NDtimes[k])
            self.NDx[k] = NDcoors[0]
            self.NDy[k] = NDcoors[1]
            self.NDz[k] = NDcoors[2]

        self.x = [value * pp.R for value in self.NDx]
        self.y = [value * pp.R for value in self.NDy]
        self.z = [value * pp.R for value in self.NDz]

    def CalculateBassetForce(self, FB, t):
        A = self.As
        X = self.roots
        C1 = self.C1
        basset_dimensional_coeff = self.basset_dimensional_coeff
        pp = self.pp
        sqrt_t = math.sqrt(t)
        FhZ = sum([(1j * A[i] / X[i] - A[i] * X[i]) * X[i] * cmath.exp(X[i] ** 2 * t) * mpmath.erfc(- X[i] * sqrt_t) for i in range(len(X))])
        FhZ *= basset_dimensional_coeff * C1 * math.sqrt(math.pi)
        FB[0] = float(FhZ.real)
        FB[1] = float(FhZ.imag)
        FB[2] = 0.0

if __name__ == "__main__":
    sim = AnalyticSimulator(pp)
    sim.CalculateTrajectory()
    #for i in range(sim.n):
        #line = [sim.times[i], sim.x[i], sim.y[i]]
        #print(line)
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sim.x, sim.y, sim.z, label='all forces')

    pp.include_history_force = 0
    sim = AnalyticSimulator(pp)
    sim.CalculateTrajectory()
    ax.plot(sim.x, sim.y, sim.z, label='no history force')
    # Create cubic bounding box to simulate equal aspect ratio

    max_range = np.array([max(sim.x)-min(sim.x), max(sim.y)-min(sim.y), max(sim.z)-min(sim.z)]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(max(sim.x)+min(sim.x))
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(max(sim.y)+min(sim.y))
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(max(sim.z)+min(sim.z))
    ax.set_xlabel('$x$ ($m$)', fontsize=20, rotation = 0)
    ax.set_ylabel('$y$ ($m$)', fontsize=20, rotation = 0)
    ax.set_zlabel('$z$ ($m$)', fontsize=20, rotation = 0)
    plt.legend(loc='upper left')
    fig.tight_layout()
    plt.savefig('spirals.eps', format='eps', dpi=1000)
    plt.show()
