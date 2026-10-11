import numpy as np
from dataclasses import dataclass

@dataclass(frozen=True)
class PlyStrength:
    R_pa_t: float       # X_t: tension parallel to fibre
    R_pa_c: float       # X_c: compression parallel to fibre
    R_tr_t: float       # Y_t: tension transverse to fibre
    R_tr_c: float       # Y_c: compression transverse to fibre
    R_trpa: float       # S_12: in-plane shear

@dataclass(frozen=True)
class Ply:
    thickness: float
    angle: float                        # [deg] from the panel x_axis about the panel z_axis
    E1: float
    E2: float
    nu12: float
    G12: float
    strength: PlyStrength | None = None  # None for isotropic material

    @property
    def Q(self) -> np.ndarray:
        """Reduced stiffness in ply coordinates (1, 2)."""
        nu21 = self.nu12 * self.E2 / self.E1
        d = 1.0 - self.nu12 * nu21
        return np.array([[self.E1 / d,             self.nu12 * self.E2 / d, 0.0],
                         [self.nu12 * self.E2 / d, self.E2 / d,             0.0],
                         [0.0,                     0.0,                     self.G12]])

    @property
    def Q_bar(self) -> np.ndarray:
        """Reduced stiffness in panel coordinates (x, y)."""
        Q11, Q12, Q22, Q66 = self.Q[0, 0], self.Q[0, 1], self.Q[1, 1], self.Q[2, 2]
        c, s = np.cos(np.radians(self.angle)), np.sin(np.radians(self.angle))
        Qb11 = Q11*c**4 + 2*(Q12 + 2*Q66)*s**2*c**2 + Q22*s**4
        Qb22 = Q11*s**4 + 2*(Q12 + 2*Q66)*s**2*c**2 + Q22*c**4
        Qb12 = (Q11 + Q22 - 4*Q66)*s**2*c**2 + Q12*(c**4 + s**4)
        Qb66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*s**2*c**2 + Q66*(s**4 + c**4)
        Qb16 = (Q11 - Q12 - 2*Q66)*c**3*s - (Q22 - Q12 - 2*Q66)*s**3*c
        Qb26 = (Q11 - Q12 - 2*Q66)*c*s**3 - (Q22 - Q12 - 2*Q66)*s*c**3
        return np.array([[Qb11, Qb12, Qb16],
                         [Qb12, Qb22, Qb26],
                         [Qb16, Qb26, Qb66]])

    @property
    def T(self) -> np.ndarray:
        """Stress transformation from panel (x, y) to ply (1, 2) coordinates: sigma_12 = T @ sigma_xy."""
        c, s = np.cos(np.radians(self.angle)), np.sin(np.radians(self.angle))
        return np.array([[c * c,  s * s,  2.0 * c * s],
                            [s * s,  c * c, -2.0 * c * s],
                            [-c * s, c * s,  c * c - s * s]])


class Laminate:
    """Classical laminate theory about the mid-plane.

    Plies are ordered from top (z = +h/2, along the panel z_axis) to bottom,
    as in Kratos SHELL_ORTHOTROPIC_LAYERS.
    """

    def __init__(self, plies: list[Ply]):
        self.plies = plies
        self.thickness = sum(ply.thickness for ply in plies)
        z = self.thickness / 2.0 - np.cumsum([0.0] + [ply.thickness for ply in plies])
        self.z_top, self.z_bottom = z[:-1], z[1:]

        self.A = sum(ply.Q_bar * (zt - zb)            for ply, zt, zb in zip(plies, self.z_top, self.z_bottom))
        self.B = sum(ply.Q_bar * (zt**2 - zb**2) / 2  for ply, zt, zb in zip(plies, self.z_top, self.z_bottom))
        self.D = sum(ply.Q_bar * (zt**3 - zb**3) / 3  for ply, zt, zb in zip(plies, self.z_top, self.z_bottom))

    @property
    def ABD(self) -> np.ndarray:
        return np.block([[self.A, self.B], [self.B, self.D]])

    def IsSymmetric(self, rtol: float = 1e-8) -> bool:
        """B = 0 within tolerance (B scales with A * thickness)."""
        return np.linalg.norm(self.B) <= rtol * np.linalg.norm(self.A) * self.thickness

    @classmethod
    def FromIsotropic(cls, E: float, nu: float, thickness: float) -> "Laminate":
        return cls([Ply(thickness, 0.0, E, E, nu, E / (2.0 * (1.0 + nu)))])

    @classmethod
    def FromOrthotropicLayers(cls, layers: np.ndarray, angle_offset: float) -> "Laminate":
        """layers: rows of SHELL_ORTHOTROPIC_LAYERS; angle_offset [deg]: material axis offset of the panel."""
        return cls([Ply(thickness=row[0], angle=row[1] + angle_offset,
                        E1=row[3], E2=row[4], nu12=row[5], G12=row[6],
                        strength=PlyStrength(*row[9:14]))
                    for row in layers])

    def PlyStresses(self, N: np.ndarray, M: np.ndarray) -> np.ndarray:
        """Ply stresses (sigma_1, sigma_2, tau_12) in ply coordinates at the top and bottom face of each ply.

        N, M: (n, 3) running loads in the panel frame.
        Returns (n, n_plies, 2, 3); third index: 0 = top face, 1 = bottom face of the ply.
        """
        strains = np.linalg.solve(self.ABD, np.hstack((N, M)).T).T     # (n, 6): eps_0, kappa
        stresses = np.empty((N.shape[0], len(self.plies), 2, 3))
        for k, ply in enumerate(self.plies):
            to_ply = ply.T @ ply.Q_bar
            for side, z in enumerate((self.z_top[k], self.z_bottom[k])):
                stresses[:, k, side] = (strains[:, :3] + z * strains[:, 3:]) @ to_ply.T
        return stresses
