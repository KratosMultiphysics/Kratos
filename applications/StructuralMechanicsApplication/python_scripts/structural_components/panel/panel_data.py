import numpy as np
from dataclasses import dataclass

@dataclass
class PanelGeometry:
    x_axis: np.ndarray
    y_axis: np.ndarray
    z_axis: np.ndarray
    a: float # Dimension along local panel x coordinate
    b: float # Dimension along local panel y coordinate
    material_axis_offset: float

    @property
    def rotation_matrix(self) -> np.ndarray:
        """Rows are the panel axes: v_panel = R @ v_global, T_panel = R @ T_global @ R.T"""
        return np.vstack((self.x_axis, self.y_axis, self.z_axis))

@dataclass
class PanelCompositeMaterial:
    #TODO: Add material properties from composite material
    element_id      : np.ndarray
    ply_id          : np.ndarray
    ply_side        : np.ndarray
    ply_thickness   : np.ndarray
    ply_angle       : np.ndarray

    strength_R_pa_t : np.ndarray
    strength_R_pa_c : np.ndarray
    strength_R_tr_t : np.ndarray
    strength_R_tr_c : np.ndarray
    strength_R_trpa : np.ndarray

    inclination_p_trtr_c : np.ndarray
    inclination_p_trpa_t : np.ndarray
    inclination_p_trpa_c : np.ndarray

    youngs_modul_E_pa    : np.ndarray
    youngs_modul_E_tr    : np.ndarray
    poissons_ratio_nu_12 : np.ndarray
    shear_modulus_G_trpa : np.ndarray

    degradationfactor_E_tr_A  : float
    degradationfactor_G_patr_A: float
    degradationfactor_E_tr_B  : float
    degradationfactor_G_patr_B: float

@dataclass
class PanelCompositeResponse:
    #TODO: Add responses that are necessary for Puck analysis
    sigma_1 : np.ndarray
    sigma_2 : np.ndarray
    tau_21  : np.ndarray

    force_x   : np.ndarray
    force_y   : np.ndarray
    force_xy  : np.ndarray

    moment_x   : np.ndarray
    moment_y   : np.ndarray
    moment_xy  : np.ndarray

    plies_element : int 
    num_elements  : int

@dataclass
class PanelLoads:
    """Running loads per element in the panel frame. Kratos sign convention: tension positive."""
    element_ids: np.ndarray # (n,)
    element_areas: np.ndarray # (n,)
    centroids: np.ndarray # (n, 2): xi = x/a, eta = y/b in [0, 1]
    N: np.ndarray # (n, 3): N_xx, N_yy, N_xy [force/length]
    M: np.ndarray # (n, 3): M_xx, M_yy, M_xy [moment/length]

    @property
    def N_average(self) -> np.ndarray:
        """Area-weighted panel average (N_xx, N_yy, N_xy), input for buckling."""
        return self.element_areas @ self.N / self.element_areas.sum()

    @property
    def N_linear(self) -> np.ndarray:
        """Area-weighted least-squares fit N(xi, eta) = c0 + c1 (xi - 1/2) + c2 (eta - 1/2).

        Returns (3, 3): rows c0, c1, c2; columns N_xx, N_yy, N_xy.
        """
        basis = np.column_stack((np.ones(len(self.N)), self.centroids - 0.5))
        sqrt_w = np.sqrt(self.element_areas)[:, None]
        coefficients, *_ = np.linalg.lstsq(sqrt_w * basis, sqrt_w * self.N, rcond=None)
        return coefficients
