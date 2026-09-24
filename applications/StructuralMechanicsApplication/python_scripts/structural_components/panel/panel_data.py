import numpy as np
from dataclasses import dataclass

@dataclass
class PanelGeometry:
    x_axis: np.ndarray
    y_axis: np.ndarray
    z_axis: np.ndarray
    a: float # Dimension along local panel x coordinate
    b: float # Dimension along local panel y coordinate
    aspect_ratio: float
    thickness: float
    curvature_type: str = "flat"

@dataclass
class PanelMaterial:
    young_modulus: float
    poisson_ratio: float

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
class PanelResponse:
    sigma_xx: float # stress along local panel x coordinate
    sigma_yy: float # stress along local panel y coordinate
    tau_xy: float
    stress_tensor: np.ndarray

@dataclass
class PanelLoadState:
    has_x_compression: bool
    has_y_compression: bool
    has_shear: bool
    is_uniaxial_compression: bool
    is_biaxial_compression: bool
    is_shear_dominant: bool

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
    