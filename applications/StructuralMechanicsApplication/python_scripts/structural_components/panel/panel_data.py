import numpy as np
from dataclasses import dataclass

@dataclass
class PanelGeometry:
    centroid: np.ndarray
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