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
