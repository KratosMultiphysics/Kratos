import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.structural_component import StructuralComponent
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.panel_buckling import PanelUniaxialBuckling, PanelBiaxialBuckling
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

class Panel(StructuralComponent):

    def __init__(self, 
                 sub_model_part, 
                 boundary_conditions: list[float]):
        super().__init__(sub_model_part, boundary_conditions)
        self.geometry = None
        self.material = None
        self.response = None
        self.load_state = None
        self.analysis_methods = self._CreateAnalysisMethods()
        self.analysis_results = []
    
    @classmethod
    def FromKratosParametersObject(cls, sub_model_part, data):
        """This methods instantiates a panel object from given kratos parameters.

        Args:
            sub_model_part (_type_): SubModelPart defined in the *.mdpa file
            data (_type_): Data from the configuration file (*.json)
        """
        boundary_conditions = [data["boundary_conditions"][i].GetDouble() for i in range(data["boundary_conditions"].size())]

        return cls(sub_model_part, 
            boundary_conditions)
    
    def Initialize(self):
        self.ExtractGeometry()
        self.ExtractMaterial()
    
    def ExtractMaterial(self) -> None:
        if self.sub_model_part.NumberOfElements() == 0:
            raise RuntimeError("Panel submodelpart contains no elements")
        
        element = next(self.sub_model_part.Elements.__iter__())
        properties = element.Properties

        E = properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
        nu = properties.GetValue(KratosMultiphysics.POISSON_RATIO)

        self.material = PanelMaterial(E, nu)

        KratosMultiphysics.Logger.PrintInfo(
        "Panel",
        f"Material extracted: E={E:.6e}, nu={nu:.6f}")   
    
    def ExtractGeometry(self) -> None:

        points = self._CreateNodeCloud()
        centroid, ez, fallback_x = self._FitAveragePlane(points)
        center_element = self._GetCentralElement(centroid)
        ex, ey = self._GetPanelAxes(center_element, ez, fallback_x)
        projected_points, _ = self._ProjectPointsToLocalPlane(points, centroid, ex, ey, ez)

        x_coords = projected_points[:, 0]
        y_coords = projected_points[:, 1]

        a = x_coords.max() - x_coords.min()
        b = y_coords.max() - y_coords.min()

        KratosMultiphysics.Logger.PrintInfo(
        "Panel Dimensions",
        f"a {a}, b {b}")
        aspect_ratio = a / b
        thickness = self._ComputeAverageThickness()
        self.geometry = PanelGeometry(centroid, ex, ey, ez, a, b, aspect_ratio, thickness)

    def ExtractResponse(self) -> None:
        self._RequireGeometry()

        total_volume = 0.0
        sigma_xx_sum = 0.0
        sigma_yy_sum = 0.0
        tau_xy_sum = 0.0

        R = self.coordinate_system

        for element in self.sub_model_part.Elements:
            area = element.GetGeometry().Area()
            thickness = element.Properties.GetValue(KratosMultiphysics.THICKNESS)
            volume = area * thickness

            stress_global = element.CalculateOnIntegrationPoints(
                SMA.SHELL_STRESS_MIDDLE_SURFACE_GLOBAL,
                self.sub_model_part.ProcessInfo
            )[0]

            sigma_global = np.array(stress_global)
            sigma_panel = R @ sigma_global @ R.T

            sigma_xx_sum += sigma_panel[0, 0] * volume
            sigma_yy_sum += sigma_panel[1, 1] * volume
            tau_xy_sum += sigma_panel[0, 1] * volume
            total_volume += volume

        if total_volume <= 0.0:
            raise RuntimeError(
                f"Panel '{self.sub_model_part.Name}' has zero total volume."
            )

        sigma_xx_panel = sigma_xx_sum / total_volume
        sigma_yy_panel = sigma_yy_sum / total_volume
        tau_xy_panel = tau_xy_sum / total_volume

        panel_stress_tensor = np.array([
            [sigma_xx_panel, tau_xy_panel, 0.0],
            [tau_xy_panel, sigma_yy_panel, 0.0],
            [0.0, 0.0, 0.0]
        ])

        self.response = PanelResponse(sigma_xx_panel, sigma_yy_panel, tau_xy_panel, panel_stress_tensor)

        KratosMultiphysics.Logger.PrintInfo(
            "Panel",
            f"Response extracted for '{self.sub_model_part.Name}': "
            f"sigma_xx={sigma_xx_panel:.6e}, "
            f"sigma_yy={sigma_yy_panel:.6e}, "
            f"tau_xy={tau_xy_panel:.6e}"
        )


    def PrepareAnalysis(self) -> None:
        self.ExtractMaterial()
        self.ExtractResponse()
        self.ClassifyLoadState()

    def RunAnalysis(self) -> None:
        self.RunMethods()

    def RunMethods(self) -> None:
        self._RequireLoadState()

        for method in self.analysis_methods:
            if method.IsApplicable(self):
                result = method.Evaluate(self)
                self.analysis_results.append(result)

                KratosMultiphysics.Logger.PrintInfo(
                "Panel",
                f"{result.method_name}: RF={result.value:.6e}"
            )

    def ClassifyLoadState(self) -> None:
        self._RequireResponse()

        sigma_xx = self.response.sigma_xx
        sigma_yy = self.response.sigma_yy
        tau_xy = self.response.tau_xy

        tolerance = 1e-12 * max(abs(sigma_xx), abs(sigma_yy), abs(tau_xy), 1.0)

        has_x_compression = sigma_xx < -tolerance
        has_y_compression = sigma_yy < -tolerance
        has_shear = abs(tau_xy) > tolerance

        is_biaxial_compression = has_x_compression and has_y_compression
        is_uniaxial_compression = (has_x_compression != has_y_compression)

        # This is just for testing purposes and needs to be changed to a "more correct" logic for panels with shear
        is_shear_dominant = has_shear and not is_biaxial_compression and not is_uniaxial_compression

        self.load_state = PanelLoadState(
            has_x_compression,
            has_y_compression,
            has_shear,
            is_uniaxial_compression,
            is_biaxial_compression,
            is_shear_dominant
        )


    def _CreateNodeCloud(self) -> np.ndarray:
        points = np.empty((self.sub_model_part.NumberOfNodes(), 3))
        for i, node in enumerate(self.sub_model_part.Nodes):
            points[i, 0] = node.X
            points[i, 1] = node.Y
            points[i, 2] = node.Z

        KratosMultiphysics.Logger.PrintInfo(
        "Panel",
        f"Node cloud created with shape {points.shape}")
        return points

    def _FitAveragePlane(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        centroid = np.mean(points, axis=0)
        centered_points = points - centroid

        _, _, vh = np.linalg.svd(centered_points, full_matrices=False)
        fallback_x = vh[0]
        ez = vh[-1]/np.linalg.norm(vh[-1])
        KratosMultiphysics.Logger.PrintInfo(
        "Panel",
        f"Average plane fitted: centroid={centroid}, normal={ez}"
    )
        return centroid, ez, fallback_x

    def _GetCentralElement(self, panel_centroid: np.ndarray):
        min_distance = float('inf')
        center_element = None

        for element in self.sub_model_part.Elements:
            coords = [np.array([node.X, node.Y, node.Z]) for node in element.GetGeometry()]
            element_centroid = np.mean(coords, axis=0)
            distance = np.linalg.norm(element_centroid - panel_centroid)

            if distance < min_distance:
                min_distance = distance
                center_element = element

        return center_element

    def _GetPanelAxes(self, center_element, ez: np.ndarray, fallback_x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        local_axis_1 = np.array(center_element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, self.sub_model_part.ProcessInfo))[0]
        ex = local_axis_1 - np.dot(local_axis_1, ez) * ez
        ex_norm = np.linalg.norm(ex)

        if ex_norm < 1e-12:
            ex = fallback_x - np.dot(fallback_x, ez) * ez
            ex_norm = np.linalg.norm(ex)

        ex /= ex_norm
        ey = np.cross(ez, ex)
        ey /= np.linalg.norm(ey)

        # Re-orthonormalize ex to avoid drift
        ex = np.cross(ey, ez)
        ex /= np.linalg.norm(ex)

        return ex, ey

    def _ProjectPointsToLocalPlane(self, points: np.ndarray, centroid: np.ndarray, ex: np.ndarray, ey: np.ndarray, ez: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        centered_points = points - centroid

        x_local = centered_points @ ex
        y_local = centered_points @ ey
        z_local = centered_points @ ez

        projected_points = np.column_stack((x_local, y_local))
        return projected_points, z_local

    def _ComputeAverageThickness(self) -> float:
        weighted_thickness_sum = 0.0
        total_area = 0.0

        for element in self.sub_model_part.Elements:
            area = element.GetGeometry().Area()
            thickness = element.Properties.GetValue(KratosMultiphysics.THICKNESS)

            weighted_thickness_sum += area * thickness
            total_area += area

        if total_area <= 0.0:
            raise RuntimeError("Panel submodelpart has zero total area")

        return weighted_thickness_sum / total_area

    @property
    def a(self):
        self._RequireGeometry()
        return self.geometry.a
    
    @property
    def b(self):
        self._RequireGeometry()
        return self.geometry.b
    
    @property
    def E(self):
        self._RequireMaterial()
        return self.material.young_modulus
    
    @property
    def nu(self):
        self._RequireMaterial()
        return self.material.poisson_ratio
    
    @property
    def t(self):
        self._RequireGeometry()
        return self.geometry.thickness
    
    @property
    def coordinate_system(self):
        self._RequireGeometry()
        return np.vstack((
            self.geometry.x_axis,
            self.geometry.y_axis,
            self.geometry.z_axis
        ))

    def _RequireGeometry(self):
        if self.geometry is None:
            raise RuntimeError(f"Panel '{self.sub_model_part.Name}' has no geometry. Call ExtractGeometry() first.")
        
    def _RequireMaterial(self):
        if self.material is None:
            raise RuntimeError(f"Panel '{self.sub_model_part.Name}' has no material. Call ExtractMaterial() first.")
        
    def _RequireResponse(self):
        if self.response is None:
            raise RuntimeError(f"Panel '{self.sub_model_part.Name}' has no response. Call ExtractResponse() first.")
        
    def _RequireLoadState(self):
        if self.load_state is None:
            raise RuntimeError(f"Panel '{self.sub_model_part.Name}' has no load state. Call ClassifyLoadState() first.")

    def _CreateAnalysisMethods(self):
        return [PanelUniaxialBuckling(), 
                PanelBiaxialBuckling()]
