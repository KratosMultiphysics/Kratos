import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods import StrengthMethods, StabilityMethods
from KratosMultiphysics.StructuralMechanicsApplication.structural_elements.structural_element import StructuralElement
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

class Panel(StructuralElement):

    def __init__(self, 
                 sub_model_part, 
                 boundary_conditions: list[float],
                  analysis_methods: list[str]):
        super().__init__(sub_model_part, boundary_conditions, analysis_methods)
    
    @classmethod
    def FromKratosParametersObject(cls, sub_model_part, data):
        """This methods instantiates a panel object from given kratos parameters.

        Args:
            sub_model_part (_type_): SubModelPart defined in the *.mdpa file
            data (_type_): Data from the configuration file (*.json)
        """
        method_list = [data["analysis_methods"][i].GetString() for i in range(data["analysis_methods"].size())]
        boundary_conditions = [data["boundary_conditions"][i].GetDouble() for i in range(data["boundary_conditions"].size())]

        return cls(sub_model_part, 
            boundary_conditions,
            method_list)
    
    def ExtractMaterial(self) -> None:
        if self.sub_model_part.NumberOfElements() == 0:
            raise RuntimeError("Panel submodelpart contains no elements")
        
        element = next(self.sub_model_part.Elements.__iter__())
        properties = element.Properties

        self.E = properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
        self.nu = properties.GetValue(KratosMultiphysics.POISSON_RATIO)
        element_thickness = []
        for elem in self.sub_model_part.Elements:
            element_thickness.append(elem.Properties.GetValue(KratosMultiphysics.THICKNESS))

        self.thickness = sum(element_thickness)/len(element_thickness)

        KratosMultiphysics.Logger.PrintInfo(
        "Panel",
        f"Material extracted: E={self.E:.6e}, nu={self.nu:.6f}")   
    
    def ExtractGeometry(self) -> None:

        points = self._CreateNodeCloud()
        centroid, ez, fallback_x = self._FitAveragePlane(points)
        center_element = self._GetCentralElement(centroid)
        ex, ey = self._GetPanelAxes(center_element, ez, fallback_x)
        projected_points, local_z = self._ProjectPointsToLocalPlane(points, centroid, ex, ey, ez)

        x_coords = projected_points[:, 0]
        y_coords = projected_points[:, 1]

        a = x_coords.max() - x_coords.min()
        b = y_coords.max() - y_coords.min()

        self.panel_centroid = centroid
        self.center_element = center_element
        self.x_axis_base_vector = ex
        self.y_axis_base_vector = ey
        self.z_axis_base_vector = ez
        self.a = a
        self.b = b
        KratosMultiphysics.Logger.PrintInfo(
        "Panel Dimensions",
        f"a {a}, b {b}")
        self.aspect_ratio = a / b
        self.thickness = self._ComputeAverageThickness()
        self.flatness_max_distance = np.max(np.abs(local_z))
        self.flatness_rms_distance = np.sqrt(np.mean(local_z**2))
        self.projected_points_2d = projected_points

    def ExtractResponse(self) -> None:
        if not hasattr(self, "x_axis_base_vector") or not hasattr(self, "y_axis_base_vector") or not hasattr(self, "z_axis_base_vector"):
            raise RuntimeError(
                f"Panel '{self.sub_model_part.Name}' has no local coordinate system. "
                "Call ExtractGeometry() before ExtractResponse()."
            )

        total_volume = 0.0
        sigma_xx_sum = 0.0
        sigma_yy_sum = 0.0
        tau_xy_sum = 0.0

        R = np.vstack((
            self.x_axis_base_vector,
            self.y_axis_base_vector,
            self.z_axis_base_vector
        ))

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

        self.sigma_xx_panel = sigma_xx_sum / total_volume
        self.sigma_yy_panel = sigma_yy_sum / total_volume
        self.tau_xy_panel = tau_xy_sum / total_volume

        self.panel_stress_tensor = np.array([
            [self.sigma_xx_panel, self.tau_xy_panel, 0.0],
            [self.tau_xy_panel, self.sigma_yy_panel, 0.0],
            [0.0, 0.0, 0.0]
        ])

        KratosMultiphysics.Logger.PrintInfo(
            "Panel",
            f"Response extracted for '{self.sub_model_part.Name}': "
            f"sigma_xx={self.sigma_xx_panel:.6e}, "
            f"sigma_yy={self.sigma_yy_panel:.6e}, "
            f"tau_xy={self.tau_xy_panel:.6e}"
        )


    def PrepareAnalysis(self) -> None:
        self.ExtractGeometry()
        self.ExtractMaterial()
        self.ExtractResponse()

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
