import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelGeometry
import numpy as np

class PanelGeometryInterpreter:

    def Interpret(self, sub_model_part) -> PanelGeometry:
        self._ValidateSubModelPart(sub_model_part)
        points = self._CreateNodeCloud(sub_model_part)
        centroid, ez, fallback_x = self._FitAveragePlane(points)
        center_element = self._GetCentralElement(centroid, sub_model_part)
        ex, ey = self._GetPanelAxes(center_element, ez, fallback_x, sub_model_part)
        projected_points, _ = self._ProjectPointsToLocalPlane(points, centroid, ex, ey, ez)

        x_coords = projected_points[:, 0]
        y_coords = projected_points[:, 1]

        a = x_coords.max() - x_coords.min()
        b = y_coords.max() - y_coords.min()

        if a <= 0.0 or b <= 0.0:
            raise RuntimeError(
                f"Panel submodelpart '{sub_model_part.Name}' has degenerate dimensions: "
                f"a={a}, b={b}."
            )

        KratosMultiphysics.Logger.PrintInfo(
        "Panel Dimensions",
        f"a {a}, b {b}")
        aspect_ratio = a / b
        thickness = self._ComputeAverageThickness(sub_model_part)
        return PanelGeometry(centroid, ex, ey, ez, a, b, aspect_ratio, thickness)

    def _CreateNodeCloud(self, sub_model_part) -> np.ndarray:
        points = np.empty((sub_model_part.NumberOfNodes(), 3))
        for i, node in enumerate(sub_model_part.Nodes):
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

    def _GetCentralElement(self, panel_centroid: np.ndarray, sub_model_part):
        min_distance = float('inf')
        center_element = None

        for element in sub_model_part.Elements:
            coords = [np.array([node.X, node.Y, node.Z]) for node in element.GetGeometry()]
            element_centroid = np.mean(coords, axis=0)
            distance = np.linalg.norm(element_centroid - panel_centroid)

            if distance < min_distance:
                min_distance = distance
                center_element = element

        if center_element is None:
            raise RuntimeError(
                f"Could not find a central element for panel submodelpart "
                f"'{sub_model_part.Name}'."
            )

        return center_element

    def _GetPanelAxes(self, center_element, ez: np.ndarray, fallback_x: np.ndarray, sub_model_part) -> tuple[np.ndarray, np.ndarray]:
        local_axis_1 = np.array(center_element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, sub_model_part.ProcessInfo))[0]
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

    def _ComputeAverageThickness(self, sub_model_part) -> float:
        weighted_thickness_sum = 0.0
        total_area = 0.0

        for element in sub_model_part.Elements:
            area = element.GetGeometry().Area()
            thickness = element.Properties.GetValue(KratosMultiphysics.THICKNESS)

            weighted_thickness_sum += area * thickness
            total_area += area

        if total_area <= 0.0:
            raise RuntimeError("Panel submodelpart has zero total area")

        return weighted_thickness_sum / total_area
    
    def _ValidateSubModelPart(self, sub_model_part):
        if sub_model_part.NumberOfNodes() == 0:
            raise RuntimeError(
                f"Panel submodelpart '{sub_model_part.Name}' contains no nodes."
            )

        if sub_model_part.NumberOfElements() == 0:
            raise RuntimeError(
                f"Panel submodelpart '{sub_model_part.Name}' contains no elements."
            )

        if sub_model_part.NumberOfNodes() < 3:
            raise RuntimeError(
                f"Panel submodelpart '{sub_model_part.Name}' needs at least 3 nodes "
                f"to define a panel plane, but has {sub_model_part.NumberOfNodes()}."
            )