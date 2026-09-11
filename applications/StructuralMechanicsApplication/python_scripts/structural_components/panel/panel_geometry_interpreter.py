import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelGeometry
import numpy as np

class PanelGeometryInterpreter:

    def Interpret(self, sub_model_part) -> PanelGeometry:
        self._ValidateSubModelPart(sub_model_part)
        points = self._GetPoints(sub_model_part)
        self.centered_points = self._GetCenteredPoints(points)
        panel_base_vectors = self._CalculatePanelCoordinateSystem(self.centered_points)
        self._CheckCurvature(sub_model_part, panel_base_vectors[-1])
        length, width = self._GetPanelDimensions(self.centered_points, panel_base_vectors)

        KratosMultiphysics.Logger.PrintInfo(
        "Panel Dimensions",
        f"a {length:.2f}, b {width:.2f}")
        aspect_ratio = length / width
        thickness = self._ComputeAverageThickness(sub_model_part)
        return PanelGeometry( panel_base_vectors[0], panel_base_vectors[1], panel_base_vectors[2], length, width, aspect_ratio, thickness)

    def _GetPoints(self, sub_model_part) -> np.ndarray:
        points = np.array([[node.X, node.Y, node.Z] for node in sub_model_part.Nodes], dtype=float)
        return points

    def _GetCenteredPoints(self, points):
        centroid = centroid = np.mean(points, axis=0)
        centered_points = points-centroid
        return centered_points

    def _CalculatePanelCoordinateSystem(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculates local panel base vectors using singular value decomposition.

        Args:
            points (np.ndarray): Contains node coordinates.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: Contains base vectors (ex, ey, ez)
        """
        centered_points = self._GetCenteredPoints(points)
        _, _, vh = np.linalg.svd(centered_points, full_matrices=False)
        ez = vh[-1]
        ex = vh[0]
        ey = np.cross(ez, ex)

        return (ex, ey, ez)

    def _GetPanelDimensions(self, centered_points, base_vectors):
        """Calculates the panel dimensions using the PCA-based panel base vectors,

        Args:
            centered_points (np.ndarray): points (node coordinates) - centroid
            base_vectors (tuple[np.ndarray, np.ndarray, np.ndarray]): (ex, ey, ez)

        Returns:
            length, width: Panel dimensions
        """
        ex = base_vectors[0]
        ey = base_vectors[1]
        length_coordinates = centered_points @ ex
        width_coordinates = centered_points @ ey

        length = np.ptp(length_coordinates)
        width = np.ptp(width_coordinates)

        return length, width

    def _ComputeAverageThickness(self, sub_model_part) -> float:
        """Loops through the elements that are part of the panel submodelpart and calculates the average thickness.

        Args:
            sub_model_part: Panel submodelpart containing the shell elements.

        Returns:
            float: Average thickness of the panel
        """
        total_area = sum(element.GetGeometry().Area() for element in sub_model_part.Elements)

        weighted_thickness = sum(
            element.GetGeometry().Area() * element.Properties.GetValue(KratosMultiphysics.THICKNESS)
            for element in sub_model_part.Elements
        )

        return weighted_thickness / total_area

    def _GetElementNormal(self, element) -> np.ndarray:
        geometry = element.GetGeometry()
        number_of_integration_points = geometry.IntegrationPointsNumber()

        normals = np.array([geometry.UnitNormal(i) for i in range(number_of_integration_points)])

        normal = np.mean(normals, axis=0)
        return normal / np.linalg.norm(normal)

    def _CheckCurvature(self, sub_model_part, ez: np.ndarray) -> None:
        """Warn if element normals deviate from the average panel normal.

        The check compares each element normal with the PCA-based panel normal direction. A large maximum angle indivates that the panel is curved or otherwise not well represented by a single flat local coordinate system.

        Args:
            sub_model_part (_type_): Panel submodelpart containing the shell elements.
            ez (np.ndarray): PCA-based average panel normal direction.
        """
        #TODO: Just gives a warning right now. Curved panels should eventually be handled by a dedicated geometry interpretation algorithm.
        angle_tolerance_degrees = 1.0

        element_normals = np.array([self._GetElementNormal(element) for element in sub_model_part.Elements])

        cos_angles = np.abs(element_normals @ ez)
        max_angle_degrees = np.degrees(np.arccos(np.clip(cos_angles, -1.0, 1.0)).max())

        if max_angle_degrees > angle_tolerance_degrees:
            KratosMultiphysics.Logger.PrintWarning(
            "Panel",
            f"Panel submodelpart '{sub_model_part.Name}' appears to be curved. "
            f"Maximum element-normal deviation from the average panel normal is "
            f"{max_angle_degrees:.2f} degrees."
        )
    
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