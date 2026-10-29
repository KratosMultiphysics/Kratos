import numpy as np
from collections import Counter
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelGeometry

class PanelGeometryInterpreter:
    """Coordinate system and dimensions of a flat, rectangular panel from its shell elements.

    z_axis: plane normal, oriented like the element normals (ply angles are defined about it)
    x_axis: along the longer edge; for square panels the edge closest to the material axis
    """
    curvature_tolerance_degrees = 1.0       # warning if element normals deviate more from the panel normal
    material_axis_tolerance_degrees = 1.0   # material axes of all elements must agree within this angle
    rectangle_tolerance = 0.01              # relative deviation of the element area from a * b

    def Interpret(self, sub_model_part) -> PanelGeometry:
        centered_points = self._GetCenteredPoints(sub_model_part)
        process_info = sub_model_part.ProcessInfo

        z_axis = self._GetPanelNormal(sub_model_part, centered_points, process_info)
        self._CheckCurvature(sub_model_part, z_axis)

        material_axis = self._GetMaterialAxis(sub_model_part, z_axis, process_info)
        x_axis = self._GetEdgeDirectionClosestTo(sub_model_part, material_axis, z_axis)
        y_axis = np.cross(z_axis, x_axis)

        a, b = self._GetPanelDimensions(centered_points, x_axis, y_axis)
        if b > a * (1.0 + 1e-6):
            # convention a >= b: rotate the in-plane axes by 90 deg about z_axis
            x_axis, y_axis = -y_axis, x_axis
            a, b = b, a

        material_axis_offset = np.arctan2(np.dot(np.cross(x_axis, material_axis), z_axis),
                                          np.dot(x_axis, material_axis))
        self._CheckRectangular(sub_model_part, a, b)

        KratosMultiphysics.Logger.PrintInfo("Panel",
            f"'{sub_model_part.Name}': a={a:.4g}, b={b:.4g}, x={np.round(x_axis, 4)}, "
            f"z={np.round(z_axis, 4)}, material axis offset={np.degrees(material_axis_offset):.2f} deg")

        return PanelGeometry(x_axis, y_axis, z_axis, a, b, material_axis_offset)

    @staticmethod
    def _GetCenteredPoints(sub_model_part) -> np.ndarray:
        points = np.array([[node.X0, node.Y0, node.Z0] for node in sub_model_part.Nodes])
        return points - points.mean(axis=0)

    @staticmethod
    def _GetPanelNormal(sub_model_part, centered_points, process_info) -> np.ndarray:
        _, _, vh = np.linalg.svd(centered_points, full_matrices=False)
        z_axis = vh[-1]
        element = next(iter(sub_model_part.Elements))
        element_normal = np.array(element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_3, process_info)[0])
        return z_axis if np.dot(z_axis, element_normal) > 0.0 else -z_axis

    def _CheckCurvature(self, sub_model_part, z_axis) -> None:
        normals = np.array([self._GetElementNormal(element) for element in sub_model_part.Elements])
        max_angle = np.degrees(np.arccos(np.clip(np.abs(normals @ z_axis), 0.0, 1.0)).max())
        if max_angle > self.curvature_tolerance_degrees:
            KratosMultiphysics.Logger.PrintWarning("Panel",
                f"'{sub_model_part.Name}' appears to be curved: element normals deviate up to {max_angle:.2f} deg.")

    @staticmethod
    def _GetElementNormal(element) -> np.ndarray:
        geometry = element.GetGeometry()
        normals = np.array([geometry.UnitNormal(i) for i in range(geometry.IntegrationPointsNumber())])
        normal = normals.mean(axis=0)
        return normal / np.linalg.norm(normal)

    def _GetMaterialAxis(self, sub_model_part, z_axis, process_info) -> np.ndarray:
        axes = np.array([element.CalculateOnIntegrationPoints(SMA.LOCAL_MATERIAL_AXIS_1, process_info)[0]
                         for element in sub_model_part.Elements])
        axes -= np.outer(axes @ z_axis, z_axis)
        axes /= np.linalg.norm(axes, axis=1)[:, None]
        if np.min(axes @ axes[0]) < np.cos(np.radians(self.material_axis_tolerance_degrees)):
            raise RuntimeError(f"Panel '{sub_model_part.Name}': elements have different material axes, "
                               f"ply angles would refer to different directions.")
        return axes[0]

    @staticmethod
    def _GetEdgeDirectionClosestTo(sub_model_part, reference, z_axis) -> np.ndarray:
        """Boundary edges are used by one element only; returns the one most parallel to the reference."""
        edge_count = Counter()
        for element in sub_model_part.Elements:
            ids = [node.Id for node in element.GetGeometry()]
            for i in range(len(ids)):
                edge_count[frozenset((ids[i], ids[(i + 1) % len(ids)]))] += 1
        boundary_edges = [tuple(edge) for edge, count in edge_count.items() if count == 1]

        directions = []
        for i, j in boundary_edges:
            node_i, node_j = sub_model_part.GetNode(i), sub_model_part.GetNode(j)
            d = np.array([node_j.X0 - node_i.X0, node_j.Y0 - node_i.Y0, node_j.Z0 - node_i.Z0])
            d -= np.dot(d, z_axis) * z_axis
            directions.append(d / np.linalg.norm(d))
        directions = np.array(directions)

        x_axis = directions[np.argmax(np.abs(directions @ reference))]
        return x_axis if np.dot(x_axis, reference) > 0.0 else -x_axis

    @staticmethod
    def _GetPanelDimensions(centered_points, x_axis, y_axis) -> tuple[float, float]:
        return float(np.ptp(centered_points @ x_axis)), float(np.ptp(centered_points @ y_axis))

    def _CheckRectangular(self, sub_model_part, a, b) -> None:
        area = sum(element.GetGeometry().Area() for element in sub_model_part.Elements)
        if abs(area - a * b) > self.rectangle_tolerance * a * b:
            raise RuntimeError(f"Panel '{sub_model_part.Name}' is not rectangular in the panel coordinate system: "
                               f"element area {area:.4g} vs a*b {a * b:.4g}.")
