from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelLoads
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

class PanelResponseExtractor:

    def ExtractRunningLoads(self, sub_model_part, geometry) -> PanelLoads:
        """Running loads N, M per element in the component frame (tension positive).

        Args:
            sub_model_part: Shell elements of the component.
            geometry: PanelGeometry with the component axes and dimensions a, b.
        """
        R = geometry.rotation_matrix
        process_info = sub_model_part.ProcessInfo
        elements = list(sub_model_part.Elements)

        N = np.array([self._InPlaneComponents(self._AverageOverIntegrationPoints(element, SMA.SHELL_FORCE_GLOBAL, process_info), R)
                      for element in elements])
        M = np.array([self._InPlaneComponents(self._AverageOverIntegrationPoints(element, SMA.SHELL_MOMENT_GLOBAL, process_info), R)
                      for element in elements])

        return PanelLoads(element_ids=np.array([element.Id for element in elements]),
                          element_areas=np.array([element.GetGeometry().Area() for element in elements]),
                          centroids=self._NormalizedCentroids(sub_model_part, elements, geometry),
                          N=N, M=M)

    @staticmethod
    def _AverageOverIntegrationPoints(element, variable, process_info) -> np.ndarray:
        values = element.CalculateOnIntegrationPoints(variable, process_info)
        return np.mean([[[value[i, j] for j in range(3)] for i in range(3)] for value in values], axis=0)

    @staticmethod
    def _InPlaneComponents(tensor_global: np.ndarray, R: np.ndarray) -> np.ndarray:
        """Rotate a 3x3 tensor from global into the component frame and return (xx, yy, xy)."""
        tensor_local = R @ tensor_global @ R.T
        return np.array([tensor_local[0, 0], tensor_local[1, 1], tensor_local[0, 1]])

    @staticmethod
    def _NormalizedCentroids(sub_model_part, elements, geometry) -> np.ndarray:
        """Element centroids as (xi, eta) in [0, 1], measured from the panel corner in the panel frame."""
        R = geometry.rotation_matrix
        nodes = np.array([[node.X0, node.Y0, node.Z0] for node in sub_model_part.Nodes]) @ R.T
        corner = nodes.min(axis=0)[:2]
        centroids = np.array([np.mean([[node.X0, node.Y0, node.Z0] for node in element.GetGeometry()], axis=0)
                              for element in elements]) @ R.T
        return (centroids[:, :2] - corner) / np.array([geometry.a, geometry.b])
