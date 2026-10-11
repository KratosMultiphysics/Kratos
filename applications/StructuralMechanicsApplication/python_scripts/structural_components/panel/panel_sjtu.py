import numpy as np
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.structural_component import StructuralComponent
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_sjtu_geometry_interpreter import PanelSJTUGeometryInterpreter
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.laminate import Laminate
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelLoads

class PanelSJTU(StructuralComponent):
    """Panel component for isotropic and composite shells.

    The panel provides data (geometry, laminate, running loads in the panel frame).
    Handbook methods decide how to evaluate it.
    """

    def __init__(self, sub_model_part, boundary_conditions: str, material_type: str, metadata):
        super().__init__(sub_model_part, metadata)
        self.material_type = material_type
        self.boundary_conditions = boundary_conditions
        self.geometry = None
        self.laminate = None
        self.loads = None
        self.analysis_methods = []
        self.analysis_results = {}
    
    @classmethod
    def FromKratosParametersObject(cls, sub_model_part, data):
        metadata = data["metadata"] if data.Has("metadata") else None
        return cls(sub_model_part,
                   data["boundary_conditions"].GetString(),
                   data["material_type"].GetString(),
                   metadata)

    def Initialize(self):
        self.ExtractGeometry()

    def PrepareAnalysis(self):
        self.ExtractMaterial()
        self.ExtractResponse()

    def ExtractGeometry(self) -> None:
        self.geometry = PanelSJTUGeometryInterpreter().Interpret(self.sub_model_part)

    def ExtractMaterial(self) -> None:
        properties = self._GetPanelProperties()
        has_layers = properties.Has(SMA.SHELL_ORTHOTROPIC_LAYERS)

        if properties.Has(SMA.SHELL_OFFSET) and properties[SMA.SHELL_OFFSET] != 0.0:
            raise RuntimeError(f"Panel '{self.name}': SHELL_OFFSET is not supported, CLT is evaluated about the mid-plane.")

        if self.material_type == "composite":
            if not has_layers:
                raise RuntimeError(f"Panel '{self.name}': material_type 'composite' requires SHELL_ORTHOTROPIC_LAYERS.")
            matrix = properties[SMA.SHELL_ORTHOTROPIC_LAYERS]
            layers = np.array([[matrix[i, j] for j in range(matrix.Size2())] for i in range(matrix.Size1())])
            self.laminate = Laminate.FromOrthotropicLayers(layers, np.degrees(self.geometry.material_axis_offset))

        else:
            if has_layers:
                raise RuntimeError(f"Panel '{self.name}': material_type 'isotropic' but SHELL_ORTHOTROPIC_LAYERS is defined.")
            self.laminate = Laminate.FromIsotropic(properties[KratosMultiphysics.YOUNG_MODULUS],
                                                properties[KratosMultiphysics.POISSON_RATIO],
                                                properties[KratosMultiphysics.THICKNESS])

        KratosMultiphysics.Logger.PrintInfo("PanelSJTU",
            f"'{self.name}': t={self.laminate.thickness:.4g}, symmetric={self.laminate.IsSymmetric()}, "
            f"D=\n{np.array2string(self.laminate.D, precision=4)}")

    def ExtractResponse(self):
        R = self.geometry.rotation_matrix
        process_info = self.sub_model_part.ProcessInfo
        elements = list(self.sub_model_part.Elements)

        N = np.array([self._InPlaneComponents(self._AverageOverIntegrationPoints(element, SMA.SHELL_FORCE_GLOBAL, process_info), R)
                      for element in elements])
        M = np.array([self._InPlaneComponents(self._AverageOverIntegrationPoints(element, SMA.SHELL_MOMENT_GLOBAL, process_info), R)
                      for element in elements])

        self.loads = PanelLoads(element_ids=np.array([element.Id for element in elements]),
                                element_areas=np.array([element.GetGeometry().Area() for element in elements]),
                                N=N, M=M)

        KratosMultiphysics.Logger.PrintInfo("PanelSJTU",
            f"'{self.name}': N_average (xx, yy, xy) = {np.array2string(self.loads.N_average, precision=4)}")

        print(np.array(elements[0].CalculateOnIntegrationPoints(SMA.SHELL_FORCE, process_info)[0]))
        print(elements[0].CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, process_info)[0])
        print(np.abs(self.loads.M).max())

    @staticmethod
    def _AverageOverIntegrationPoints(element, variable, process_info) -> np.ndarray:
        values = element.CalculateOnIntegrationPoints(variable, process_info)
        return np.mean([[[value[i, j] for j in range(3)] for i in range(3)] for value in values], axis=0)

    @staticmethod
    def _InPlaneComponents(tensor_global: np.ndarray, R: np.ndarray) -> np.ndarray:
        """Rotate a 3x3 tensor from global into the panel frame and return (xx, yy, xy)."""
        tensor_panel = R @ tensor_global @ R.T
        return np.array([tensor_panel[0, 0], tensor_panel[1, 1], tensor_panel[0, 1]])

    def _GetPanelProperties(self):
        """A Panel is one laminate: all elements must share the same Properties.
        """
        properties = {element.Properties.Id: element.Properties for element in self.sub_model_part.Elements}
        if len(properties) != 1:
            raise RuntimeError(f"Panel '{self.name}' has {len(properties)} different Properties, expected one laminate.")
        return next(iter(properties.values()))