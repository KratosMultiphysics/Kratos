from KratosMultiphysics.StructuralMechanicsApplication.structural_components.structural_component import StructuralComponent
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_geometry_interpreter import PanelGeometryInterpreter
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_material_extractor import PanelMaterialExtractor
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_response_extractor import PanelResponseExtractor



class Panel(StructuralComponent):
    """Flat, rectangular panel of shell elements.

    Provides geometry, laminate and running loads in the panel coordinate system.
    The data is extracted by dedicated classes; handbook methods evaluate it.
    """

    def __init__(self, sub_model_part, boundary_conditions: str, material_type: str, metadata):
        super().__init__(sub_model_part, metadata)
        self.boundary_conditions = boundary_conditions
        self.material_type = material_type
        self.geometry = None
        self.laminate = None
        self.loads = None

    @classmethod
    def FromKratosParametersObject(cls, sub_model_part, data):
        metadata = data["metadata"] if data.Has("metadata") else None
        return cls(sub_model_part,
                   data["boundary_conditions"].GetString(),
                   data["material_type"].GetString(),
                   metadata)

    def Initialize(self) -> None:
        """Once per analysis: the mesh does not change."""
        self.geometry = PanelGeometryInterpreter().Interpret(self.sub_model_part)

    def PrepareAnalysis(self) -> None:
        """After the solution: properties may change between evaluations, loads change per load case."""
        self.laminate = PanelMaterialExtractor().ExtractLaminate(self.sub_model_part, 
                                                                 self.material_type,
                                                                 self.geometry.material_axis_offset)
        
        self.loads = PanelResponseExtractor().ExtractRunningLoads(self.sub_model_part, 
                                                                  self.geometry)
        print(self.loads.N_average)
        print(self.loads.N_linear)

