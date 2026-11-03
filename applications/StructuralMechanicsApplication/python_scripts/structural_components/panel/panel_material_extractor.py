import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.laminate import Laminate

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

class PanelMaterialExtractor:

    def ExtractLaminate(self, sub_model_part, material_type: str, material_axis_offset: float) -> Laminate:
        """Laminate (CLT) of a shell component, ply angles in the component frame.

        Args:
            sub_model_part: Shell elements of the component; all must share one Properties.
            material_type: "isotropic" or "composite" (user input, checked against the Properties).
            material_axis_offset: Angle [rad] from the component x_axis to LOCAL_MATERIAL_AXIS_1.
        """
        properties = self._GetUniqueProperties(sub_model_part)
        has_layers = properties.Has(SMA.SHELL_ORTHOTROPIC_LAYERS)

        if properties.Has(SMA.SHELL_OFFSET) and properties[SMA.SHELL_OFFSET] != 0.0:
            raise RuntimeError(f"'{sub_model_part.Name}': SHELL_OFFSET is not supported, CLT is evaluated about the mid-plane.")

        if material_type == "composite":
            if not has_layers:
                raise RuntimeError(f"'{sub_model_part.Name}': material_type 'composite' requires SHELL_ORTHOTROPIC_LAYERS.")
            matrix = properties[SMA.SHELL_ORTHOTROPIC_LAYERS]
            layers = np.array([[matrix[i, j] for j in range(matrix.Size2())] for i in range(matrix.Size1())])
            laminate = Laminate.FromOrthotropicLayers(layers, np.degrees(material_axis_offset))
        else:
            if has_layers:
                raise RuntimeError(f"'{sub_model_part.Name}': material_type 'isotropic' but SHELL_ORTHOTROPIC_LAYERS is defined.")
            laminate = Laminate.FromIsotropic(properties[KratosMultiphysics.YOUNG_MODULUS],
                                              properties[KratosMultiphysics.POISSON_RATIO],
                                              properties[KratosMultiphysics.THICKNESS])

        KratosMultiphysics.Logger.PrintInfo("Panel",
            f"'{sub_model_part.Name}': t={laminate.thickness:.4g}, symmetric={laminate.IsSymmetric()}")
        return laminate

    @staticmethod
    def _GetUniqueProperties(sub_model_part):
        """A component is one laminate: all elements must share the same Properties."""
        properties = {element.Properties.Id: element.Properties for element in sub_model_part.Elements}
        if len(properties) != 1:
            raise RuntimeError(f"'{sub_model_part.Name}' has {len(properties)} different Properties, expected one laminate.")
        return next(iter(properties.values()))
