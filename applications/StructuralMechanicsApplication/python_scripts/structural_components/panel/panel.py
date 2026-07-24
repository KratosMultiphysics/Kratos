import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.structural_component import StructuralComponent
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.panel_buckling import PanelUniaxialBuckling, PanelBiaxialBuckling
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import (PanelGeometry, 
                                                                                                      PanelMaterial,
                                                                                                      PanelResponse,
                                                                                                      PanelLoadState)
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_geometry_interpreter import PanelGeometryInterpreter

class Panel(StructuralComponent):

    def __init__(self, 
                 sub_model_part, 
                 boundary_conditions: list[float]):
        super().__init__(sub_model_part, boundary_conditions)
        self.sub_model_part = sub_model_part
        self.geometry = None
        self.material = None
        self.response = None
        self.load_state = None
        self.analysis_methods = self._CreateAnalysisMethods()
    
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

    def PrepareAnalysis(self) -> None:
        self.ExtractMaterial()
        self.ExtractResponse()
        self.ClassifyLoadState()

    def RunAnalysis(self) -> None:
        self.RunMethods()
        self.StoreResults()
    
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
        self.geometry = PanelGeometryInterpreter().Interpret(self.sub_model_part)

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

    def RunMethods(self) -> None:
        self._RequireLoadState()
        self.analysis_results = []
        for method in self.analysis_methods:
            if method.IsApplicable(self):
                result = method.Evaluate(self)
                self.analysis_results.append(result)

                KratosMultiphysics.Logger.PrintInfo(
                "Panel",
                f"{result.method_name}: RF={result.value:.6e}"
            )
                
    def StoreResults(self) -> None:
        for result in self.analysis_results:
            if result.output_variable is None:
                continue

            self.sub_model_part.SetValue(result.output_variable, result.value)

            KratosMultiphysics.Logger.PrintInfo(
                "Panel",
                f"Stored/read back {result.method_name}: "
                f"{result.output_variable.Name()}={self.sub_model_part.GetValue(result.output_variable):.6e} "
                f"on '{self.sub_model_part.Name}'"
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
