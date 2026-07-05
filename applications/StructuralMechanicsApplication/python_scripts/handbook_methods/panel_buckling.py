import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.method_base import HandbookMethod
import KratosMultiphysics.StructuralMechanicsApplication as SMA

class PanelUniaxialBuckling(HandbookMethod):
    name = "panel_uniaxial_buckling"
    category = "stability"

    def IsApplicable(self, panel) -> bool:
        panel._RequireLoadState()
        return panel.load_state.is_uniaxial_compression

    def Evaluate(self, panel) -> AnalysisResult:

        if panel.load_state.has_x_compression:
            applied_stress = abs(panel.response.sigma_xx)
            buckling_length = panel.a
        else:
            applied_stress = abs(panel.response.sigma_yy)
            buckling_length = panel.b

        sigma_crit = (
            panel.E * np.pi**2
            / (12.0 * (1.0 - panel.nu**2))
            * (panel.t / buckling_length)**2
        )

        rf = sigma_crit / applied_stress

        return AnalysisResult(self.name, self.category, rf)

class PanelBiaxialBuckling(HandbookMethod):
    name = "panel_biaxial_buckling"
    category = "stability"

    def IsApplicable(self, panel) -> bool:
        panel._RequireLoadState()
        return panel.load_state.is_biaxial_compression

    def Evaluate(self, panel) -> AnalysisResult:
        sigma_x = abs(panel.response.sigma_xx)
        sigma_y = abs(panel.response.sigma_yy)

        if sigma_x <= 0.0:
            raise RuntimeError(
                f"Panel '{panel.sub_model_part.Name}' has no x-compression for biaxial buckling."
            )

        beta = sigma_y / sigma_x

        m = 1
        n = 1

        a = panel.a
        b = panel.b
        t = panel.t
        E = panel.E
        nu = panel.nu

        D = E * t**3 / (12.0 * (1.0 - nu**2))

        sigma_x_crit = (
            D * np.pi**2 * ((m / a)**2 + (n / b)**2)**2
            / (t * ((m / a)**2 + beta * (n / b)**2))
        )

        rf = sigma_x_crit / sigma_x

        return AnalysisResult(self.name, self.category, rf, output_variable=SMA.RESPONSE_VALUE)
