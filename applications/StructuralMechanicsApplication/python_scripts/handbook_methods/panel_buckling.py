import numpy as np
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.method_base import HandbookMethod
import KratosMultiphysics.StructuralMechanicsApplication as SMA

class PanelUniaxialBuckling(HandbookMethod):
    """Currently limited to isotropic simply supported panels.

    Returns:
        AnalysisResult: Constains method name, category and the reserve factor
    """
    name = "panel_uniaxial_buckling"
    category = "stability"

    def IsApplicable(self, panel) -> bool:
        panel._RequireLoadState()
        return panel.load_state.is_uniaxial_compression

    def Evaluate(self, panel) -> AnalysisResult:
        KratosMultiphysics.Logger.PrintInfo("::[PanelUniaxialBuckling]::", "Evaluate Critical Buckling Load")
        if panel.load_state.has_x_compression:
            applied_stress = abs(panel.response.sigma_xx)
            buckling_length = panel.a
            width = panel.b
        else:
            applied_stress = abs(panel.response.sigma_yy)
            buckling_length = panel.b
            width = panel.a

        #TODO: Move this to panel material
        D = panel.E * panel.t**3 / (12.0 * (1.0 - panel.nu**2))

        D11 = D
        D22 = D
        D12 = panel.nu * D
        D66 = 0.5 * (1.0 - panel.nu) * D

        #Taken from Mittelstedt for simply supported plates
        m = (buckling_length/width) * np.pow((D22/D11), 1/4)
        m_floor = max(1, int(np.floor(m)))
        m_ceil = max(1, int(np.ceil(m)))

        print(f"m: {m} \n m_floor: {m_floor} \n m_ceil: {m_ceil}")

        N_floor = np.pi**2 * (D11 * m_floor**2 / buckling_length**2 + 2.0 * (D12 + 2.0 * D66) / width**2 + D22 * buckling_length**2 / (width**4 * m_floor**2))
        N_ceil = np.pi**2 * (D11 * m_ceil**2 / buckling_length**2 + 2.0 * (D12 + 2.0 * D66) / width**2 + D22 * buckling_length**2 / (width**4 * m_ceil**2))

        print(f"N_floor: {N_floor}, \n N_ceil: {N_ceil}")

        N_crit = min(N_floor, N_ceil)

        sigma_crit = N_crit/panel.t

        rf = sigma_crit / applied_stress

        return AnalysisResult(self.name, self.category, rf, output_variable=SMA.RESPONSE_VALUE)

class PanelBiaxialBuckling(HandbookMethod):
    name = "panel_biaxial_buckling"
    category = "stability"

    def IsApplicable(self, panel) -> bool:
        panel._RequireLoadState()
        return panel.load_state.is_biaxial_compression

    def Evaluate(self, panel) -> AnalysisResult:
        KratosMultiphysics.Logger.PrintInfo("::[PanelBiaxialBuckling]::", "Evaluate Critical Buckling Load")
        sigma_x = abs(panel.response.sigma_xx)
        sigma_y = abs(panel.response.sigma_yy)

        if sigma_x <= 0.0:
            raise RuntimeError(
                f"Panel '{panel.sub_model_part.Name}' has no x-compression for biaxial buckling."
            )

        beta = sigma_y / sigma_x

        a = panel.a
        b = panel.b
        t = panel.t
        E = panel.E
        nu = panel.nu

        D = E * t**3 / (12.0 * (1.0 - nu**2))

        max_half_waves = 20

        mode_pairs = ([(1, n) for n in range(1, max_half_waves + 1)] + [(m, 1) for m in range(1, max_half_waves + 1)])

        sigma_x_crit, m_crit, n_crit = min(
            (
                D * np.pi**2 * ((m / a)**2 + (n / b)**2)**2/ (t * ((m / a)**2 + beta * (n / b)**2)), m, n,
                ) 
                for m, n in mode_pairs
                )

        rf = sigma_x_crit / sigma_x

        KratosMultiphysics.Logger.PrintInfo("::[PanelBiaxialBuckling]::", "m_crit:", m_crit, "n_crit:", n_crit)

        return AnalysisResult(self.name, 
                              self.category, 
                              rf, 
                              output_variable=SMA.RESPONSE_VALUE,
                              metadata={
                                  "m": m_crit,
                                  "n": n_crit
                              })
