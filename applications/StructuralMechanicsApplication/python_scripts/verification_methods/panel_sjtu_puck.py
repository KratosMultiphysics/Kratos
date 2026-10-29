import numpy as np
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.method_base import VerificationMethod
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.puck import PuckAnalysis

class PuckFirstPlyFailure(VerificationMethod):
    """Puck first-ply failure from CLT ply stresses, evaluated per element at the top and bottom face of each ply.

    RF = min(RF_FF, RF_IFF). No degradation, no last-ply failure.
    """
    name = "puck"
    category = "strength"
    default_inclinations = {"p_trtr_c": 0.25, "p_trpa_t": 0.25, "p_trpa_c": 0.25}
    _puck = PuckAnalysis()      # reuses the vectorized PuckFF / PuckIFF

    def CheckCompatibility(self, panel) -> None:
        if panel.material_type != "composite":
            raise ValueError(f"{self.name}: only applicable to composite panels, '{panel.name}' is {panel.material_type}.")

    def IsApplicable(self, panel) -> bool:
        return True

    def Evaluate(self, panel) -> AnalysisResult:
        laminate = panel.laminate
        stresses = laminate.PlyStresses(panel.loads.N, panel.loads.M)          # (n_elements, n_plies, 2, 3)
        sigma_1, sigma_2, tau_21 = stresses[..., 0], stresses[..., 1], stresses[..., 2]

        def Strength(name):
            return np.array([getattr(ply.strength, name) for ply in laminate.plies])[None, :, None]
        p = self._InclinationParameters(panel)

        with np.errstate(divide="ignore", invalid="ignore"):
            _, f_E_FF, mode_FF = self._puck.PuckFF(sigma_1, Strength("R_pa_t"), Strength("R_pa_c"))
            _, f_E_IFF, mode_IFF, theta_fp = self._puck.PuckIFF(sigma_2, tau_21, Strength("R_tr_t"), Strength("R_tr_c"),
                                                                Strength("R_trpa"), p["p_trtr_c"], p["p_trpa_t"], p["p_trpa_c"])

        rf_ff, rf_iff = self._ReserveFactor(f_E_FF.max()), self._ReserveFactor(f_E_IFF.max())
        governing_f_E, governing_mode = (f_E_FF, mode_FF) if rf_ff <= rf_iff else (f_E_IFF, mode_IFF)
        element, ply, side = np.unravel_index(np.argmax(governing_f_E), governing_f_E.shape)
        # worst element and face per ply, FF or IFF; mirrored plies combined (symmetric laminate -> one value per ply pair)
        f_E_per_ply = np.maximum(f_E_FF, f_E_IFF).max(axis=(0, 2))
        f_E_per_pair = np.maximum(f_E_per_ply, f_E_per_ply[::-1])[: (len(f_E_per_ply) + 1) // 2]
        rf_per_ply = [self._ReserveFactor(f_E) for f_E in f_E_per_pair]

        metadata = {"rf_ff": rf_ff, "rf_iff": rf_iff,
                    "rf_per_ply": rf_per_ply,
                    "failure_mode": str(governing_mode[element, ply, side]),
                    "element_id": int(panel.loads.element_ids[element]),
                    "ply": int(ply) + 1, "ply_angle": laminate.plies[ply].angle,
                    "side": ("top", "bottom")[side],
                    "theta_fp": float(theta_fp[element, ply, side]) if rf_iff < rf_ff else None,
                    }
        rf = min(rf_ff, rf_iff)

        KratosMultiphysics.Logger.PrintInfo("Puck", f"'{panel.name}': RF={rf:.4g}, {metadata}")
        return AnalysisResult(self.name, self.category, rf, metadata=metadata)

    @staticmethod
    def _ReserveFactor(f_E: float) -> float:
        """Puck stress exposures scale linearly with the load, so RF = 1 / f_E."""
        return float(1.0 / f_E) if f_E > 0.0 else np.inf

    def _InclinationParameters(self, panel) -> dict:
        p = dict(self.default_inclinations)
        if panel.metadata is not None and panel.metadata.Has("puck_parameters"):
            for key in p:
                if panel.metadata["puck_parameters"].Has(key):
                    p[key] = panel.metadata["puck_parameters"][key].GetDouble()
        return p
