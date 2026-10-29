import numpy as np
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.method_base import VerificationMethod
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.plate_buckling_solutions import (NavierSpeciallyOrthotropic,
                                                                                                       RitzSSSS)

class PanelBuckling(VerificationMethod):
    """Buckling of flat rectangular panels under the panel-averaged running loads.

    Navier (exact) for specially orthotropic plates without shear, Ritz otherwise.
    """
    name = "panel_buckling"
    category = "stability"
    ritz_terms = (10, 20, 30, 40, 50)   # Ritz converges from above: more terms -> lower, safer RF
    ritz_tolerance = 5e-3               # stop when the RF changes by less than 0.5 % between refinements
    relative_tolerance = 1e-6   # below this, N_xy and D16/D26 are treated as numerical noise

    def CheckCompatibility(self, panel) -> None:
        if panel.boundary_conditions != "SSSS":
            raise ValueError(f"{self.name}: boundary conditions '{panel.boundary_conditions}' are not supported, only SSSS.")

    def IsApplicable(self, panel) -> bool:
        return True

    def Evaluate(self, panel) -> AnalysisResult:
        if not panel.laminate.IsSymmetric():
            raise RuntimeError(f"{self.name}: laminate of '{panel.name}' is not symmetric (B != 0).")

        D = panel.laminate.D
        a, b = panel.geometry.a, panel.geometry.b
        N0x, N0y, N0xy = -panel.loads.N_average       # Kratos: tension positive, solutions: compression positive

        has_shear = abs(N0xy) > self.relative_tolerance * max(abs(N0x), abs(N0y), abs(N0xy))
        has_coupling = max(abs(D[0, 2]), abs(D[1, 2])) > self.relative_tolerance * D[0, 0]

        if not has_shear and not has_coupling:
            rf, m, n = NavierSpeciallyOrthotropic(D, a, b, N0x, N0y)
            metadata = {"solution": "navier", "m": m, "n": n}
        else:
            rf, W, ritz_info = self._RitzConverged(D, a, b, N0x, N0y, N0xy)
            metadata = {"solution": "ritz", **ritz_info}
            if W is not None:
                m, n = np.unravel_index(np.argmax(np.abs(W)), W.shape)
                metadata.update({"m": int(m) + 1, "n": int(n) + 1,
                                 "mode_share": float(W[m, n]**2 / np.sum(W**2)), "coefficients": W})

        log_metadata = {key: value for key, value in metadata.items() if key != "coefficients"}
        KratosMultiphysics.Logger.PrintInfo("PanelBuckling", f"'{panel.name}': RF={rf:.4g}, {log_metadata}")
        return AnalysisResult(self.name, self.category, rf, metadata=metadata)

    def _RitzConverged(self, D, a, b, N0x, N0y, N0xy) -> tuple[float, np.ndarray | None, dict]:
        previous, _ = RitzSSSS(D, a, b, N0x, N0y, N0xy, self.ritz_terms[0])
        for n_terms in self.ritz_terms[1:]:
            current, W = RitzSSSS(D, a, b, N0x, N0y, N0xy, n_terms)
            change = abs(previous - current) / current if np.isfinite(current) else 0.0
            if change <= self.ritz_tolerance:
                return current, W, {"n_terms": n_terms, "relative_change": change, "converged": True}
            previous = current
        KratosMultiphysics.Logger.PrintWarning("PanelBuckling",
            f"Ritz not converged with {n_terms} terms (last change {change:.2%}).")
        return current, W, {"n_terms": n_terms, "relative_change": change, "converged": False}
