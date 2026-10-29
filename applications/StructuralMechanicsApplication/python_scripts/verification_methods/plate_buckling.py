import numpy as np
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.method_base import VerificationMethod
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.plate_buckling_solutions import (NavierSpeciallyOrthotropic,
                                                                                                       RitzSSSS)

class PlateBuckling(VerificationMethod):
    """Buckling of flat, rectangular plates simply supported on all edges (SSSS).

    Navier (exact) for specially orthotropic plates under constant N_xx, N_yy; Ritz otherwise
    (shear, bending-twisting coupling D16/D26, linear in-plane load gradients).
    The solution used is stored in metadata["solution"].

    Works on any component that provides the plate interface:
        boundary_conditions   "SSSS"
        geometry              a, b: dimensions along the component x and y axes
        laminate              D (3x3 bending stiffness), IsSymmetric()
        loads                 N_linear (3x3): rows c0, c1 (along x), c2 (along y); columns N_xx, N_yy, N_xy;
                              Kratos sign convention (tension positive)
    Load gradients not covered by the solutions (N_xx along x, N_yy along y, N_xy) are ignored.
    """
    name = "plate_buckling"
    category = "stability"
    plate_interface = ("boundary_conditions", "geometry", "laminate", "loads")
    ritz_terms = (10, 20, 30, 40, 50)   # Ritz converges from above: more terms -> lower, safer RF
    ritz_tolerance = 5e-3               # stop when the RF changes by less than 0.5 % between refinements
    relative_tolerance = 1e-6           # below this, N_xy and D16/D26 are treated as numerical noise
    gradient_tolerance = 1e-3           # below this, load gradients are treated as numerical noise

    def CheckCompatibility(self, component) -> None:
        missing = [attribute for attribute in self.plate_interface if not hasattr(component, attribute)]
        if missing:
            raise ValueError(f"{self.name}: component '{component.name}' does not provide {missing}.")
        if component.boundary_conditions != "SSSS":
            raise ValueError(f"{self.name}: boundary conditions '{component.boundary_conditions}' are not supported, only SSSS.")

    def IsApplicable(self, component) -> bool:
        return True

    def Evaluate(self, component) -> AnalysisResult:
        if not component.laminate.IsSymmetric():
            raise RuntimeError(f"{self.name}: laminate of '{component.name}' is not symmetric (B != 0).")

        D = component.laminate.D
        a, b = component.geometry.a, component.geometry.b
        coefficients = -component.loads.N_linear        # solutions use compression positive
        N0x, N0y, N0xy = coefficients[0]
        dN0x = coefficients[2, 0]                       # N_xx varying across the width (in-plane bending)
        dN0y = coefficients[1, 1]                       # N_yy varying along the length

        load_scale = max(abs(N0x), abs(N0y), abs(N0xy), abs(dN0x), abs(dN0y))
        has_shear = abs(N0xy) > self.relative_tolerance * load_scale
        has_gradient = max(abs(dN0x), abs(dN0y)) > self.gradient_tolerance * load_scale
        has_coupling = max(abs(D[0, 2]), abs(D[1, 2])) > self.relative_tolerance * D[0, 0]

        if not (has_shear or has_gradient or has_coupling):
            rf, m, n = NavierSpeciallyOrthotropic(D, a, b, N0x, N0y)
            metadata = {"solution": "navier", "m": m, "n": n}
        else:
            rf, W, ritz_info = self._RitzConverged(D, a, b, N0x, N0y, N0xy, dN0x, dN0y)
            metadata = {"solution": "ritz", **ritz_info}
            if W is not None:
                m, n = np.unravel_index(np.argmax(np.abs(W)), W.shape)
                metadata.update({"m": int(m) + 1, "n": int(n) + 1,
                                 "mode_share": float(W[m, n]**2 / np.sum(W**2)), "coefficients": W})

        log_metadata = {key: value for key, value in metadata.items() if key != "coefficients"}
        KratosMultiphysics.Logger.PrintInfo("PlateBuckling", f"'{component.name}': RF={rf:.4g}, {log_metadata}")
        return AnalysisResult(self.name, self.category, float(rf), metadata=metadata)

    def _RitzConverged(self, D, a, b, N0x, N0y, N0xy, dN0x, dN0y) -> tuple[float, np.ndarray | None, dict]:
        previous, _ = RitzSSSS(D, a, b, N0x, N0y, N0xy, self.ritz_terms[0], dN0x=dN0x, dN0y=dN0y)
        for n_terms in self.ritz_terms[1:]:
            current, W = RitzSSSS(D, a, b, N0x, N0y, N0xy, n_terms, dN0x=dN0x, dN0y=dN0y)
            change = abs(previous - current) / current if np.isfinite(current) else 0.0
            if change <= self.ritz_tolerance:
                return current, W, {"n_terms": n_terms, "relative_change": change, "converged": True}
            previous = current
        KratosMultiphysics.Logger.PrintWarning("PlateBuckling",
            f"Ritz not converged with {n_terms} terms (last change {change:.2%}).")
        return current, W, {"n_terms": n_terms, "relative_change": change, "converged": False}
