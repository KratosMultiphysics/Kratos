#  ____             _       |
# |  _ \ _   _  ___| | __   |
# | |_) | | | |/ __| |/ /   |
# |  __/| |_| | (__|   <    |
# |_|    \__,_|\___|_|\_\   ANALYSIS
#
#  Main authors:    Lucas Rimpl
#  Co-authors:      Tobias Siemer
#

#The current test model using these material data does not produce a fiber failure in the “First Fiber Fracture” degradation analysis. 
#Therefore, this is not fully covered by the test

import dataclasses
import numpy as np
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.method_base import VerificationMethod
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.puck_criteria import (FibreFracture,
                                                                                               InterFibreFracture,
                                                                                               ReserveFactor)
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.laminate import Laminate


class PuckAnalysis(VerificationMethod):
    """Puck failure analysis of a laminate, with progressive degradation after inter-fibre fracture.

    Ply stresses from CLT (panel.laminate, panel.loads), evaluated per element at the top and bottom face of
    each ply. Results are flat arrays in the order element, ply (top to bottom), face (top, bottom).

    Degradation per element: plies with inter-fibre fracture in Mode A or B get reduced transverse and shear
    moduli (factors from metadata); the stresses are recalculated with the degraded laminate until no further
    ply fails. Degradation is irreversible. Mode C (risk of delamination) stops the degradation.
    """
    name = "puck"
    category = "strength"
    criteria = ("first_ply_failure", "first_fiber_failure", "last_ply_failure")
    default_inclinations = {"p_trtr_c": 0.25, "p_trpa_t": 0.25, "p_trpa_c": 0.25}
    default_degradation_factors = {"E_tr_A": 0.0, "G_patr_A": 0.0, "E_tr_B": 0.5, "G_patr_B": 0.5}
    singular_tolerance = 1e10

    def IsApplicable(self, structural_component) -> bool:
        return structural_component.material_type == "composite"

    def Evaluate(self, structural_component) -> AnalysisResult:
        criterion = structural_component.metadata["criterion"].GetString()
        if criterion not in self.criteria:
            raise ValueError(f"{self.name}: unknown failure criterion '{criterion}', expected one of {self.criteria}.")
        inclinations = self._MetadataValues(structural_component.metadata, "puck_parameters", self.default_inclinations)
        factors = self._MetadataValues(structural_component.metadata, "degradation_factors", self.default_degradation_factors)

        laminate, loads = structural_component.laminate, structural_component.loads
        stresses = laminate.PlyStresses(loads.N, loads.M)                       # (element, ply, face, 3)
        f_E_FF, mode_FF, f_E_IFF, mode_IFF, theta_fp = self._Criteria(laminate, stresses, inclinations)
        RF_FF, RF_IFF = ReserveFactor(f_E_FF), ReserveFactor(f_E_IFF)

        RF_FF_degraded = np.zeros_like(RF_FF)                                   # 0: no degradation step performed
        RF_IFF_degraded = np.zeros_like(RF_IFF)
        for element in range(len(loads.element_ids)):
            result = self._DegradeElement(laminate, loads.N[element], loads.M[element], RF_FF[element],
                                          RF_IFF[element], mode_IFF[element], criterion, inclinations, factors)
            if result is not None:
                RF_FF_degraded[element], RF_IFF_degraded[element] = result

        metadata = {"RF_FF": RF_FF.ravel(), "RF_IFF": RF_IFF.ravel(),
                    "Puck_DegradationFF": RF_FF_degraded.ravel(), "Puck_DegradationIFF": RF_IFF_degraded.ravel(),
                    "failure_mode_FF": mode_FF.ravel(), "failure_mode_IFF": mode_IFF.ravel(), "theta_fp": theta_fp.ravel()}
        return AnalysisResult(self.name, self.category, float(np.nanmin(RF_IFF)), metadata=metadata)

    @staticmethod
    def _Criteria(laminate, stresses, inclinations):
        """Puck criteria for stresses (..., ply, face, 3) of the given laminate."""
        def PerPly(name):
            return np.array([getattr(ply.strength, name) for ply in laminate.plies])[:, None]   # broadcasts to (ply, face)
        f_E_FF, mode_FF = FibreFracture(stresses[..., 0], PerPly("R_pa_t"), PerPly("R_pa_c"))
        f_E_IFF, mode_IFF, theta_fp = InterFibreFracture(stresses[..., 1], stresses[..., 2],
                                                         PerPly("R_tr_t"), PerPly("R_tr_c"), PerPly("R_trpa"),
                                                         inclinations["p_trtr_c"], inclinations["p_trpa_t"], inclinations["p_trpa_c"])
        return f_E_FF, mode_FF, f_E_IFF, mode_IFF, theta_fp

    def _DegradeElement(self, laminate, N, M, RF_FF, RF_IFF, mode_IFF, criterion, inclinations, factors):
        """Progressive degradation of one element. RF_* and mode_IFF: (ply, face).

        Returns (RF_FF, RF_IFF), each (ply, face), after the last degradation step, or None if no step was performed.
        """
        degradation = [None] * len(laminate.plies)          # per ply: None (intact), "A" or "B"
        result = None
        while True:
            if np.max(RF_IFF) < 1.0 or np.min(RF_IFF) >= 1.0:       # all faces failed, or nothing to degrade
                break
            if criterion == "first_fiber_failure" and np.any(RF_FF < 1.0):
                break
            if criterion == "first_ply_failure" and np.any(RF_IFF < 1.0):
                break

            critical_face = np.argmin(RF_IFF, axis=1)
            plies = range(len(laminate.plies))
            ply_RF = RF_IFF[plies, critical_face]
            ply_mode = mode_IFF[plies, critical_face]
            if np.any((ply_RF < 1.0) & (ply_mode == "IFF Mode C")):
                break                                               # delamination risk: not allowed

            new_degradation = [state if state is not None or rf >= 1.0 else mode[-1]
                               for state, rf, mode in zip(degradation, ply_RF, ply_mode)]
            if new_degradation == degradation:
                break                                               # no further ply failed
            degradation = new_degradation

            degraded = self._DegradedLaminate(laminate, degradation, factors)
            if not self._IsLoadCarrying(degraded, self.singular_tolerance):
                KratosMultiphysics.Logger.PrintWarning("Puck", "Degraded laminate has no stiffness left in some "
                                                       "deformation mode: total failure of the element.")
                failed = np.zeros_like(RF_IFF)
                return failed, failed
            stresses = degraded.PlyStresses(N[None, :], M[None, :])[0]
            f_E_FF, _, f_E_IFF, mode_IFF, _ = self._Criteria(degraded, stresses, inclinations)
            RF_FF, RF_IFF = ReserveFactor(f_E_FF), ReserveFactor(f_E_IFF)
            result = (RF_FF, RF_IFF)
        return result

    @staticmethod
    def _DegradedLaminate(laminate, degradation, factors) -> Laminate:
        plies = []
        for ply, state in zip(laminate.plies, degradation):
            if state is not None:
                ply = dataclasses.replace(ply, E2=factors[f"E_tr_{state}"] * ply.E2, G12=factors[f"G_patr_{state}"] * ply.G12)
            plies.append(ply)
        return Laminate(plies)

    @staticmethod
    def _MetadataValues(metadata, group: str, defaults: dict) -> dict:
        values = dict(defaults)
        if metadata is not None and metadata.Has(group):
            for key in values:
                if metadata[group].Has(key):
                    values[key] = metadata[group][key].GetDouble()
        return values

    @staticmethod
    def _IsLoadCarrying(laminate, tolerance: float) -> bool:
        """False if the laminate has (almost) no stiffness left in some deformation mode.

        A and D are normalized with the thickness (A / h, 12 D / h^3), so that the condition number
        compares stiffnesses of the same unit instead of membrane against bending stiffness.
        """
        h = laminate.thickness
        W = np.diag([h**-0.5] * 3 + [np.sqrt(12.0 / h**3)] * 3)
        return np.linalg.cond(W @ laminate.ABD @ W) < tolerance