"""Generator of ``AugmentedLagrangianMethodFrictionalMortarContactCondition`` (thesis section 4.3.4).

Run ``python3 generate_frictional_mortar_condition.py [--combinations 2,2,2 ...] [--normal-variation false|true|both]``.
The documented version of this generator is the notebook ``ALM_frictional_mortar_condition.ipynb`` in this folder.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import mortar_condition_generator as generator


def frictional_functional(s, node, branch):
    """Galerkin functional of one slave node in one active-set branch (thesis eqs. 4.60-4.62, 4.72).

    ``branch`` is one of ``inactive``, ``slip_objective``, ``slip_non_objective``, ``stick_objective``,
    ``stick_non_objective``. The multiplier, penalty and slip quantities are those of ``SymbolSet``.
    """
    rv_galerkin = 0
    if branch == "inactive":
        # The multiplier is penalised to zero (gap zone): normal and tangential components
        rv_galerkin -= s.ScaleFactor**2 / s.PenaltyParameter[node] * s.LMNormal[node] * s.wLMNormal[node]
        rv_galerkin -= s.ScaleFactor**2 / (s.PenaltyParameter[node] * s.TangentFactor) * (s.LMTangent.row(node)).dot(s.wLMTangent.row(node))
    else:
        # Normal contact constraint weighted with the normal test multiplier
        rv_galerkin += s.ScaleFactor * s.NormalGap[node] * s.wLMNormal[node]

        augmented_normal_contact_pressure = (s.ScaleFactor * s.LMNormal[node] + s.PenaltyParameter[node] * s.NormalGap[node])
        augmented_lm = (s.ScaleFactor * s.LM.row(node) + s.PenaltyParameter[node] * s.NormalGap[node] * s.NormalSlave.row(node))

        if branch in ("slip_objective", "slip_non_objective"):
            # Coulomb limit: the tangential multiplier is driven to -mu * p_n * tau
            augmented_tangent_contact_pressure = - s.mu[node] * augmented_normal_contact_pressure * s.TangentSlave.row(node)
            rv_galerkin -= s.ScaleFactor**2 / s.PenaltyParameter[node] * (s.wLMTangent.row(node).dot(s.LMTangent.row(node) - augmented_tangent_contact_pressure / s.ScaleFactor))
        else:
            # Stick: the tangential slip is constrained to zero and augments the multiplier
            tangent_slip = s.TangentSlipObjective if branch == "stick_objective" else s.TangentSlipNonObjective
            augmented_lm += s.TangentFactor * s.PenaltyParameter[node] * tangent_slip.row(node)
            rv_galerkin += s.ScaleFactor * (tangent_slip.row(node)).dot(s.wLMTangent.row(node))

        # Virtual work of the augmented contact traction
        rv_galerkin += s.DynamicFactor[node] * (augmented_lm).dot(s.Dw1Mw2.row(node))

    return rv_galerkin


if __name__ == "__main__":
    generator.Main(generator.ALM_FRICTIONAL, frictional_functional, os.path.dirname(os.path.abspath(__file__)))
