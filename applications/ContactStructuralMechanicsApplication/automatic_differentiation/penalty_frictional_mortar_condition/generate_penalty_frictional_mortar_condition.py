"""Generator of ``PenaltyMethodFrictionalMortarContactCondition`` (thesis section 4.3.4, penalty variant).

Run ``python3 generate_penalty_frictional_mortar_condition.py [--combinations 2,2,2 ...] [--normal-variation false|true|both]``.
The documented version of this generator is the notebook ``penalty_frictional_mortar_condition.ipynb`` in this folder.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import mortar_condition_generator as generator


def penalty_frictional_functional(s, node, branch):
    """Galerkin functional of one slave node in one active-set branch (thesis eqs. 4.57-4.58, 4.71).

    ``branch`` is one of ``inactive``, ``slip``, ``stick``. Only the objective slip is used (there is no
    multiplier equation, so the objective/non-objective switch of the ALM condition is not needed).
    """
    rv_galerkin = 0
    if branch == "inactive":
        rv_galerkin += 0
    else:
        # Normal penalty traction
        augmented_normal_contact_pressure = s.PenaltyParameter[node] * s.NormalGap[node]
        rv_galerkin += s.DynamicFactor[node] * augmented_normal_contact_pressure * s.NormalwGap[node]

        if branch == "slip":
            # Coulomb traction -mu * p_n * tau
            augmented_tangent_contact_pressure = - s.mu[node] * augmented_normal_contact_pressure * s.TangentSlave.row(node)
        else:
            # Tangential penalty traction
            augmented_tangent_contact_pressure = s.TangentFactor * s.PenaltyParameter[node] * s.TangentSlipObjective.row(node)
        rv_galerkin += s.DynamicFactor[node] * augmented_tangent_contact_pressure.dot(s.TangentwSlipObjective.row(node))

    return rv_galerkin


if __name__ == "__main__":
    generator.Main(generator.PENALTY_FRICTIONAL, penalty_frictional_functional, os.path.dirname(os.path.abspath(__file__)))
