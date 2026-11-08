"""Generator of ``PenaltyMethodFrictionlessMortarContactCondition`` (thesis section 4.3.3, penalty variant).

Run ``python3 generate_penalty_frictionless_mortar_condition.py [--combinations 2,2,2 ...] [--normal-variation false|true|both]``.
The documented version of this generator is the notebook ``penalty_frictionless_mortar_condition.ipynb`` in this folder.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import mortar_condition_generator as generator


def penalty_frictionless_functional(s, node, branch):
    """Galerkin functional of one slave node (thesis eqs. 4.57a, 4.71a). ``branch`` is ``inactive`` or ``active``."""
    rv_galerkin = 0
    if branch == "active":
        augmented_contact_pressure = (s.PenaltyParameter[node] * s.NormalGap[node])
        rv_galerkin += s.DynamicFactor[node] * (augmented_contact_pressure * s.NormalSlave.row(node)).dot(s.Dw1Mw2.row(node))
    else:
        rv_galerkin += 0
    return rv_galerkin


if __name__ == "__main__":
    generator.Main(generator.PENALTY_FRICTIONLESS, penalty_frictionless_functional, os.path.dirname(os.path.abspath(__file__)))
