"""Generator of ``AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition`` (thesis section 4.3.3, vector multiplier).

Run ``python3 generate_frictionless_components_mortar_condition.py [--combinations 2,2,2 ...] [--normal-variation false|true|both]``.
The documented version of this generator is the notebook ``ALM_frictionless_components_mortar_condition.ipynb`` in this folder.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import mortar_condition_generator as generator


def frictionless_components_functional(s, node, branch):
    """Galerkin functional of one slave node (thesis eqs. 4.36-4.44 with a vector multiplier whose
    tangential part is penalised away). ``branch`` is ``inactive`` or ``active``."""
    rv_galerkin = 0
    if branch == "active":
        augmented_lm = (s.ScaleFactor * s.LM.row(node) + s.PenaltyParameter[node] * s.NormalGap[node] * s.NormalSlave.row(node))
        rv_galerkin += s.DynamicFactor[node] * (augmented_lm).dot(s.Dw1Mw2.row(node))
        rv_galerkin += s.ScaleFactor * s.NormalGap[node] * s.wLMNormal[node]
        rv_galerkin -= s.ScaleFactor**2 / s.PenaltyParameter[node] * (s.wLMTangent.row(node).dot(s.LMTangent.row(node)))
    else:
        rv_galerkin -= s.ScaleFactor**2 / s.PenaltyParameter[node] * (s.wLM.row(node).dot(s.LM.row(node)))
    return rv_galerkin


if __name__ == "__main__":
    generator.Main(generator.ALM_FRICTIONLESS_COMPONENTS, frictionless_components_functional, os.path.dirname(os.path.abspath(__file__)))
