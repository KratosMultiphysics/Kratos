"""Guard test of the symbolic generation of the mortar contact conditions.

The condition local systems are generated with sympy by the scripts / notebooks of
``automatic_differentiation/`` (see ``mortar_condition_generator.py``). This test regenerates
the smallest geometry of the ALM frictionless family into a temporary folder and checks that the
pipeline still runs with the installed sympy version, that the emitted code has the expected
structure and no symbolic leftovers, and that it is numerically equivalent to the committed file.
"""
import os
import sys
import tempfile

import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import sympy  # noqa: F401
    sympy_available = True
except ImportError:
    sympy_available = False

_APPLICATION_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_AD_DIR = os.path.join(_APPLICATION_DIR, "automatic_differentiation")


@KratosUnittest.skipUnless(sympy_available, "sympy is not available")
@KratosUnittest.skipUnless(os.path.isdir(_AD_DIR), "the automatic_differentiation folder of the source tree is not available")
class TestSymbolicGeneration(KratosUnittest.TestCase):

    def test_alm_frictionless_2d2n(self):
        for path in (_AD_DIR, os.path.join(_AD_DIR, "ALM_frictionless_mortar_condition")):
            if path not in sys.path:
                sys.path.insert(0, path)
        import mortar_condition_generator as generator
        import compare_generated_conditions
        from generate_frictionless_mortar_condition import frictionless_functional

        spec = generator.ALM_FRICTIONLESS
        template_dir = os.path.join(_AD_DIR, "ALM_frictionless_mortar_condition")
        with tempfile.TemporaryDirectory() as output_dir:
            output_path = generator.Generate(spec, frictionless_functional, template_dir, output_dir, combinations=((2, 2, 2),), normal_variations=(False, True), log=lambda *args: None)
            with open(output_path) as generated_file:
                generated = generated_file.read()

        # Structure: one LHS body per normal variation, one RHS body (false) and one forwarder (true)
        self.assertEqual(generated.count("::CalculateLocalLHS("), 2)
        self.assertEqual(generated.count("::StaticCalculateLocalRHS("), 3)  # two specialisations + the forwarding call
        self.assertIn("<2,2, false, 2>::StaticCalculateLocalRHS(", generated)
        self.assertIn("<2,2, true, 2>::StaticCalculateLocalRHS(", generated)
        self.assertIn("DeltaNormalSlave[", generated)
        self.assertIn("DeltaDOperator[", generated)
        generated_blocks = generated.split("BEGIN AD REPLACEMENT")[1:]
        self.assertEqual(len(generated_blocks), 2)
        for block in generated_blocks:
            block = block.split("END AD REPLACEMENT")[0]
            for leftover in ("Derivative(", "//subsvar_", "Not supported", "TNumNodes", "MatrixSize"):
                self.assertNotIn(leftover, block)

        # Numerical equivalence with the committed file (2D2N blocks only)
        committed_path = os.path.join(_APPLICATION_DIR, "custom_conditions", spec.output_file)
        with open(committed_path) as committed_file:
            committed = committed_file.read()
        self.assertTrue(compare_generated_conditions.Compare(committed, generated, samples=2, tolerance=1.0e-9, require_same_blocks=False))


if __name__ == "__main__":
    KratosUnittest.main()
