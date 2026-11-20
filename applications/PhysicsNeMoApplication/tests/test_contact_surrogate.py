"""Contact, the last compiled application nothing here exercised.

ContactStructuralMechanicsApplication is built in the reference environment
but no bridge had ever been pointed at it. This drives its smallest ALM
frictionless patch test and shows the contact pressure travelling through
the same tensor bridge and deployment path every other field uses.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_contact = kratos_utils.CheckIfApplicationsAvailable(
    "ContactStructuralMechanicsApplication", "StructuralMechanicsApplication",
    "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_contact, "Missing ContactStructuralMechanics applications.")
class TestContactCase(KratosUnittest.TestCase):
    def setUp(self):
        import contact_case

        if not contact_case.IsAvailable():
            self.skipTest("ContactStructuralMechanics test data is not present.")
        self.pressures, self.model_part = contact_case.SolveContact(Kratos.Model())

    def test_TheBlocksTouchAndCarryPressure(self):
        self.assertEqual(self.model_part.NumberOfNodes(), 8)   # two four-node blocks
        self.assertTrue(numpy.isfinite(self.pressures).all())
        in_contact = numpy.abs(self.pressures) > 1e-12
        self.assertGreater(int(in_contact.sum()), 0)
        # and the nodes away from the interface carry none
        self.assertLess(int(in_contact.sum()), self.model_part.NumberOfNodes())

    def test_TheImposedDisplacementIsTransmitted(self):
        import contact_case

        displacements = contact_case.Displacements(self.model_part)
        self.assertTrue(numpy.isfinite(displacements).all())
        self.assertGreater(numpy.abs(displacements).max(), 0.0)

    def test_TheFixtureIsReusedNotCopied(self):
        import contact_case

        self.assertIn("ContactStructuralMechanicsApplication",
                      str(contact_case._PARAMETERS_FILE))


@KratosUnittest.skipUnless(have_contact and have_torch,
                           "Missing ContactStructuralMechanics applications or torch.")
class TestContactFieldsThroughTheBridge(KratosUnittest.TestCase):
    """The contact pressure is an ordinary nodal field as far as the bridge
    is concerned - which is the point: no contact-specific plumbing."""

    def setUp(self):
        import contact_case

        if not contact_case.IsAvailable():
            self.skipTest("ContactStructuralMechanics test data is not present.")
        self.model = Kratos.Model()
        self.pressures, self.model_part = contact_case.SolveContact(self.model)
        self.checkpoint = Path("test_contact_surrogate.pt")

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def test_TheContactPressureGathersThroughTheTensorBridge(self):
        import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
        from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
            GetTensorAdaptor)

        adaptor = GetTensorAdaptor(
            self.model_part, "node_historical",
            CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)
        gathered = torch_bridge.KratosTensorToTorch(adaptor)
        numpy.testing.assert_allclose(
            gathered.numpy().reshape(-1), self.pressures, atol=1e-12)

    def test_ASurrogateDeploysOntoTheContactInterface(self):
        """A model writing a contact pressure needs nothing contact-aware:
        it is the same InferenceProcess contract as any other field."""
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            inference_process)

        class Doubler(torch.nn.Module):
            def forward(self, x):
                return 2.0 * x[..., :1]

        torch.jit.script(Doubler()).save(str(self.checkpoint))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Structure",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "torchscript",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "DISPLACEMENT_Y", "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "NODAL_PAUX",    "data_location" : "node_non_historical" } ]
            }
        }""" % self.checkpoint)
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        written = numpy.array([
            node.GetValue(Kratos.NODAL_PAUX) for node in self.model_part.Nodes])
        expected = 2.0 * numpy.array([
            node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y)
            for node in self.model_part.Nodes])
        numpy.testing.assert_allclose(written, expected, atol=1e-12)


if __name__ == '__main__':
    KratosUnittest.main()
