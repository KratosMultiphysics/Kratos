"""Tests for self-supervised geometry pretraining with AeroJEPA: a shape
family costs no solves, and what transfers is the trunk's encoder."""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.training import aerojepa_pretraining

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        from physicsnemo.experimental.models.aerojepa import AeroJEPA  # noqa: F401
    have_aerojepa = True
except ImportError:
    have_aerojepa = False

try:
    from physicsnemo.mesh import Mesh  # noqa: F401
    have_mesh = True
except ImportError:
    have_mesh = False


@KratosUnittest.skipUnless(have_torch and have_aerojepa and have_mesh,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestAeroJepaPretraining(KratosUnittest.TestCase):
    def setUp(self):
        self.model = aerojepa_pretraining.CreateAeroJepaModel(Kratos.Parameters("{}"))
        self.samples = aerojepa_pretraining.CreateGeometryFamilyDataset(
            count=6, seed=0, n_context=128, n_query=32, resolution=16)

    def test_TheFamilyCostsNoSolves(self):
        """Every label here comes from the geometry itself - that is what
        makes the pretraining self-supervised and free."""
        self.assertEqual(len(self.samples), 6)
        sample = self.samples[0]
        self.assertEqual(sample["context_pos"].shape[1], 3)
        self.assertEqual(sample["query_sdf"].shape, (32, 1))
        # a signed distance really is signed: points inside and outside
        self.assertLess(float(sample["query_sdf"].min()), 0.0)
        self.assertGreater(float(sample["query_sdf"].max()), 0.0)
        self.assertEqual({s["kind"] for s in self.samples}, {0, 1})

    def test_PretrainingImproves(self):
        history = aerojepa_pretraining.PretrainOnGeometry(
            self.model, self.samples, Kratos.Parameters("""{
                "epochs"  : 6,
                "device"  : "cpu",
                "seed"    : 0
            }"""))
        self.assertEqual(len(history), 6)
        self.assertTrue(numpy.isfinite(history).all())
        self.assertLess(min(history[1:]), history[0])

    def test_TheEmbeddingIsAGeometryDescriptor(self):
        """What EncodeGeometry delivers: a fixed-width vector that DEPENDS on
        the shape, which is what a downstream head consumes.

        It does not yet deliver a linearly separable shape space. Measured on
        this family after six epochs, the distance between the sphere and box
        means (0.74) is SMALLER than the spread within each class (2.00), so
        the classes overlap. That is the honest state of a few epochs on six
        shapes with no JEPA objective upstream, and the test records it
        rather than tuning until a separation claim passes.
        """
        aerojepa_pretraining.PretrainOnGeometry(
            self.model, self.samples,
            Kratos.Parameters('{"epochs": 6, "device": "cpu", "seed": 0}'))
        embeddings = [aerojepa_pretraining.EncodeGeometry(self.model, sample)
                      for sample in self.samples]
        self.assertEqual(embeddings[0].ndim, 1)
        self.assertTrue(all(e.shape == embeddings[0].shape for e in embeddings))
        self.assertTrue(all(torch.isfinite(e).all() for e in embeddings))
        # different geometries give different descriptors
        distinct = torch.stack(embeddings)
        self.assertGreater(float((distinct - distinct.mean(0)).norm(dim=1).min()), 0.0)

    def test_EmptySamplesRaise(self):
        with self.assertRaisesRegex(ValueError, "at least one sample"):
            aerojepa_pretraining.PretrainOnGeometry(
                self.model, [], Kratos.Parameters('{"epochs": 1}'))

    def test_TheTokenizerNeedsItsVoxelSize(self):
        """The trap: the default tokenizer strategy requires a voxel size and
        raises deep inside the forward, long after the model was built - so
        the bridge makes it a settings key with a default."""
        model = aerojepa_pretraining.CreateAeroJepaModel(Kratos.Parameters("{}"))
        sample = self.samples[0]
        prediction = model(
            context_pos=sample["context_pos"], context_feat=sample["context_feat"],
            gen_params=sample["gen_params"], query_pos=sample["query_pos"],
            query_sdf=torch.zeros_like(sample["query_sdf"]))
        self.assertEqual(tuple(prediction.shape), (32, 1))


if __name__ == '__main__':
    KratosUnittest.main()
