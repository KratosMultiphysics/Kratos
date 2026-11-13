"""Tests for mesh-dataset curation: coherent random augmentations, seeding,
multi-series mixing, and zero-padded trajectory exports."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.export import mesh_export_process
from KratosMultiphysics.PhysicsNeMoApplication.training import torch_dataset
try:
    import torch
    import physicsnemo.mesh
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def _MakeMesh(seed=0):
    """Tetrahedral mesh carrying a vector field, a tensor field, a scalar and
    a non-spatial feature block (the shape upstream's bare `True` rejects)."""
    generator = torch.Generator().manual_seed(seed)
    points = torch.rand((8, 3), generator=generator, dtype=torch.float64)
    cells = torch.tensor([[0, 1, 2, 3], [4, 5, 6, 7]], dtype=torch.int64)
    return physicsnemo.mesh.Mesh(points=points, cells=cells, point_data={
        "VELOCITY": points.clone(),                                  # (N, 3) vector
        "STRESS": torch.eye(3, dtype=torch.float64).repeat(8, 1, 1),  # (N, 3, 3) tensor
        "PRESSURE": points[:, 0].clone(),                             # (N,) scalar
        "EMBEDDING": torch.arange(8 * 5, dtype=torch.float64).reshape(8, 5),
    })


class _StubReader:
    """Duck-typed physicsnemo mesh reader over in-memory meshes."""

    def __init__(self, meshes):
        self._meshes = list(meshes)

    def __len__(self):
        return len(self._meshes)

    def __getitem__(self, index):
        return self._meshes[index].clone(), {"index": index}


def _StubDataset(meshes, transforms=None, seed=-1):
    from physicsnemo.datapipes.mesh_dataset import MeshDataset
    dataset = MeshDataset(_StubReader(meshes), transforms=transforms)
    if seed >= 0:
        dataset.set_generator(torch.Generator().manual_seed(seed))
    return dataset


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python modules: torch, physicsnemo.")
class TestMeshAugmentations(KratosUnittest.TestCase):

    def test_RotationTransformsVectorsCoherently(self):
        transforms = torch_dataset.MakeMeshAugmentations(
            rotation=True, vector_fields=["VELOCITY"], tensor_fields=["STRESS"])
        self.assertEqual(len(transforms), 2)  # dtype cast + the rotation

        mesh = _MakeMesh()
        dataset = _StubDataset([mesh], transforms=transforms, seed=0)
        rotated, _ = dataset[0]

        # VELOCITY was a copy of the coordinates: it must still be, exactly
        self.assertTrue(torch.allclose(rotated.point_data["VELOCITY"], rotated.points))
        self.assertFalse(torch.allclose(rotated.points.double(), mesh.points, atol=1e-6))
        # rank-2 tensor: R I R^T == I, and the scalar/feature block untouched
        self.assertTrue(torch.allclose(
            rotated.point_data["STRESS"].double(), mesh.point_data["STRESS"], atol=1e-6))
        self.assertTrue(torch.allclose(
            rotated.point_data["PRESSURE"].double(), mesh.point_data["PRESSURE"], atol=1e-6))
        self.assertTrue(torch.allclose(
            rotated.point_data["EMBEDDING"].double(), mesh.point_data["EMBEDDING"], atol=1e-6))

    def test_ScaleIsCoherentAndTranslationLeavesVectorsAlone(self):
        mesh = _MakeMesh()
        scaled, _ = _StubDataset([mesh], transforms=torch_dataset.MakeMeshAugmentations(
            scale={"low": 2.0, "high": 2.0}, vector_fields=["VELOCITY"]), seed=1)[0]
        self.assertTrue(torch.allclose(scaled.points.double(), 2.0 * mesh.points, atol=1e-6))
        self.assertTrue(torch.allclose(
            scaled.point_data["VELOCITY"].double(), 2.0 * mesh.point_data["VELOCITY"], atol=1e-6))

        translated, _ = _StubDataset([mesh], transforms=torch_dataset.MakeMeshAugmentations(
            translation={"low": 1.0, "high": 1.0}, vector_fields=["VELOCITY"]), seed=1)[0]
        self.assertTrue(torch.allclose(translated.points.double(), mesh.points + 1.0, atol=1e-6))
        # translation must NOT touch vector values
        self.assertTrue(torch.allclose(
            translated.point_data["VELOCITY"].double(), mesh.point_data["VELOCITY"], atol=1e-6))

    def test_AxisAlignedRotationAndValidation(self):
        transforms = torch_dataset.MakeMeshAugmentations(
            rotation={"mode": "axis_aligned", "axes": ["z"], "low": 0.0, "high": 0.0},
            vector_fields=["VELOCITY"])
        mesh = _MakeMesh()
        rotated, _ = _StubDataset([mesh], transforms=transforms, seed=0)[0]
        self.assertTrue(torch.allclose(rotated.points.double(), mesh.points, atol=1e-6))

        with self.assertRaisesRegex(ValueError, "rotation mode"):
            torch_dataset.MakeMeshAugmentations(rotation={"mode": "quaternion"})
        with self.assertRaisesRegex(ValueError, "scale"):
            torch_dataset.MakeMeshAugmentations(scale="wide")

    def test_NoAugmentationsRequestedIsEmpty(self):
        self.assertEqual(torch_dataset.MakeMeshAugmentations(), [])

    def test_RandomnessRedrawnPerItemAndSeedReproducible(self):
        meshes = [_MakeMesh()]
        kwargs = dict(transforms=torch_dataset.MakeMeshAugmentations(
            rotation=True, vector_fields=["VELOCITY"]))

        dataset = _StubDataset(meshes, seed=7, **kwargs)
        first, _ = dataset[0]
        second, _ = dataset[0]
        self.assertFalse(torch.allclose(first.points, second.points))  # redrawn per item

        again = _StubDataset(meshes, seed=7, **kwargs)
        repeated, _ = again[0]
        self.assertTrue(torch.allclose(first.points, repeated.points))  # same seed, same draw

        epoched = _StubDataset(meshes, seed=7, **kwargs)
        epoched.set_epoch(3)
        self.assertFalse(torch.allclose(epoched[0][0].points, first.points))


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python modules: torch, physicsnemo.")
class TestMeshDatasetCuration(KratosUnittest.TestCase):

    def setUp(self):
        self.output_path = Path("test_mesh_curation")
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Curation")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        props = self.model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                 (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]):
            node = self.model_part.CreateNewNode(i + 1, *xyz)
            node.SetSolutionStepValue(Kratos.VELOCITY, [node.X, node.Y, node.Z])
        self.model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], props)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def _Export(self, path, steps, zero_pad_steps=0):
        process = mesh_export_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Curation",
                "list_of_fields"  : [ { "variable_name" : "VELOCITY",
                                        "data_location" : "node_historical" } ],
                "output_path"     : "%s",
                "zero_pad_steps"  : %d
            }
        }""" % (path, zero_pad_steps)), self.model)
        process.ExecuteInitialize()
        for step in steps:
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()

    def test_ZeroPaddedExportSortsInTimeOrder(self):
        directory = self.output_path / "padded"
        self._Export(directory, range(1, 12), zero_pad_steps=5)
        names = sorted(p.name for p in directory.glob("*.pmsh"))
        self.assertEqual(names[0], "mesh_00001.pmsh")
        self.assertEqual(names[-1], "mesh_00011.pmsh")

        unpadded = self.output_path / "plain"
        self._Export(unpadded, range(1, 12))
        plain_names = sorted(p.name for p in unpadded.glob("*.pmsh"))
        # the hazard the padding fixes: lexicographic order is not time order
        self.assertEqual(plain_names[1], "mesh_10.pmsh")

    def test_AugmentedDatasetOverExportedSeries(self):
        directory = self.output_path / "series"
        self._Export(directory, range(1, 4), zero_pad_steps=4)

        dataset = torch_dataset.CreateAugmentedMeshDataset(
            directory, rotation=True, scale={"low": 1.0, "high": 1.0},
            vector_fields=["VELOCITY"], seed=0)
        self.assertEqual(len(dataset), 3)
        mesh, metadata = dataset[0]
        self.assertIn("index", metadata)
        # exported VELOCITY equals the coordinates, so coherence survives export
        self.assertTrue(torch.allclose(
            mesh.point_data["VELOCITY"].double(), mesh.points.double(), atol=1e-6))

    def test_MultiMeshDatasetMixesSeries(self):
        first = self.output_path / "a"
        second = self.output_path / "b"
        self._Export(first, range(1, 4), zero_pad_steps=4)
        self._Export(second, range(1, 3), zero_pad_steps=4)

        mixed = torch_dataset.CreateMultiMeshDataset(
            [first, second], seed=0, rotation=True, vector_fields=["VELOCITY"])
        self.assertEqual(len(mixed), 5)
        _, metadata = mixed[4]
        self.assertEqual(metadata["dataset_index"], 1)
        _, first_metadata = mixed[0]
        self.assertEqual(first_metadata["dataset_index"], 0)

    def test_MultiMeshDatasetAcceptsReadyDatasetsAndValidates(self):
        directory = self.output_path / "ready"
        self._Export(directory, range(1, 3), zero_pad_steps=4)
        ready = torch_dataset.CreateMeshDataset(directory)
        mixed = torch_dataset.CreateMultiMeshDataset([ready], output_strict=False)
        self.assertEqual(len(mixed), 2)
        with self.assertRaisesRegex(ValueError, "at least one source"):
            torch_dataset.CreateMultiMeshDataset([])


if __name__ == '__main__':
    KratosUnittest.main()
