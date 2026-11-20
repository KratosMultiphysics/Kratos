"""MPI tests: mesh topology gathering, MPI-aware mesh export, CAE export.

The 3D fixture mirrors the 2D slab of test_mpi_dataset_export: every rank
builds its own x-slab of a structured unit-cube hexahedron mesh (plus a
"Skin" sub-part of quad conditions on the z = 0 face), with consistent
PARTITION_INDEX values and ParallelFillCommunicator. The hex+quad mix
exercises the variable-length connectivity encoding and both generic-name
tables. Rebuilding the fixture with world size 1 gives the serial reference
every rank-0 result must reproduce exactly.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes.export.cae_dataset_export_process import CaeDatasetExportProcess
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

try:
    import physicsnemo.mesh  # noqa: F401
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def _AnalyticValue(node):
    return 1.0 + 2.0 * node.X + 3.0 * node.Y + 4.0 * node.Z


def _BuildSlabModelPart(model, rank, size, divisions=3, name="Main"):
    """Builds this rank's x-slab of the structured hex cube (rank=0/size=1
    builds the whole serial reference)."""
    model_part = model.CreateModelPart(name)
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)
    properties = model_part.CreateNewProperties(1)

    n = divisions

    def cell_owner(cell_x):
        return min(cell_x * size // n, size - 1)

    def node_owner(node_x):
        return cell_owner(max(node_x - 1, 0))

    def node_id(ix, iy, iz):
        return ix * (n + 1) * (n + 1) + iy * (n + 1) + iz + 1

    my_columns = [cx for cx in range(n) if cell_owner(cx) == rank]
    needed = set()
    for cx in my_columns:
        for cy in range(n):
            for cz in range(n):
                for dx in (0, 1):
                    for dy in (0, 1):
                        for dz in (0, 1):
                            needed.add((cx + dx, cy + dy, cz + dz))
    for ix, iy, iz in sorted(needed):
        node = model_part.CreateNewNode(node_id(ix, iy, iz), ix / n, iy / n, iz / n)
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, node_owner(ix))

    skin_condition_ids = []
    for cx in my_columns:
        for cy in range(n):
            for cz in range(n):
                corners = [node_id(cx, cy, cz), node_id(cx + 1, cy, cz),
                           node_id(cx + 1, cy + 1, cz), node_id(cx, cy + 1, cz),
                           node_id(cx, cy, cz + 1), node_id(cx + 1, cy, cz + 1),
                           node_id(cx + 1, cy + 1, cz + 1), node_id(cx, cy + 1, cz + 1)]
                model_part.CreateNewElement(
                    "Element3D8N", cx * n * n + cy * n + cz + 1, corners, properties)
            # one quad condition per owned column on the z = 0 face
            condition_id = 10000 + cx * n + cy
            model_part.CreateNewCondition(
                "SurfaceCondition3D4N", condition_id,
                [node_id(cx, cy, 0), node_id(cx + 1, cy, 0),
                 node_id(cx + 1, cy + 1, 0), node_id(cx, cy + 1, 0)], properties)
            skin_condition_ids.append(condition_id)

    skin = model_part.CreateSubModelPart("Skin")
    skin.AddConditions(skin_condition_ids)
    skin.AddNodes([node_id(ix, iy, 0) for (ix, iy, iz) in needed if iz == 0])

    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.PRESSURE, _AnalyticValue(node))
    return model_part


def _CreateDistributedModelPart(model, divisions=3):
    import KratosMultiphysics.mpi as KratosMPI

    data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
    model_part = _BuildSlabModelPart(
        model, data_communicator.Rank(), data_communicator.Size(), divisions)
    KratosMPI.ParallelFillCommunicator(model_part, data_communicator).Execute()
    return model_part, data_communicator


def _CreateSerialReference(divisions=3):
    model = Kratos.Model()
    return model, _BuildSlabModelPart(model, 0, 1, divisions)


class TestMpiGatherModelPart(KratosUnittest.TestCase):
    def test_ShadowReproducesSerialTessellation(self):
        model = Kratos.Model()
        model_part, data_communicator = _CreateDistributedModelPart(model)

        gathered = distributed_utils.GatherModelPartToRank0(
            model_part, [("PRESSURE", "node_historical")], "Elements")

        if data_communicator.Rank() != 0:
            self.assertIsNone(gathered.model_part)
            self.assertIsNone(gathered.model)
            return
        reference_model, reference_part = _CreateSerialReference()
        shadow_provenance = domain_mesh_builder.BuildProvenance(gathered.model_part)
        reference_provenance = domain_mesh_builder.BuildProvenance(reference_part)
        numpy.testing.assert_allclose(
            shadow_provenance.simplex_points, reference_provenance.simplex_points)
        numpy.testing.assert_array_equal(
            shadow_provenance.simplex_cells, reference_provenance.simplex_cells)

        for node in gathered.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.PRESSURE), _AnalyticValue(node), places=12)

    def test_SkinConditionsGather(self):
        model = Kratos.Model()
        model_part, data_communicator = _CreateDistributedModelPart(model)
        skin = model_part.GetSubModelPart("Skin")

        gathered = distributed_utils.GatherModelPartToRank0(
            skin, [("PRESSURE", "node_historical")], "Conditions")

        if data_communicator.Rank() != 0:
            self.assertIsNone(gathered.model_part)
            return
        reference_model, reference_part = _CreateSerialReference()
        shadow_provenance = domain_mesh_builder.BuildProvenance(gathered.model_part, "Conditions")
        reference_provenance = domain_mesh_builder.BuildProvenance(
            reference_part.GetSubModelPart("Skin"), "Conditions")
        numpy.testing.assert_allclose(
            shadow_provenance.simplex_points, reference_provenance.simplex_points)
        numpy.testing.assert_array_equal(
            shadow_provenance.simplex_cells, reference_provenance.simplex_cells)
        self.assertAlmostEqual(
            float(shadow_provenance.ComputeSimplexMeasures().sum()), 1.0, places=12)

    def test_TwoDimensionalFixtureGathers(self):
        # the 2D triangle slab from the dataset-export tests exercises the
        # Element2D3N branch of the name table
        from test_mpi_dataset_export import _CreateDistributedModelPart as create_2d

        model = Kratos.Model()
        model_part, data_communicator = create_2d(model)
        gathered = distributed_utils.GatherModelPartToRank0(model_part, [], "Elements")
        if data_communicator.Rank() != 0:
            self.assertIsNone(gathered.model_part)
            return
        local_count = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        # collective SumAll happens on all ranks below
        self.assertEqual(gathered.model_part.NumberOfElements(), 2 * 8 * 8)

    def tearDown(self):
        Kratos.Testing.GetDefaultDataCommunicator().Barrier()


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestMpiMeshExport(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_mpi_mesh_export")

    def tearDown(self):
        data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
        data_communicator.Barrier()
        if data_communicator.Rank() == 0:
            KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))
        data_communicator.Barrier()

    def _RunExport(self, model, model_part_name, source_container):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export.mesh_export_process import MeshExportProcess

        process = MeshExportProcess(model, Kratos.Parameters("""{
            "model_part_name"  : "%s",
            "list_of_fields"   : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
            "source_container" : "%s",
            "output_path"      : "%s"
        }""" % (model_part_name, source_container, self.output_path)))
        process.ExecuteInitialize()
        process.ExecuteFinalizeSolutionStep()

    def test_DistributedElementExportMatchesSerial(self):
        model = Kratos.Model()
        model_part, data_communicator = _CreateDistributedModelPart(model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        self._RunExport(model, "Main", "Elements")
        data_communicator.Barrier()

        if data_communicator.Rank() != 0:
            return
        loaded = domain_mesh_builder.LoadMesh(self.output_path / "mesh_1.pmsh")

        reference_model, reference_part = _CreateSerialReference()
        reference_part.ProcessInfo[Kratos.STEP] = 1
        reference_mesh, _ = domain_mesh_builder.BuildMesh(
            reference_part, [(Kratos.PRESSURE, "node_historical")])
        numpy.testing.assert_allclose(
            numpy.asarray(loaded.points.cpu()), numpy.asarray(reference_mesh.points.cpu()))
        numpy.testing.assert_array_equal(
            numpy.asarray(loaded.cells.cpu()), numpy.asarray(reference_mesh.cells.cpu()))
        numpy.testing.assert_allclose(
            numpy.asarray(loaded.point_data["PRESSURE"].cpu()),
            numpy.asarray(reference_mesh.point_data["PRESSURE"].cpu()))

    def test_DistributedConditionExportWritesRankZeroOnly(self):
        model = Kratos.Model()
        model_part, data_communicator = _CreateDistributedModelPart(model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        self._RunExport(model, "Main.Skin", "Conditions")
        data_communicator.Barrier()

        if data_communicator.Rank() == 0:
            loaded = domain_mesh_builder.LoadMesh(self.output_path / "mesh_1.pmsh")
            self.assertEqual(numpy.asarray(loaded.cells.cpu()).shape, (18, 3))  # 9 quads x 2 tris


class TestMpiCaeDatasetExport(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_mpi_cae_export")

    def tearDown(self):
        data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
        data_communicator.Barrier()
        if data_communicator.Rank() == 0:
            KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))
        data_communicator.Barrier()

    def test_DistributedCaeExport(self):
        model = Kratos.Model()
        model_part, data_communicator = _CreateDistributedModelPart(model)
        model_part.ProcessInfo[Kratos.STEP] = 1

        process = CaeDatasetExportProcess(model, Kratos.Parameters("""{
            "model_part_name"         : "Main",
            "surface_model_part_name" : "Main.Skin",
            "surface_fields"          : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
            "volume_fields"           : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
            "output_path"             : "%s"
        }""" % self.output_path))
        process.ExecuteInitialize()
        process.ExecuteFinalizeSolutionStep()

        _, gathered_pressure = distributed_utils.GatherFieldToRank0(
            model_part, "PRESSURE", "node_historical")
        local_count = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        global_count = data_communicator.SumAll(local_count)
        data_communicator.Barrier()

        if data_communicator.Rank() != 0:
            return
        with numpy.load(self.output_path / "case_1.npz") as data:
            self.assertAlmostEqual(float(numpy.array(data["stl_areas"]).sum()), 1.0, places=6)
            self.assertEqual(numpy.array(data["volume_mesh_centers"]).shape, (global_count, 3))
            numpy.testing.assert_allclose(
                numpy.array(data["volume_fields"])[:, 0],
                numpy.asarray(gathered_pressure, dtype=numpy.float64).reshape(-1),
                rtol=1e-6)


if __name__ == '__main__':
    KratosUnittest.main()
