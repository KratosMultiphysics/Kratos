"""MPI tests: field gathering and MPI-aware dataset export.

Run with e.g.:
    mpiexec -np 2 python3 test_PhysicsNeMoApplication_mpi.py

The distributed model part is built without Metis: every rank creates its
own slab of a structured triangle mesh (owned plus interface nodes, with
consistent PARTITION_INDEX values) and ParallelFillCommunicator wires the
communication meshes.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes.export.dataset_export_process import DatasetExportProcess


def _AnalyticValue(node):
    return 1.0 + 2.0 * node.X + 3.0 * node.Y


def _CreateDistributedModelPart(model, divisions=8):
    """Each rank builds its own x-slab of a structured triangle mesh (owned
    plus interface nodes, consistent PARTITION_INDEX everywhere) and
    ParallelFillCommunicator computes the communication meshes."""
    import KratosMultiphysics.mpi as KratosMPI

    data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
    rank, size = data_communicator.Rank(), data_communicator.Size()

    model_part = model.CreateModelPart("Main")
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)
    properties = model_part.CreateNewProperties(1)

    n = divisions

    def cell_owner(cell_x):
        return min(cell_x * size // n, size - 1)

    def node_owner(node_x):
        # a node belongs to the owner of the cell column at its left
        return cell_owner(max(node_x - 1, 0))

    def node_id(node_x, node_y):
        return node_x * (n + 1) + node_y + 1

    my_columns = [cell_x for cell_x in range(n) if cell_owner(cell_x) == rank]
    needed = set()
    for cell_x in my_columns:
        for cell_y in range(n):
            for dx, dy in ((0, 0), (1, 0), (1, 1), (0, 1)):
                needed.add((cell_x + dx, cell_y + dy))
    for node_x, node_y in sorted(needed):
        node = model_part.CreateNewNode(node_id(node_x, node_y), node_x / n, node_y / n, 0.0)
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, node_owner(node_x))
    for cell_x in my_columns:
        for cell_y in range(n):
            quad = [node_id(cell_x, cell_y), node_id(cell_x + 1, cell_y),
                    node_id(cell_x + 1, cell_y + 1), node_id(cell_x, cell_y + 1)]
            base = 2 * (cell_x * n + cell_y)  # globally unique element ids
            model_part.CreateNewElement("Element2D3N", base + 1, quad[:3], properties)
            model_part.CreateNewElement("Element2D3N", base + 2, [quad[0], quad[2], quad[3]], properties)

    KratosMPI.ParallelFillCommunicator(model_part, data_communicator).Execute()

    for node in model_part.Nodes:  # local AND ghost nodes carry the field
        node.SetSolutionStepValue(Kratos.PRESSURE, _AnalyticValue(node))
    return model_part, data_communicator


class TestMpiGatherHelpers(KratosUnittest.TestCase):
    def test_NodalFieldGatheredSortedAndComplete(self):
        model = Kratos.Model()
        model_part, data_communicator = _CreateDistributedModelPart(model)

        ids, values = distributed_utils.GatherFieldToRank0(
            model_part, "PRESSURE", "node_historical")

        local_count = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        global_count = data_communicator.SumAll(local_count)
        if data_communicator.Rank() == 0:
            self.assertEqual(len(ids), global_count)
            self.assertTrue((numpy.diff(ids) > 0).all())  # ids unique + ascending
            # every gathered value is the analytic field of its node id
            id_to_value = dict(zip(ids.tolist(), numpy.asarray(values).reshape(len(ids), -1)[:, 0]))
            for node in model_part.GetCommunicator().LocalMesh().Nodes:
                self.assertAlmostEqual(id_to_value[node.Id], _AnalyticValue(node), places=12)
        else:
            self.assertIsNone(ids)
            self.assertIsNone(values)


class TestMpiDatasetExport(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_mpi_dataset_export")

    def tearDown(self):
        data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
        data_communicator.Barrier()
        if data_communicator.Rank() == 0:
            KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))
        data_communicator.Barrier()

    def test_DistributedExportMatchesSerialLayout(self):
        model = Kratos.Model()
        model_part, data_communicator = _CreateDistributedModelPart(model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        model_part.ProcessInfo[Kratos.TIME] = 1.0

        process = DatasetExportProcess(model, Kratos.Parameters("""{
            "model_part_name" : "Main",
            "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
            "output_path"     : "test_mpi_dataset_export"
        }"""))
        process.ExecuteInitialize()
        process.ExecuteFinalizeSolutionStep()
        data_communicator.Barrier()

        sample_file = self.output_path / "sample_1.npz"
        if data_communicator.Rank() == 0:
            self.assertTrue(sample_file.is_file())
            with numpy.load(sample_file) as data:
                exported = numpy.asarray(data["PRESSURE__node_historical"])
                self.assertEqual(int(data["STEP"]), 1)
        else:
            exported = None

        local_count = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        global_count = data_communicator.SumAll(local_count)
        # collective on ALL ranks (never inside a rank guard)
        ids, gathered = distributed_utils.GatherFieldToRank0(
            model_part, "PRESSURE", "node_historical")
        if data_communicator.Rank() == 0:
            # one row per global node, id-sorted: the serial layout invariant
            self.assertEqual(exported.reshape(global_count, -1).shape[0], global_count)
            reference = numpy.sort(exported.reshape(global_count, -1)[:, 0])
            # values match the analytic field evaluated on the id-sorted lattice
            numpy.testing.assert_allclose(
                exported.reshape(global_count, -1),
                numpy.asarray(gathered).reshape(global_count, -1))
            self.assertGreater(reference.max(), reference.min())


if __name__ == '__main__':
    KratosUnittest.main()
