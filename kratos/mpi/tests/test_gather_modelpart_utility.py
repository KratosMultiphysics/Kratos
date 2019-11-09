from collections import namedtuple
from KratosMultiphysics.KratosUnittest import isclose
from math import exp
from itertools import chain


import KratosMultiphysics as Kratos
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KratosUnittest as UnitTest


class TestMesh:
    """An embedded test mdpa."""

    def __init__(self):
        MeshNode = namedtuple('MeshNode', ['Id', 'X', 'Y', 'Z'])
        self.nodes = [
            MeshNode(1, 0.00000, 0.00000, 0.00000),
            MeshNode(2, 1.00000, 1.00000, 0.00000),
            MeshNode(3, 1.00000, 1.00000, 0.00000),
            MeshNode(4, 0.00000, 0.00000, 0.00000),
            MeshNode(5, 2.00000, 1.00000, 0.00000),
            MeshNode(6, 2.00000, 0.00000, 0.00000),
            MeshNode(7, 3.00000, 1.00000, 0.00000),
            MeshNode(8, 3.00000, 0.00000, 0.00000),
        ]
        MeshElement = namedtuple(
            'MeshElement', ['Name', 'Id', 'NodeIds', 'PropertyId'])
        self.elements = [
            MeshElement('Element2D3N', 1, [1, 3, 4], 0),
            MeshElement('Element2D3N', 2, [1, 2, 3], 0),
            MeshElement('Element2D3N', 3, [2, 6, 3], 0),
            MeshElement('Element2D3N', 4, [2, 5, 6], 0),
            MeshElement('Element2D3N', 5, [5, 8, 6], 0),
            MeshElement('Element2D3N', 6, [5, 7, 8], 0),
        ]

    def LocalElements(self):
        """Return a sequence of elements for creation on the current process.

        Currently there is at most one element per process for testing.
        """
        mpi_comm = Kratos.DataCommunicator.GetDefault()
        if mpi_comm.Rank() < len(self.elements):
            yield self.elements[mpi_comm.Rank()]
        else:
            return

    def LocalNodes(self):
        """Return a set of nodes for creation on the current process."""
        connectivities = (el.NodeIds for el in self.LocalElements())
        nodes = (self.nodes[i-1] for i in chain(*connectivities))
        return set(nodes)

    def PartitionIndex(self, node_id):
        """Return the partition index for the given node id."""
        def predicate(el): return node_id in el.NodeIds
        try:
            return next(filter(predicate, self.elements)).Id - 1
        except StopIteration:
            return 0


def CreateLocalModelPart(model):
    model_part = model.CreateModelPart("LocalPart", 1)
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)
    KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
    model_part.SetBufferSize(1)
    prop = model_part.CreateNewProperties(0)
    mesh = TestMesh()
    for node in mesh.LocalNodes():
        model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
    for el in mesh.LocalElements():
        model_part.CreateNewElement(el.Name, el.Id, el.NodeIds, prop)
    for node in model_part.Nodes:
        partition_index = mesh.PartitionIndex(node.Id)
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, partition_index)
    return model_part


def SetPressureOnLocalNodes(model_part, unary_func):
    mpi_comm = Kratos.DataCommunicator.GetDefault()
    for node in model_part.Nodes:
        partition_index = node.GetSolutionStepValue(Kratos.PARTITION_INDEX)
        if partition_index == mpi_comm.Rank():
            node.SetSolutionStepValue(Kratos.PRESSURE, unary_func(node))


def SetPressureOnAllNodes(model_part, unary_func):
    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.PRESSURE, unary_func(node))


class TestGatherModelPartUtility(UnitTest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateLocalModelPart(self.model)
        self.gather_model_part = self.model.CreateModelPart("GatherPart", 1)
        self.master_rank = 0
        self.mesh_id = 0
        self.gather_model_part_util = KratosMPI.GatherModelPartUtility(
            self.master_rank, self.model_part, self.mesh_id, self.gather_model_part)

    def testGatherOnMaster(self):
        mpi_comm = Kratos.DataCommunicator.GetDefault()
        def test_func(node): return exp(-float(node.Id))
        SetPressureOnLocalNodes(
            self.model_part, unary_func=test_func)
        self.gather_model_part_util.GatherOnMaster(Kratos.PRESSURE)
        if mpi_comm.Rank() == self.master_rank:
            node_checks = (
                isclose(node.GetSolutionStepValue(Kratos.PRESSURE),
                        test_func(node)) for node in self.gather_model_part.Nodes
            )
            ok = min(node_checks, default=True)
            self.assertTrue(ok)
        mpi_comm.Barrier()

    def testScatterFromMaster(self):
        mpi_comm = Kratos.DataCommunicator.GetDefault()
        def test_func(node): return exp(float(node.Id))
        if mpi_comm.Rank() == self.master_rank:
            SetPressureOnAllNodes(
                self.gather_model_part, unary_func=test_func)
        self.gather_model_part_util.ScatterFromMaster(Kratos.PRESSURE)
        node_checks = (
            isclose(node.GetSolutionStepValue(Kratos.PRESSURE),
                    test_func(node)) for node in self.model_part.Nodes
        )
        ok = min(node_checks, default=True)
        self.assertTrue(ok)
        mpi_comm.Barrier()


def main():
    UnitTest.main()


if __name__ == '__main__':
    main()
