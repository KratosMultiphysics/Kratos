"""Distributed co-simulation surrogates: the wrapper on many ranks.

Run with:
    mpiexec -np 2 python3 test_PhysicsNeMoApplication_mpi.py

Covers the two distributed mesh-import paths ("partition_mdpa" true/false),
proves per-rank inference reproduces the serial prediction node for node,
and runs a full CoSimulationAnalysis in the distributed<->distributed and
subgroup (N-rank <-> M-rank) configurations.
"""

from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

have_cosim = KratosUtilities.CheckIfApplicationsAvailable("CoSimulationApplication")
have_mapping = KratosUtilities.CheckIfApplicationsAvailable("MappingApplication")
have_metis = KratosUtilities.CheckIfApplicationsAvailable("MetisApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

_FIXTURE_DIR = Path(__file__).parent / "cosim_cases"
# The meshed fixture: Metis partitions a connectivity graph, so the 8-node
# element-free surrogate_interface.mdpa cannot be used here.
_MDPA_FILE = str(_FIXTURE_DIR / "surrogate_interface_meshed")


def _SaveAffineModel(path, scale, offset):
    """TorchScript toy: y = scale * x + offset, elementwise on (N, 3)."""
    class Affine(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.scale = scale
            self.offset = offset

        def forward(self, x):
            return self.scale * x + self.offset

    torch.jit.script(Affine()).save(str(path))


def _SolverBlock(checkpoint, time_step=0.0, distributed=True, partition_mdpa=True,
                 num_processes=0, data_communicator_name=""):
    settings = Kratos.Parameters("""{
        "type" : "KratosMultiphysics.PhysicsNeMoApplication.deployment.cosim_surrogate_solver_wrapper",
        "solver_wrapper_settings" : {
            "mdpa_file"       : "%s",
            "time_step"       : %f,
            "distributed"     : %s,
            "partition_mdpa"  : %s,
            "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
            "input_fields"    : [ { "variable_name" : "FORCE",        "data_location" : "node_historical" } ],
            "output_fields"   : [ { "variable_name" : "DISPLACEMENT", "data_location" : "node_historical" } ]
        },
        "data" : {
            "load" : { "model_part_name" : "Surrogate", "variable_name" : "FORCE",        "dimension" : 3 },
            "disp" : { "model_part_name" : "Surrogate", "variable_name" : "DISPLACEMENT", "dimension" : 3 }
        }
    }""" % (_MDPA_FILE, time_step, "true" if distributed else "false",
            "true" if partition_mdpa else "false", checkpoint))
    if num_processes > 0:
        mpi_settings = Kratos.Parameters("""{
            "num_processes" : 0, "data_communicator_name" : ""
        }""")
        mpi_settings["num_processes"].SetInt(num_processes)
        mpi_settings["data_communicator_name"].SetString(data_communicator_name)
        settings.AddValue("mpi_settings", mpi_settings)
    return settings


def _CreateWrapper(settings, model, name):
    from KratosMultiphysics.PhysicsNeMoApplication.deployment.cosim_surrogate_solver_wrapper import (
        CoSimSurrogateSolverWrapper)
    return CoSimSurrogateSolverWrapper(settings, model, name)


@KratosUnittest.skipUnless(have_cosim, "Missing required application: CoSimulationApplication.")
@KratosUnittest.skipUnless(have_metis, "Missing required application: MetisApplication.")
@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestMpiDistributedSurrogateWrapper(KratosUnittest.TestCase):
    """The wrapper itself on a distributed model part."""

    def setUp(self):
        self.data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
        self.checkpoint = Path(f"test_mpi_cosim_affine_{self.data_communicator.Rank()}.pt")
        _SaveAffineModel(self.checkpoint, 2.0, 1.0)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))
        self.data_communicator.Barrier()

    def _Solve(self, partition_mdpa, seed_with_ids=False):
        model = Kratos.Model()
        wrapper = _CreateWrapper(
            _SolverBlock(self.checkpoint, partition_mdpa=partition_mdpa),
            model, "surrogate")
        for node in wrapper.model_part.Nodes:
            # seed_with_ids encodes the node Id in the field, so a gathered
            # row can be checked against its own id without gathering the
            # coordinates through the same helper - which would permute
            # identically and hide exactly the bug being looked for
            force = ([float(node.Id), 0.0, 0.0] if seed_with_ids
                     else [node.X, 2.0 * node.Y, 3.0 * node.Z])
            node.SetSolutionStepValue(Kratos.FORCE, force)
        wrapper.Initialize()
        wrapper.SolveSolutionStep()
        return wrapper

    def test_OwnedNodesPartitionTheGlobalMesh(self):
        for partition_mdpa in (True, False):
            with self.subTest(partition_mdpa=partition_mdpa):
                wrapper = self._Solve(partition_mdpa)
                communicator = wrapper.model_part.GetCommunicator()
                self.assertTrue(wrapper.model_part.IsDistributed())
                owned = communicator.LocalMesh().NumberOfNodes()
                # exact disjoint cover: no node owned twice, none unowned
                self.assertEqual(self.data_communicator.SumAll(owned),
                                 communicator.GlobalNumberOfNodes())

    def test_DistributedInferenceMatchesSerial(self):
        # The real claim: partitioning changes nothing about the answer.
        for partition_mdpa in (True, False):
            with self.subTest(partition_mdpa=partition_mdpa):
                wrapper = self._Solve(partition_mdpa)
                communicator = wrapper.model_part.GetCommunicator()
                for node in communicator.LocalMesh().Nodes:
                    displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
                    expected = (2.0 * node.X + 1.0,
                                2.0 * 2.0 * node.Y + 1.0,
                                2.0 * 3.0 * node.Z + 1.0)
                    for value, reference in zip(displacement, expected):
                        self.assertAlmostEqual(value, reference, places=6)

    def test_GhostsAreSynchronizedAfterSolve(self):
        # Owned rows are authoritative; SynchronizeVariable must leave every
        # ghost carrying its owner's value, or the coupling reads garbage.
        wrapper = self._Solve(partition_mdpa=True)
        ghost_mesh = wrapper.model_part.GetCommunicator().GhostMesh()
        for node in ghost_mesh.Nodes:
            displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            expected = (2.0 * node.X + 1.0,
                        2.0 * 2.0 * node.Y + 1.0,
                        2.0 * 3.0 * node.Z + 1.0)
            for value, reference in zip(displacement, expected):
                self.assertAlmostEqual(value, reference, places=6)

    def test_GatheredPredictionReproducesSerialLayout(self):
        from KratosMultiphysics.PhysicsNeMoApplication.distributed import graph_partition_utils
        import numpy

        wrapper = self._Solve(partition_mdpa=True, seed_with_ids=True)
        model_part = wrapper.model_part
        communicator = model_part.GetCommunicator()
        owned_ids = {node.Id for node in communicator.LocalMesh().Nodes}
        node_ids = numpy.array([node.Id for node in model_part.Nodes], dtype=numpy.int64)
        owned_mask = numpy.array([node.Id in owned_ids for node in model_part.Nodes])
        values = numpy.array(
            [list(node.GetSolutionStepValue(Kratos.DISPLACEMENT)) for node in model_part.Nodes])

        ids, gathered = graph_partition_utils.GatherOwnedPredictionsToRank0(
            model_part, node_ids, owned_mask, values, self.data_communicator)

        # GlobalNumberOfNodes is collective - every rank must reach it, so it
        # cannot live inside the rank-zero branch below.
        global_nodes = communicator.GlobalNumberOfNodes()
        if self.data_communicator.Rank() == 0:
            self.assertEqual(len(ids), global_nodes)
            self.assertEqual(sorted(ids.tolist()), ids.tolist())  # serial id order
            # ...and every VALUE must belong to its own id. The model is
            # 2x + 1 and the input carried the node Id, so row k must be
            # 2*ids[k] + 1. Asserting only that the ids are sorted lets any
            # permutation of the values through.
            numpy.testing.assert_allclose(
                gathered[:, 0], 2.0 * ids.astype(float) + 1.0, rtol=1e-9)
        else:
            self.assertIsNone(ids)

    def test_CouplingInterfaceDataUsesTheOwnedLayout(self):
        wrapper = self._Solve(partition_mdpa=True)
        wrapper.Initialize()
        local_size = wrapper.data_dict["disp"].Size()
        global_nodes = wrapper.model_part.GetCommunicator().GlobalNumberOfNodes()
        # 3 components per node, ghost-free, summing to the global mesh
        self.assertEqual(self.data_communicator.SumAll(local_size), 3 * global_nodes)

    def test_NonPartitionSafeInterfaceIsRejected(self):
        settings = _SolverBlock(self.checkpoint)
        settings["solver_wrapper_settings"]["model_interface"].SetString("transolver") \
            if settings["solver_wrapper_settings"].Has("model_interface") \
            else settings["solver_wrapper_settings"].AddString("model_interface", "transolver")
        with self.assertRaisesRegex(ValueError, "not partition-safe"):
            _CreateWrapper(settings, Kratos.Model(), "surrogate")

    def test_SerialDefaultIsUnchangedUnderMpi(self):
        # "distributed" defaults to false: the wrapper must still take the
        # rank-zero communicator and build a non-distributed model part.
        settings = _SolverBlock(self.checkpoint, distributed=False)
        wrapper = _CreateWrapper(settings, Kratos.Model(), "surrogate")
        self.assertEqual(wrapper.data_communicator.IsDefinedOnThisRank(),
                         self.data_communicator.Rank() == 0)
        if self.data_communicator.Rank() == 0:
            self.assertFalse(wrapper.model_part.IsDistributed())
        else:
            self.assertIsNone(wrapper.model_part)


@KratosUnittest.skipUnless(have_cosim, "Missing required application: CoSimulationApplication.")
@KratosUnittest.skipUnless(have_mapping, "Missing required application: MappingApplication.")
@KratosUnittest.skipUnless(have_metis, "Missing required application: MetisApplication.")
@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestMpiDistributedSurrogateCoupledLoop(KratosUnittest.TestCase):
    """A full CoSimulationAnalysis with distributed surrogate wrappers."""

    def setUp(self):
        self.data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
        rank = self.data_communicator.Rank()
        self.affine = Path(f"test_mpi_cosim_loop_affine_{rank}.pt")
        self.identity = Path(f"test_mpi_cosim_loop_identity_{rank}.pt")
        _SaveAffineModel(self.affine, 2.0, 1.0)
        _SaveAffineModel(self.identity, 1.0, 0.0)

    def tearDown(self):
        for path in (self.affine, self.identity):
            KratosUtilities.DeleteFileIfExisting(str(path))
        self.data_communicator.Barrier()

    @staticmethod
    def _AnalysisParameters():
        return Kratos.Parameters("""{
            "problem_data" : {
                "start_time"    : 0.0,
                "end_time"      : 2.0,
                "echo_level"    : 0,
                "print_colors"  : false,
                "parallel_type" : "MPI"
            },
            "solver_settings" : {
                "type"       : "coupled_solvers.gauss_seidel_weak",
                "echo_level" : 0,
                "data_transfer_operators" : {
                    "mapper" : {
                        "type"            : "kratos_mapping",
                        "mapper_settings" : { "mapper_type" : "nearest_neighbor" }
                    }
                },
                "coupling_sequence" : [
                    {
                        "name"             : "predictor",
                        "input_data_list"  : [],
                        "output_data_list" : []
                    },
                    {
                        "name"            : "receiver",
                        "input_data_list" : [{
                            "data"                   : "load",
                            "from_solver"            : "predictor",
                            "from_solver_data"       : "disp",
                            "data_transfer_operator" : "mapper"
                        }],
                        "output_data_list" : []
                    }
                ],
                "solvers" : {}
            }
        }""")

    def _RunLoop(self, predictor_block, receiver_block):
        from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import (
            CoSimulationAnalysis)

        parameters = self._AnalysisParameters()
        parameters["solver_settings"]["solvers"].AddValue("predictor", predictor_block)
        parameters["solver_settings"]["solvers"].AddValue("receiver", receiver_block)

        analysis = CoSimulationAnalysis(parameters)
        analysis.Initialize()

        wrappers = analysis._GetSolver().solver_wrappers
        predictor = wrappers["predictor"]
        if predictor.IsDefinedOnThisRank():
            for node in predictor.model_part.Nodes:
                node.SetSolutionStepValue(Kratos.FORCE, [node.X, 2.0 * node.Y, 3.0 * node.Z])

        analysis.RunSolutionLoop()
        analysis.Finalize()
        return wrappers

    def _AssertReceiverGotThePredictorsAnswer(self, receiver):
        if not receiver.IsDefinedOnThisRank():
            return
        nodes = receiver.model_part.GetCommunicator().LocalMesh().Nodes \
            if receiver.model_part.IsDistributed() else receiver.model_part.Nodes
        checked = 0
        for node in nodes:
            displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            expected = (2.0 * node.X + 1.0,
                        2.0 * 2.0 * node.Y + 1.0,
                        2.0 * 3.0 * node.Z + 1.0)
            for value, reference in zip(displacement, expected):
                self.assertAlmostEqual(value, reference, places=6)
            checked += 1
        self.assertGreater(checked, 0)

    def test_DistributedToDistributed(self):
        wrappers = self._RunLoop(
            _SolverBlock(self.affine, time_step=1.0),
            _SolverBlock(self.identity))
        self._AssertReceiverGotThePredictorsAnswer(wrappers["receiver"])

    def test_DistributedToRankZero(self):
        # Mixed group sizes: the receiver stays on the serial (rank-zero)
        # path while the predictor spans the world. Note the reverse of this
        # - BOTH sides rank-zero - deadlocks in CoSimulation's kratos_mapping
        # operator, which is why the wrapper documents keeping one side
        # distributed.
        wrappers = self._RunLoop(
            _SolverBlock(self.affine, time_step=1.0),
            _SolverBlock(self.identity, distributed=False))
        self._AssertReceiverGotThePredictorsAnswer(wrappers["receiver"])

    def test_SubgroupOfRanks(self):
        # The receiver lives on the first N-1 ranks only; the remaining rank
        # holds an UndefinedSolver and still enters every collective.
        size = self.data_communicator.Size()
        if size < 3:
            self.skipTest("needs at least 3 ranks to leave a rank out of the subgroup")
        wrappers = self._RunLoop(
            _SolverBlock(self.affine, time_step=1.0),
            _SolverBlock(self.identity, num_processes=size - 1,
                         data_communicator_name="surrogate_receiver_group"))
        receiver = wrappers["receiver"]
        self.assertEqual(receiver.IsDefinedOnThisRank(),
                         self.data_communicator.Rank() < size - 1)
        self._AssertReceiverGotThePredictorsAnswer(receiver)


if __name__ == '__main__':
    KratosUnittest.main()
