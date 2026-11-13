"""Tests for the scalable GNN variants: the bistride multiscale hierarchy,
proximity "world" edges, and all three model interfaces through the process."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import graph_inference_process
from test_grid_bridge import CreateStructuredTetModelPart

sys.path.insert(0, str(Path(__file__).parent / "kratos_solver_cases"))

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.nn.module.gnn_layers.graph_types import PYG_AVAILABLE
    from physicsnemo.models.meshgraphnet import (
        BiStrideMeshGraphNet, HybridMeshGraphNet, MeshGraphKAN)
    have_scalable_gnns = PYG_AVAILABLE
except ImportError:
    have_scalable_gnns = False

have_cuda = have_torch and torch.cuda.is_available()

_MISSING = "Missing required python modules: physicsnemo with torch_geometric/torch_scatter."


def _PathGraph(n):
    """(edge_index, positions) of a bidirectional 1D chain of n nodes."""
    pairs = numpy.array([[i, i + 1] for i in range(n - 1)], dtype=numpy.int64).T
    edge_index = numpy.concatenate([pairs, pairs[::-1]], axis=1)
    positions = numpy.stack(
        [numpy.arange(n, dtype=float), numpy.zeros(n), numpy.zeros(n)], axis=1)
    return edge_index, positions


class TestBistrideHierarchy(KratosUnittest.TestCase):
    """Pure numpy/scipy: no torch or physicsnemo needed."""

    def test_PathGraphMatchesReferenceHierarchy(self):
        edge_index, positions = _PathGraph(16)
        ms_edges, ms_ids = graph_bridge.BuildBistrideHierarchy(
            edge_index, 16, positions, num_levels=2)

        # the reference BSMS output for this graph: parity selection halving
        # each level, two-hop edges between the survivors
        self.assertEqual([edges.shape for edges in ms_edges],
                         [(2, 30), (2, 14), (2, 6)])
        self.assertEqual([ids.shape for ids in ms_ids], [(8,), (4,)])
        self.assertEqual(ms_ids[0].tolist(), [1, 3, 5, 7, 9, 11, 13, 15])
        self.assertEqual(ms_ids[1].tolist(), [1, 3, 5, 7])

    def test_LengthContractAndLocalRenumbering(self):
        edge_index, positions = _PathGraph(16)
        levels = 2
        ms_edges, ms_ids = graph_bridge.BuildBistrideHierarchy(
            edge_index, 16, positions, num_levels=levels)

        # the model's contract: the two lists differ in length by one
        self.assertEqual(len(ms_edges), levels + 1)
        self.assertEqual(len(ms_ids), levels)

        node_counts = [16] + [ids.size for ids in ms_ids]
        for level, (edges, count) in enumerate(zip(ms_edges, node_counts)):
            self.assertTrue((edges >= 0).all(), msg=f"level {level}")
            self.assertLess(int(edges.max()), count,
                            msg=f"level {level} edges must index its own node space")
        for level, ids in enumerate(ms_ids):
            self.assertLess(int(ids.max()), node_counts[level],
                            msg=f"ms_ids[{level}] indexes the PARENT level")
            self.assertEqual(len(set(ids.tolist())), ids.size)
        # levels shrink monotonically
        self.assertEqual(node_counts, sorted(node_counts, reverse=True))

    def test_DisconnectedComponentsEachGetRepresentation(self):
        # two disjoint chains: every component must survive coarsening
        first, positions_first = _PathGraph(8)
        second, positions_second = _PathGraph(8)
        edge_index = numpy.concatenate([first, second + 8], axis=1)
        positions = numpy.concatenate(
            [positions_first, positions_second + numpy.array([100.0, 0.0, 0.0])])

        ms_edges, ms_ids = graph_bridge.BuildBistrideHierarchy(
            edge_index, 16, positions, num_levels=1)
        kept = ms_ids[0]
        self.assertTrue((kept < 8).any(), msg="first component vanished")
        self.assertTrue((kept >= 8).any(), msg="second component vanished")

    def test_IsolatedNodeSurvivesAndKeepsASelfConsistentLevel(self):
        # an isolated node is its own component; upstream's message passing
        # divides by node degree, so it must not be dropped or left degenerate
        edge_index, positions = _PathGraph(8)
        positions = numpy.concatenate([positions, numpy.array([[50.0, 0.0, 0.0]])])
        ms_edges, ms_ids = graph_bridge.BuildBistrideHierarchy(
            edge_index, 9, positions, num_levels=1)
        self.assertIn(8, ms_ids[0].tolist())

    def test_Validation(self):
        edge_index, positions = _PathGraph(16)
        with self.assertRaisesRegex(ValueError, "num_levels"):
            graph_bridge.BuildBistrideHierarchy(edge_index, 16, positions, num_levels=0)
        with self.assertRaisesRegex(ValueError, "collapsed"):
            graph_bridge.BuildBistrideHierarchy(edge_index, 16, positions, num_levels=6)


class TestWorldEdges(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "World", 3, historical_variables=(Kratos.PRESSURE,))

    def test_RadiusEdgesMatchBruteForce(self):
        radius = 0.4
        edge_index, edge_features = graph_bridge.BuildWorldEdges(
            self.model_part, Kratos.Parameters(
                '{"type": "radius", "radius": %f, "backend": "numpy"}' % radius))

        positions = graph_bridge.NodePositions(self.model_part)
        distances = numpy.linalg.norm(
            positions[:, None, :] - positions[None, :, :], axis=2)
        numpy.fill_diagonal(distances, numpy.inf)
        expected = set(map(tuple, numpy.argwhere(distances < radius).tolist()))
        produced = set(map(tuple, edge_index.T.tolist()))
        self.assertEqual(produced, expected)

        # symmetric, and features follow the mesh-graph convention
        self.assertEqual(produced, {(j, i) for i, j in produced})
        self.assertEqual(edge_features.shape, (edge_index.shape[1], 4))
        relative = positions[edge_index[1]] - positions[edge_index[0]]
        numpy.testing.assert_allclose(edge_features[:, :3], relative, atol=1e-12)
        numpy.testing.assert_allclose(
            edge_features[:, 3], numpy.linalg.norm(relative, axis=1), atol=1e-12)

    def test_ConcatenationPutsMeshEdgesFirst(self):
        _, mesh_edge_index, _, _ = graph_bridge.BuildGraph(self.model_part)
        world_edge_index, _ = graph_bridge.BuildWorldEdges(
            self.model_part, Kratos.Parameters(
                '{"type": "radius", "radius": 0.4, "backend": "numpy"}'))
        combined = graph_bridge.ConcatenateEdgeSets(mesh_edge_index, world_edge_index)

        # the model splits the concatenated set positionally, so the mesh
        # block must occupy exactly the leading columns
        self.assertEqual(combined.shape[1],
                         mesh_edge_index.shape[1] + world_edge_index.shape[1])
        numpy.testing.assert_array_equal(
            combined[:, :mesh_edge_index.shape[1]], mesh_edge_index)
        numpy.testing.assert_array_equal(
            combined[:, mesh_edge_index.shape[1]:], world_edge_index)


@KratosUnittest.skipUnless(have_scalable_gnns, _MISSING)
class TestScalableGnnsThroughProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.checkpoint = "test_scalable_gnn.mdlus"
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Scalable", 3,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(self.checkpoint)

    def _Save(self, model):
        model.save(self.checkpoint)

    def _RunProcess(self, model_interface, extra_settings=""):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Scalable",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "model_interface" : "%s",
                "input_fields"    : [ { "variable_name" : "PRESSURE",
                                        "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE",
                                        "data_location" : "node_historical" } ]
                %s
            }
        }""" % (self.checkpoint, model_interface, extra_settings))
        process = graph_inference_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        return numpy.array([node.GetSolutionStepValue(Kratos.TEMPERATURE)
                            for node in self.model_part.Nodes])

    def _AssertWritten(self, values):
        self.assertEqual(values.shape, (self.model_part.NumberOfNodes(),))
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(float(values.std()), 0.0)

    def test_BistrideThroughProcess(self):
        torch.manual_seed(0)
        self._Save(BiStrideMeshGraphNet(
            input_dim_nodes=1, input_dim_edges=4, output_dim=1, processor_size=2,
            hidden_dim_processor=16, hidden_dim_node_encoder=8, hidden_dim_edge_encoder=8,
            hidden_dim_node_decoder=8, num_layers_bistride=1,
            num_mesh_levels=1, bistride_pos_dim=3))
        self._AssertWritten(self._RunProcess("bistride", ', "multiscale_levels": 1'))

    def test_HybridThroughProcess(self):
        torch.manual_seed(1)
        self._Save(HybridMeshGraphNet(
            input_dim_nodes=1, input_dim_edges=4, output_dim=1, processor_size=2,
            hidden_dim_processor=16, hidden_dim_node_encoder=8,
            hidden_dim_edge_encoder=8, hidden_dim_node_decoder=8))
        self._AssertWritten(self._RunProcess(
            "hybrid", ', "world_edges": { "type": "radius", "radius": 0.4, "backend": "numpy" }'))

    def test_MeshGraphKanThroughProcess(self):
        torch.manual_seed(2)
        self._Save(MeshGraphKAN(
            input_dim_nodes=1, input_dim_edges=4, output_dim=1, processor_size=2,
            hidden_dim_processor=16, hidden_dim_edge_encoder=8,
            hidden_dim_node_decoder=8, num_harmonics=3))
        self._AssertWritten(self._RunProcess("meshgraphkan"))

    def test_UnknownInterfaceRejected(self):
        with self.assertRaisesRegex(ValueError, "Unsupported model interface"):
            self._RunProcess("bistride_v2")

    def test_MultiscaleLevelsValidated(self):
        with self.assertRaisesRegex(ValueError, "multiscale_levels"):
            self._RunProcess("bistride", ', "multiscale_levels": 0')

    def test_BistrideWithoutMultiscaleTablesRaises(self):
        torch.manual_seed(0)
        model = BiStrideMeshGraphNet(
            input_dim_nodes=1, input_dim_edges=4, output_dim=1, processor_size=2,
            hidden_dim_processor=16, hidden_dim_node_encoder=8, hidden_dim_edge_encoder=8,
            hidden_dim_node_decoder=8, num_layers_bistride=1,
            num_mesh_levels=1, bistride_pos_dim=3)
        node_features, edge_index, edge_features, node_ids = graph_bridge.BuildGraph(
            self.model_part, (("PRESSURE", "node_historical"),))
        graph = graph_bridge.ToPyGGraph(edge_index, len(node_ids),
                                        graph_bridge.NodePositions(self.model_part))
        with self.assertRaisesRegex(ValueError, "multiscale"):
            graph_inference_process.RunGraphForward(
                model, torch.device("cpu"), "bistride", node_features, edge_features, graph)

    def test_HybridWithoutWorldEdgesRaises(self):
        torch.manual_seed(1)
        model = HybridMeshGraphNet(
            input_dim_nodes=1, input_dim_edges=4, output_dim=1, processor_size=2,
            hidden_dim_processor=16, hidden_dim_node_encoder=8,
            hidden_dim_edge_encoder=8, hidden_dim_node_decoder=8)
        node_features, edge_index, edge_features, node_ids = graph_bridge.BuildGraph(
            self.model_part, (("PRESSURE", "node_historical"),))
        graph = graph_bridge.ToPyGGraph(edge_index, len(node_ids))
        with self.assertRaisesRegex(ValueError, "world edge"):
            graph_inference_process.RunGraphForward(
                model, torch.device("cpu"), "hybrid", node_features, edge_features, graph)

    @KratosUnittest.skipUnless(have_cuda, "CUDA is not available.")
    def test_ScalableInterfacesOnCuda(self):
        # regression: every cached tensor (graph, pos, ms_edges AND ms_ids -
        # the model moves ms_edges but not ms_ids) must reach the device
        torch.manual_seed(0)
        self._Save(BiStrideMeshGraphNet(
            input_dim_nodes=1, input_dim_edges=4, output_dim=1, processor_size=2,
            hidden_dim_processor=16, hidden_dim_node_encoder=8, hidden_dim_edge_encoder=8,
            hidden_dim_node_decoder=8, num_layers_bistride=1,
            num_mesh_levels=1, bistride_pos_dim=3))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"   : "Scalable",
                "model_settings"    : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cuda"
                },
                "model_interface"   : "bistride",
                "multiscale_levels" : 1,
                "input_fields"      : [ { "variable_name" : "PRESSURE",
                                          "data_location" : "node_historical" } ],
                "output_fields"     : [ { "variable_name" : "TEMPERATURE",
                                          "data_location" : "node_historical" } ]
            }
        }""" % self.checkpoint)
        process = graph_inference_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self._AssertWritten(numpy.array(
            [node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes]))


have_fluid = KratosUtilities.CheckIfApplicationsAvailable(
    "FluidDynamicsApplication", "LinearSolversApplication")


@KratosUnittest.skipUnless(have_scalable_gnns and have_fluid,
                           "Missing physicsnemo/PyG or FluidDynamics/LinearSolvers.")
class TestScalableGnnOnRealFluidSolve(KratosUnittest.TestCase):
    """The external-aero pairing on a real FluidDynamics solve: a lid-driven
    cavity, its solved velocity/pressure field, and a bistride surrogate
    trained to reproduce the pressure from the velocity."""

    def test_BistrideSurrogateOnCavityFlow(self):
        import fluid_case
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
        model = Kratos.Model()
        analysis = fluid_case.CreateFluidAnalysis(model, lid_velocity=1.0, divisions=8)
        analysis.Run()
        model_part = model["FluidModelPart"]

        # a real solved cavity: the lid drives the flow, the interior circulates
        velocity = numpy.array([[node.GetSolutionStepValue(Kratos.VELOCITY_X),
                                 node.GetSolutionStepValue(Kratos.VELOCITY_Y)]
                                for node in model_part.Nodes])
        pressure = numpy.array([node.GetSolutionStepValue(Kratos.PRESSURE)
                                for node in model_part.Nodes])
        self.assertAlmostEqual(float(numpy.abs(velocity).max()), 1.0, places=6)
        self.assertGreater(float(numpy.abs(velocity[:, 1]).max()), 0.0)
        self.assertGreater(float(pressure.std()), 0.0)
        self.assertEqual(next(iter(model_part.Elements)).Info().split()[0], "VMS")

        node_features, edge_index, edge_features, node_ids = graph_bridge.BuildGraph(
            model_part, (("VELOCITY", "node_historical"),))
        positions = graph_bridge.NodePositions(model_part)
        ms_edges, ms_ids = graph_bridge.BuildBistrideHierarchy(
            edge_index, len(node_ids), positions, num_levels=2)
        self.assertEqual(len(ms_edges), 3)
        self.assertLess(ms_ids[1].size, ms_ids[0].size)

        torch.manual_seed(0)
        surrogate = BiStrideMeshGraphNet(
            input_dim_nodes=3, input_dim_edges=4, output_dim=1, processor_size=2,
            hidden_dim_processor=32, hidden_dim_node_encoder=16, hidden_dim_edge_encoder=16,
            hidden_dim_node_decoder=16, num_layers_bistride=1,
            num_mesh_levels=2, bistride_pos_dim=3)

        graph = graph_bridge.ToPyGGraph(edge_index, len(node_ids), positions)
        graph.pos = graph.pos.float()
        inputs = (torch.tensor(node_features, dtype=torch.float32),
                  torch.tensor(edge_features, dtype=torch.float32))
        multiscale = ([torch.from_numpy(edges) for edges in ms_edges],
                      [torch.from_numpy(ids) for ids in ms_ids])
        target = torch.tensor(pressure, dtype=torch.float32).reshape(-1, 1)

        optimizer = torch.optim.Adam(surrogate.parameters(), lr=5e-3)
        losses = []
        for _ in range(40):
            optimizer.zero_grad()
            prediction = surrogate(inputs[0], inputs[1], graph, *multiscale)
            loss = torch.nn.functional.mse_loss(prediction, target)
            loss.backward()
            optimizer.step()
            losses.append(loss.item())
        self.assertLess(min(losses[1:]), losses[0])

        # and it deploys through the process unchanged
        checkpoint = "test_cavity_bistride.mdlus"
        try:
            surrogate.save(checkpoint)
            # written to the non-historical database: no new historical
            # variable can be added to an already-populated model part
            process = graph_inference_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"   : "FluidModelPart",
                    "model_settings"    : { "checkpoint_file" : "%s",
                                            "checkpoint_type" : "physicsnemo",
                                            "device"          : "cpu" },
                    "model_interface"   : "bistride",
                    "multiscale_levels" : 2,
                    "input_fields"      : [ { "variable_name" : "VELOCITY",
                                              "data_location" : "node_historical" } ],
                    "output_fields"     : [ { "variable_name" : "PRESSURE",
                                              "data_location" : "node_non_historical" } ]
                }
            }""" % checkpoint), model)
            process.ExecuteInitialize()
            model_part.ProcessInfo[Kratos.STEP] = 1
            process.ExecuteFinalizeSolutionStep()
            predicted = numpy.array([node.GetValue(Kratos.PRESSURE)
                                     for node in model_part.Nodes])
            self.assertTrue(numpy.isfinite(predicted).all())
            self.assertGreater(float(predicted.std()), 0.0)
        finally:
            KratosUtilities.DeleteFileIfExisting(checkpoint)


if __name__ == '__main__':
    KratosUnittest.main()
