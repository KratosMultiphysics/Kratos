from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import graph_inference_process
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.nn.module.gnn_layers.graph_types import PYG_AVAILABLE
    from physicsnemo.models.meshgraphnet import MeshGraphNet
    have_meshgraphnet = PYG_AVAILABLE
except ImportError:
    have_meshgraphnet = False

have_cuda = have_torch and torch.cuda.is_available()


def _CreateTetPairModelPart(model):
    model_part = model.CreateModelPart("Main")
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
    props = model_part.CreateNewProperties(1)
    for i, xyz in enumerate([
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 1.0)]):
        node = model_part.CreateNewNode(i + 1, *xyz)
        node.SetSolutionStepValue(Kratos.PRESSURE, float(node.Id))
    model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], props)
    model_part.CreateNewElement("Element3D4N", 2, [2, 5, 3, 4], props)
    return model_part


@KratosUnittest.skipUnless(have_meshgraphnet, "Missing required python modules: physicsnemo with torch_geometric/torch_scatter.")
class TestGraphInferenceProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_meshgraphnet.mdlus")
        model = MeshGraphNet(
            input_dim_nodes=1, input_dim_edges=4, output_dim=1,
            processor_size=2,
            hidden_dim_processor=8, hidden_dim_node_encoder=8,
            hidden_dim_edge_encoder=8, hidden_dim_node_decoder=8)
        model.save(str(self.checkpoint))
        self.model = Kratos.Model()
        self.model_part = _CreateTetPairModelPart(self.model)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_MeshGraphNetThroughProcess(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "test_meshgraphnet.mdlus",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        process = graph_inference_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # untrained net: plumbing/shape check

        values = [node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes]
        self.assertTrue(numpy.isfinite(values).all())
        self.assertTrue(any(abs(v) > 0.0 for v in values))

    def _Settings(self, update_edge_features=False):
        return Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "test_meshgraphnet.mdlus",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "update_edge_features" : """ + ("true" if update_edge_features else "false") + """
            }
        }""")

    def test_OutputNormalizationFromTheCardIsApplied(self):
        """Same checkpoint, same mesh, with and without a card: the written
        field must differ by exactly the card's affine map. The raw run is
        the reference because the meshgraphnet interface takes a PyG graph,
        which TorchScript cannot script into an affine stand-in."""
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry

        def Run():
            process = graph_inference_process.Factory(self._Settings(), self.model)
            process.ExecuteInitialize()
            self.model_part.ProcessInfo[Kratos.STEP] = 1
            process.ExecuteFinalizeSolutionStep()
            return numpy.array([
                node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])

        raw = Run()
        mean, std = 5.0, 3.0
        model_registry.SaveModelCard(self.checkpoint, {
            "output_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")
        numpy.testing.assert_allclose(Run(), std * raw + mean, rtol=1e-10, atol=1e-12)

    def test_InputNormalizationFromTheCardIsApplied(self):
        """A card standardizing the raw field must give exactly the output
        of the same checkpoint fed the pre-standardized field without one."""
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        mean, std = 2.0, 4.0

        def Run():
            process = graph_inference_process.Factory(self._Settings(), self.model)
            process.ExecuteInitialize()
            self.model_part.ProcessInfo[Kratos.STEP] = 1
            process.ExecuteFinalizeSolutionStep()
            return numpy.array([
                node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])

        for node in self.model_part.Nodes:  # pre-standardized, no card
            node.SetSolutionStepValue(Kratos.PRESSURE, (float(node.Id) - mean) / std)
        reference = Run()
        for node in self.model_part.Nodes:  # raw field, card does the work
            node.SetSolutionStepValue(Kratos.PRESSURE, float(node.Id))
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")
        numpy.testing.assert_allclose(Run(), reference, rtol=1e-10, atol=1e-12)
        self.assertGreater(numpy.abs(reference).max(), 0.0)

    def test_TheGraphIsNotReExtractedEveryStep(self):
        """The edge set is topology, extracted once in ExecuteInitialize.

        RunInference used to call BuildGraph per step for its node features
        and discard the edge index it had already cached - three orders of
        magnitude more work than the gather it wanted.
        """
        process = graph_inference_process.Factory(self._Settings(), self.model)
        process.ExecuteInitialize()

        calls = []
        original = graph_inference_process.graph_bridge.BuildGraph

        def Counting(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)

        graph_inference_process.graph_bridge.BuildGraph = Counting
        try:
            values = []
            for step in (1, 2, 3):
                self.model_part.ProcessInfo[Kratos.STEP] = step
                process.ExecuteFinalizeSolutionStep()
                values.append([n.GetSolutionStepValue(Kratos.TEMPERATURE)
                               for n in self.model_part.Nodes])
        finally:
            graph_inference_process.graph_bridge.BuildGraph = original

        self.assertEqual(calls, [], "the graph was re-extracted during stepping")
        # a static mesh and static inputs must give a repeatable answer
        numpy.testing.assert_allclose(values[1], values[0], rtol=0.0, atol=0.0)
        numpy.testing.assert_allclose(values[2], values[0], rtol=0.0, atol=0.0)

    def test_MovingNodesUpdateEdgeFeaturesFromTheCachedIndex(self):
        """update_edge_features is geometry, so it must track moved nodes.

        The cached edge *index* is reused; the features are recomputed from
        current coordinates. Serving the cached features here would be
        silently wrong on a deforming mesh.
        """
        process = graph_inference_process.Factory(
            self._Settings(update_edge_features=True), self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        before = [n.GetSolutionStepValue(Kratos.TEMPERATURE) for n in self.model_part.Nodes]

        for node in self.model_part.Nodes:      # deform, keeping topology
            node.X = node.X * 3.0
            node.Y = node.Y * 3.0
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        process.ExecuteFinalizeSolutionStep()
        after = [n.GetSolutionStepValue(Kratos.TEMPERATURE) for n in self.model_part.Nodes]

        self.assertFalse(
            numpy.allclose(after, before),
            "moving the nodes did not change the prediction, so the edge "
            "features were served from the cache instead of recomputed")

    def test_AGrownMeshRebuildsTheEdgeIndex(self):
        """Topology is cached once, but the mesh can change under a run.

        AdaptiveRemeshProcess ships in this application. Without a guard the
        cached edge index keeps describing the old mesh, and the scatter
        either raises or writes the wrong layout.
        """
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "test_meshgraphnet.mdlus",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        process = graph_inference_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        cached = len(process._node_ids)

        # grow the mesh, as a remesh would
        properties = self.model_part.GetProperties()[1]
        next_node = max(node.Id for node in self.model_part.Nodes) + 1
        self.model_part.CreateNewNode(next_node, 0.5, 0.5, 1.0)
        first_two = [node.Id for node in self.model_part.Nodes][:2]
        next_element = max(e.Id for e in self.model_part.Elements) + 1
        self.model_part.CreateNewElement(
            "Element2D3N", next_element, first_two + [next_node], properties)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0)

        self.model_part.ProcessInfo[Kratos.STEP] = 2
        process.ExecuteFinalizeSolutionStep()

        self.assertEqual(len(process._node_ids), cached + 1)
        values = [node.GetSolutionStepValue(Kratos.TEMPERATURE)
                  for node in self.model_part.Nodes]
        self.assertEqual(len(values), cached + 1)
        self.assertTrue(numpy.isfinite(values).all())

    @KratosUnittest.skipUnless(have_cuda, "Requires a CUDA device.")
    def test_MeshGraphNetThroughProcessOnCuda(self):
        # regression: the cached PyG graph (built once in RunInference) must
        # follow the model to its resolved device, or the first CUDA forward
        # pass fails with a device-mismatch RuntimeError in scatter_add
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "test_meshgraphnet.mdlus",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cuda"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        process = graph_inference_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        values = [node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes]
        self.assertTrue(numpy.isfinite(values).all())


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestGraphBridgePyGConversion(KratosUnittest.TestCase):
    def test_ToPyGGraph(self):
        try:
            import torch_geometric  # noqa: F401
        except ImportError:
            self.skipTest("Missing required python module: torch_geometric.")
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
        model = Kratos.Model()
        model_part = _CreateTetPairModelPart(model)
        _, edge_index, _, node_ids = graph_bridge.BuildGraph(model_part)
        graph = graph_bridge.ToPyGGraph(edge_index, len(node_ids))
        self.assertEqual(graph.num_nodes, 5)
        self.assertEqual(tuple(graph.edge_index.shape), (2, 18))


if __name__ == '__main__':
    KratosUnittest.main()
