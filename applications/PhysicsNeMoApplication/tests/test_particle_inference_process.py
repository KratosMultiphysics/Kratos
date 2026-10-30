from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import torch_geometric  # noqa: F401
    import torch_scatter  # noqa: F401
    import physicsnemo.models.meshgraphnet  # noqa: F401
    have_mgn = True
except ImportError:
    have_mgn = False

have_cuda = have_torch and torch.cuda.is_available()


def _CreateParticleModelPart(model, n=12, seed=0, buffer_size=2):
    model_part = model.CreateModelPart("Particles")
    for variable in (Kratos.VELOCITY, Kratos.ACCELERATION, Kratos.DISPLACEMENT):
        model_part.AddNodalSolutionStepVariable(variable)
    model_part.SetBufferSize(buffer_size)
    rng = numpy.random.default_rng(seed)
    for i, xyz in enumerate(rng.random((n, 3))):
        model_part.CreateNewNode(i + 1, *xyz)
    model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.1
    return model_part


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestParticleInferenceProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_particle_stub_model.pt")
        self.model = Kratos.Model()
        self.model_part = _CreateParticleModelPart(self.model)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _CreateProcess(self, extra=""):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import particle_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Particles",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "model_interface" : "tensor",
                "connectivity"    : { "type" : "radius", "radius" : 0.4, "backend" : "numpy" },
                "history_size"    : 2
                %s
            }
        }""" % (self.checkpoint, extra))
        return particle_inference_process.Factory(settings, self.model)

    def test_FreeFallMatchesClosedForm(self):
        # stub predicting constant gravity: after k steps of semi-implicit
        # Euler from rest, v_k = -g k dt and z_k = z0 - g dt^2 k(k+1)/2
        class Gravity(torch.nn.Module):
            def forward(self, nodes, edges, edge_index):
                acceleration = torch.zeros(nodes.shape[0], 3, dtype=nodes.dtype)
                acceleration[:, 2] = -9.81
                return acceleration

        torch.jit.script(Gravity()).save(str(self.checkpoint))
        process = self._CreateProcess()

        z0 = numpy.array([node.Z for node in self.model_part.Nodes])
        g, dt, steps = 9.81, 0.1, 10
        # history_size=2: the first due step only warms the history up
        for step in range(1, steps + 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
            if step == 1:  # warm-up: nothing moved yet
                self.assertTrue(all(node.Z == z for node, z in
                                    zip(self.model_part.Nodes, z0)))

        z = numpy.array([node.Z for node in self.model_part.Nodes])
        vz = numpy.array([
            node.GetSolutionStepValue(Kratos.VELOCITY)[2] for node in self.model_part.Nodes])
        expected_z = z0 - g * dt * dt * steps * (steps + 1) / 2.0
        self.assertTrue(numpy.allclose(vz, -g * dt * steps, atol=1e-10))
        self.assertTrue(numpy.allclose(z, expected_z, atol=1e-9))
        # DISPLACEMENT tracks the motion
        disp_z = numpy.array([
            node.GetSolutionStepValue(Kratos.DISPLACEMENT)[2] for node in self.model_part.Nodes])
        self.assertTrue(numpy.allclose(disp_z, z - z0, atol=1e-12))

    def test_InputNormalizationFromTheCardIsApplied(self):
        """CreateParticleTrajectoryDataset(normalize=True) standardizes the
        FEATURES too, so the network expects standardized velocity
        histories - fed raw, an 18% position drift was measured on the
        Examples case. The stub echoes its last velocity slot as the
        acceleration, so the process must have handed it (v - mean) / std
        exactly, over the whole K*3-wide history."""
        class EchoLastVelocity(torch.nn.Module):
            def forward(self, nodes, edges, edge_index):
                return nodes[:, -3:]

        torch.jit.script(EchoLastVelocity()).save(str(self.checkpoint))
        mean, std = [1.0, 2.0, 3.0] * 2, [2.0, 4.0, 8.0] * 2  # history_size 2 -> width 6
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": mean, "std": std}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting,
                        str(self.checkpoint) + ".card.json")
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, [5.0, 6.0, 7.0])

        process = self._CreateProcess()
        for step in (1, 2):  # step 1 warms the history up
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        acceleration = numpy.array([
            node.GetSolutionStepValue(Kratos.ACCELERATION) for node in self.model_part.Nodes])
        numpy.testing.assert_allclose(
            acceleration, numpy.tile([2.0, 1.0, 0.5], (acceleration.shape[0], 1)), rtol=1e-12)

    def test_NormalizedModelIsDeNormalizedFromTheCard(self):
        """A model trained on standardized accelerations must not be
        deployed raw.

        CreateParticleTrajectoryDataset(normalize=True) standardizes the
        targets, so the network predicts (a - mean) / std. The process
        integrates its output TWICE - v += dt*a, then x += dt*v - so a
        missing de-normalization compounds straight into node positions.
        """
        mean, std = [0.0, 0.0, -9.81], [1.0, 1.0, 2.0]

        class Standardized(torch.nn.Module):
            """Predicts gravity in STANDARDIZED units: (-9.81 + 9.81) / 2 = 0."""

            def forward(self, nodes, edges, edge_index):
                return torch.zeros(nodes.shape[0], 3, dtype=nodes.dtype)

        torch.jit.script(Standardized()).save(str(self.checkpoint))
        model_registry.SaveModelCard(self.checkpoint, {
            "output_normalization": {"type": "mean_std", "mean": mean, "std": std}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting,
                        str(self.checkpoint) + ".card.json")

        process = self._CreateProcess()
        z0 = numpy.array([node.Z for node in self.model_part.Nodes])
        g, dt, steps = 9.81, 0.1, 5
        for step in range(1, steps + 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()

        # de-normalized, the standardized 0 IS -9.81, so free fall follows
        acceleration_z = numpy.array([
            node.GetSolutionStepValue(Kratos.ACCELERATION)[2]
            for node in self.model_part.Nodes])
        self.assertTrue(numpy.allclose(acceleration_z, -g, atol=1e-9))

        vz = numpy.array([node.GetSolutionStepValue(Kratos.VELOCITY)[2]
                          for node in self.model_part.Nodes])
        z = numpy.array([node.Z for node in self.model_part.Nodes])
        self.assertTrue(numpy.allclose(vz, -g * dt * steps, atol=1e-9))
        self.assertTrue(numpy.allclose(
            z, z0 - g * dt * dt * steps * (steps + 1) / 2.0, atol=1e-8))

    def test_WithoutACardTheOutputIsUntouched(self):
        # the identity path: every configuration written before the card
        # carried a normalization must behave exactly as before
        class Constant(torch.nn.Module):
            def forward(self, nodes, edges, edge_index):
                return torch.full((nodes.shape[0], 3), 2.0, dtype=nodes.dtype)

        torch.jit.script(Constant()).save(str(self.checkpoint))
        process = self._CreateProcess()
        for step in (1, 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        acceleration = numpy.array([
            node.GetSolutionStepValue(Kratos.ACCELERATION)[2]
            for node in self.model_part.Nodes])
        self.assertTrue(numpy.allclose(acceleration, 2.0, atol=1e-12))

    def test_NodeTypeOneHot(self):
        # stub echoing the one-hot width through the acceleration: features
        # are (N, 2*3 + 2); return the last one-hot column as a_z
        class Echo(torch.nn.Module):
            def forward(self, nodes, edges, edge_index):
                acceleration = torch.zeros(nodes.shape[0], 3, dtype=nodes.dtype)
                acceleration[:, 2] = nodes[:, -1]
                return acceleration

        torch.jit.script(Echo()).save(str(self.checkpoint))
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.FLAG_VARIABLE, float(node.Id % 2))
        process = self._CreateProcess(extra=""",
                "node_type_variable" : "FLAG_VARIABLE",
                "num_node_types"     : 2""")
        for step in (1, 2):  # step 1 warms the history up
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.ACCELERATION)[2],
                float(node.Id % 2), places=12)

    def test_WrongOutputShapeRaises(self):
        class Flat(torch.nn.Module):
            def forward(self, nodes, edges, edge_index):
                return nodes.sum(dim=-1)

        torch.jit.script(Flat()).save(str(self.checkpoint))
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # warm-up
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        with self.assertRaisesRegex(ValueError, "acceleration"):
            process.ExecuteFinalizeSolutionStep()

    def test_MissingDeltaTimeRaises(self):
        class Gravity(torch.nn.Module):
            def forward(self, nodes, edges, edge_index):
                return torch.zeros(nodes.shape[0], 3, dtype=nodes.dtype)

        torch.jit.script(Gravity()).save(str(self.checkpoint))
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.0
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # warm-up
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        with self.assertRaisesRegex(ValueError, "DELTA_TIME"):
            process.ExecuteFinalizeSolutionStep()


@KratosUnittest.skipUnless(have_torch and have_mgn,
                           "Missing required python modules: torch, torch_geometric, torch_scatter, physicsnemo.")
class TestMeshGraphNetParticles(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_particle_mgn.mdlus")
        self.model = Kratos.Model()
        self.model_part = _CreateParticleModelPart(self.model, n=15, seed=3)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_RealMeshGraphNet(self):
        from physicsnemo.models.meshgraphnet import MeshGraphNet
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import particle_inference_process
        torch.manual_seed(0)
        mgn = MeshGraphNet(
            input_dim_nodes=6, input_dim_edges=4, output_dim=3,
            processor_size=2, hidden_dim_node_encoder=16,
            hidden_dim_edge_encoder=16, hidden_dim_node_decoder=16,
            hidden_dim_processor=16)
        mgn.save(str(self.checkpoint))

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Particles",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "connectivity"    : { "type" : "radius", "radius" : 0.5 },
                "history_size"    : 2
            }
        }""" % self.checkpoint)
        process = particle_inference_process.Factory(settings, self.model)
        for step in (1, 2):  # step 1 warms the history up
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()

        accelerations = numpy.array([
            node.GetSolutionStepValue(Kratos.ACCELERATION)
            for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(accelerations).all())
        self.assertGreater(numpy.abs(accelerations).max(), 0.0)

    @KratosUnittest.skipUnless(have_cuda, "Requires a CUDA device.")
    def test_RealMeshGraphNetOnCuda(self):
        from physicsnemo.models.meshgraphnet import MeshGraphNet
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import particle_inference_process
        # regression: the PyG graph rebuilt every step must follow the
        # model to its resolved device, or the CUDA forward pass fails with
        # a device-mismatch RuntimeError in scatter_add
        torch.manual_seed(0)
        mgn = MeshGraphNet(
            input_dim_nodes=6, input_dim_edges=4, output_dim=3,
            processor_size=2, hidden_dim_node_encoder=16,
            hidden_dim_edge_encoder=16, hidden_dim_node_decoder=16,
            hidden_dim_processor=16)
        mgn.save(str(self.checkpoint))

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Particles",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cuda"
                },
                "connectivity"    : { "type" : "radius", "radius" : 0.5 },
                "history_size"    : 2
            }
        }""" % self.checkpoint)
        process = particle_inference_process.Factory(settings, self.model)
        for step in (1, 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()

        accelerations = numpy.array([
            node.GetSolutionStepValue(Kratos.ACCELERATION)
            for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(accelerations).all())


if __name__ == '__main__':
    KratosUnittest.main()
