from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.models.transolver
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False

try:
    import warp  # noqa: F401 - FIGConvUNet's neighbor search backend
    import physicsnemo.models.figconvnet.figconvunet
    have_figconvnet = True
except ImportError:
    have_figconvnet = False

try:
    import warnings
    with warnings.catch_warnings():
        # physicsnemo.experimental warns (no API-stability guarantee) on import
        warnings.simplefilter("ignore")
        import physicsnemo.experimental.models.geotransolver  # noqa: F401
    have_geotransolver = True
except ImportError:
    have_geotransolver = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestPointCloudInferenceProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_point_cloud_toy_model.pt")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 10.0 * node.X)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _CreateProcess(self, model_interface="generic", normalize=True, pass_geometry=True):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import point_cloud_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"       : "Main",
                "model_interface"       : "%s",
                "normalize_coordinates" : %s,
                "pass_geometry"         : %s,
                "model_settings"        : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"          : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"         : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % (model_interface, "true" if normalize else "false",
                "true" if pass_geometry else "false", self.checkpoint))
        return point_cloud_inference_process.Factory(settings, self.model)

    def test_GenericInterface(self):
        # model input rows are (x, y, z, PRESSURE); summing them makes the
        # expected nodal output exactly checkable (unit cube: normalized
        # coordinates coincide with the raw ones).
        class SumFeatures(torch.nn.Module):
            def forward(self, x):  # (1, N, 4) -> (1, N, 1)
                return x.sum(dim=-1, keepdim=True)

        torch.jit.script(SumFeatures()).save(str(self.checkpoint))
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                node.X + node.Y + node.Z + 10.0 * node.X, places=6)

    def test_InputNormalizationFromTheCardIsApplied(self):
        # the FIELD features are standardized by the card; the coordinates
        # keep their own convention (normalize_coordinates) and are not
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry

        class SumFeatures(torch.nn.Module):
            def forward(self, x):  # (1, N, 4) -> (1, N, 1)
                return x.sum(dim=-1, keepdim=True)

        torch.jit.script(SumFeatures()).save(str(self.checkpoint))
        mean, std = 4.0, 2.0
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                node.X + node.Y + node.Z + (10.0 * node.X - mean) / std, places=6)

    def test_TransolverInterfaceContract(self):
        # a toy two-argument model exercising the (fx, embedding) call path
        class TwoArgs(torch.nn.Module):
            def forward(self, fx, embedding):  # (1,N,1), (1,N,3) -> (1,N,1)
                return fx + embedding.sum(dim=-1, keepdim=True)

        torch.jit.script(TwoArgs()).save(str(self.checkpoint))
        process = self._CreateProcess(model_interface="transolver")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                10.0 * node.X + node.X + node.Y + node.Z, places=6)

    def test_FigConvNetInterfaceContract(self):
        # a toy two-argument model returning the (features, scalar) tuple
        from typing import Tuple

        class TupleModel(torch.nn.Module):
            def forward(self, vertices, features) -> Tuple[torch.Tensor, torch.Tensor]:
                # (1,N,3), (1,N,1) -> ((1,N,1), (1,1))
                point = features + vertices.sum(dim=-1, keepdim=True)
                return point, torch.ones((1, 1), dtype=features.dtype)

        torch.jit.script(TupleModel()).save(str(self.checkpoint))
        process = self._CreateProcess(model_interface="figconvnet")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                10.0 * node.X + node.X + node.Y + node.Z, places=6)
        self.assertEqual(process.last_scalar_prediction, 1.0)

    def test_FlareInterfaceContract(self):
        # flare shares transolver's (fx, embedding) call contract
        class TwoArgs(torch.nn.Module):
            def forward(self, fx, embedding):  # (1,N,1), (1,N,3) -> (1,N,1)
                return fx + embedding.sum(dim=-1, keepdim=True)

        torch.jit.script(TwoArgs()).save(str(self.checkpoint))
        process = self._CreateProcess(model_interface="flare")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                10.0 * node.X + node.X + node.Y + node.Z, places=6)

    def test_GeoTransolverInterfaceContract(self):
        # a toy model exercising the keyword call path
        # model(local_embedding, local_positions=..., geometry=...)
        from typing import Optional

        class GeoToy(torch.nn.Module):
            def forward(self, local_embedding, local_positions,
                        geometry: Optional[torch.Tensor] = None):
                out = local_embedding + local_positions.sum(dim=-1, keepdim=True)
                if geometry is not None:
                    out = out + geometry.sum(dim=-1, keepdim=True)
                return out

        torch.jit.script(GeoToy()).save(str(self.checkpoint))

        process = self._CreateProcess(model_interface="geotransolver")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                10.0 * node.X + 2.0 * (node.X + node.Y + node.Z), places=6)

        # pass_geometry=false must forward geometry=None
        process = self._CreateProcess(model_interface="geotransolver", pass_geometry=False)
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                10.0 * node.X + node.X + node.Y + node.Z, places=6)

    def test_UnknownInterfaceRaises(self):
        with self.assertRaisesRegex(ValueError, "model interface"):
            self._CreateProcess(model_interface="voxel")

    def test_WrongOutputRankRaises(self):
        class Flat(torch.nn.Module):
            def forward(self, x):
                return x.sum(dim=-1).sum(dim=-1)

        torch.jit.script(Flat()).save(str(self.checkpoint))
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, r"\(1, N, C_out\)"):
            process.ExecuteFinalizeSolutionStep()


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestTransolverThroughProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_point_cloud_transolver.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + 2.0 * node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_RealTransolver(self):
        from physicsnemo.models.transolver import Transolver
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import point_cloud_inference_process
        torch.manual_seed(0)
        transolver = Transolver(
            functional_dim=1, out_dim=1, embedding_dim=3, n_layers=1,
            n_hidden=16, n_head=2, slice_num=4, use_te=False)
        transolver.save(str(self.checkpoint))

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_interface" : "transolver",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % self.checkpoint)
        process = point_cloud_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        values = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)


@KratosUnittest.skipUnless(have_torch and have_geotransolver,
                           "Missing required python modules: torch, physicsnemo (experimental).")
class TestGeoTransolverThroughProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_point_cloud_geotransolver.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + 2.0 * node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_RealGeoTransolver(self):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from physicsnemo.experimental.models.geotransolver import GeoTransolver
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import point_cloud_inference_process
        torch.manual_seed(0)
        # use_te=False: transformer_engine is an optional GPU-only extra
        geotransolver = GeoTransolver(
            functional_dim=1, out_dim=1, geometry_dim=3, n_layers=1,
            n_hidden=16, n_head=2, slice_num=4, use_te=False)
        geotransolver.save(str(self.checkpoint))

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_interface" : "geotransolver",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % self.checkpoint)
        process = point_cloud_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # .mdlus load re-imports the experimental module
            process.ExecuteFinalizeSolutionStep()

        values = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)


@KratosUnittest.skipUnless(have_torch and have_figconvnet,
                           "Missing required python modules: torch, warp, physicsnemo.")
class TestFIGConvUNetThroughProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_point_cloud_figconvunet.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_RealFIGConvUNet(self):
        from physicsnemo.models.figconvnet import FIGConvUNet
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import point_cloud_inference_process
        torch.manual_seed(0)
        figconvunet = FIGConvUNet(
            in_channels=1, out_channels=1, kernel_size=3,
            hidden_channels=[4, 4, 4], num_levels=2, mlp_channels=[8, 8],
            resolution_memory_format_pairs=[("b_xc_y_z", (2, 8, 8)),
                                            ("b_yc_x_z", (8, 2, 8)),
                                            ("b_zc_x_y", (8, 8, 2))],
            has_input_features=True, pooling_layers=[1])
        figconvunet.save(str(self.checkpoint))

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_interface" : "figconvnet",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % self.checkpoint)
        process = point_cloud_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        values = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)
        self.assertIsInstance(process.last_scalar_prediction, float)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestCoordinateNormalizationOffTheUnitCube(KratosUnittest.TestCase):
    """normalize_coordinates on a box that is NOT [0,1]^3.

    Every fixture in this file - and in twenty others - is the unit cube,
    where min-max normalization is exactly the identity and the branch is
    therefore never really executed. A comment in this very file says so.
    That blind spot is how a length-scale bug once survived the whole
    suite, so these run on an anisotropically stretched box.
    """

    _EXTENT = (2.0, 5.0, 0.5)

    def setUp(self):
        self.checkpoint = Path("test_point_cloud_scale_model.pt")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE),
            extent=self._EXTENT)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_EachAxisIsNormalizedIndependently(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.point_cloud_inference_process import (
            GatherPointCloudCoordinates)

        raw = GatherPointCloudCoordinates(self.model_part, normalize=False)
        normalized = GatherPointCloudCoordinates(self.model_part, normalize=True)

        # the raw box really is anisotropic, so this is not a unit cube in
        # disguise and the assertion below has content
        numpy.testing.assert_allclose(raw.max(axis=0), numpy.array(self._EXTENT), rtol=1e-12)
        numpy.testing.assert_allclose(normalized.min(axis=0), [0.0, 0.0, 0.0], atol=1e-12)
        numpy.testing.assert_allclose(normalized.max(axis=0), [1.0, 1.0, 1.0], rtol=1e-12)
        # per axis, not by one global scale
        numpy.testing.assert_allclose(
            normalized, raw / numpy.array(self._EXTENT), rtol=1e-12)

    def test_ADegenerateAxisCollapsesToZeroNotNaN(self):
        # a planar mesh has zero extent on one axis; the guard divides by 1
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.point_cloud_inference_process import (
            GatherPointCloudCoordinates)

        flat = Kratos.Model()
        part = CreateStructuredTetModelPart(
            flat, "Flat", divisions=2, historical_variables=(Kratos.PRESSURE,),
            extent=(1.0, 1.0, 0.0))
        normalized = GatherPointCloudCoordinates(part, normalize=True)
        self.assertTrue(numpy.all(numpy.isfinite(normalized)))
        numpy.testing.assert_allclose(normalized[:, 2], 0.0, atol=1e-12)

    def _RunOn(self, extent, normalize):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import point_cloud_inference_process
        class GeometrySensitive(torch.nn.Module):
            """Reads the coordinates the process feeds it, so the written
            field reports exactly what the model saw."""

            def forward(self, x):          # (1, N, 4) -> (1, N, 1)
                return x[..., :3].sum(dim=-1, keepdim=True)

        torch.jit.script(GeometrySensitive()).save(str(self.checkpoint))

        model = Kratos.Model()
        part = CreateStructuredTetModelPart(
            model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE), extent=extent)
        for node in part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"       : "Main",
                "model_interface"       : "generic",
                "normalize_coordinates" : %s,
                "model_settings"        : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"          : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"         : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % ("true" if normalize else "false", self.checkpoint))
        process = point_cloud_inference_process.Factory(settings, model)
        part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        return numpy.array([node.GetSolutionStepValue(Kratos.TEMPERATURE)
                            for node in part.Nodes])

    def test_ScaleCovariance(self):
        """With normalization ON the answer must not depend on the box size;
        with it OFF it must. Neither direction was asserted anywhere."""
        small = self._RunOn((1.0, 1.0, 1.0), normalize=True)
        large = self._RunOn((2.0, 5.0, 0.5), normalize=True)
        numpy.testing.assert_allclose(small, large, rtol=1e-6)

        raw_small = self._RunOn((1.0, 1.0, 1.0), normalize=False)
        raw_large = self._RunOn((2.0, 5.0, 0.5), normalize=False)
        self.assertGreater(numpy.abs(raw_small - raw_large).max(), 1e-3)


if __name__ == '__main__':
    KratosUnittest.main()
