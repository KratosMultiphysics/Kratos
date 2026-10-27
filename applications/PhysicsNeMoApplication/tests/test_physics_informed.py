import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.physics import physics_informed
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.sym.eq.phy_informer  # noqa: F401 - bundled since 2.1
    have_sym = True
except ImportError:
    have_sym = False


@KratosUnittest.skipUnless(have_torch and have_sym,
                           "Missing required python modules: torch, physicsnemo (.sym).")
class TestPhysicsInformed(KratosUnittest.TestCase):
    def test_ManufacturedSolutionAutodiff(self):
        # u = x^2 + y^2 + z^2 -> -D lap(u) = -6 D, exactly
        informer, names, _ = physics_informed.CreatePhysicsInformer(Kratos.Parameters("""{
            "pde"           : "builtin:diffusion",
            "pde_arguments" : { "D" : 0.5 },
            "grad_method"   : "autodiff"
        }"""))
        self.assertEqual(names, ["diffusion"])

        coordinates = torch.rand(1, 3, 20, requires_grad=True)  # (B, 3, N)
        u = coordinates.square().sum(dim=1, keepdim=True)       # (B, 1, N)
        residual = informer.forward({"u": u, "coordinates": coordinates})["diffusion"]
        self.assertEqual(list(residual.shape), [1, 1, 20])
        self.assertTrue(torch.allclose(residual, torch.full_like(residual, -3.0), atol=1e-5))

    def test_LossTermExactValue(self):
        # a model computing u = x^2 + y^2 + z^2 from its coordinate channels:
        # the diffusion residual is -6 D everywhere -> loss = (6 D)^2 * weight
        class Quadratic(torch.nn.Module):
            def forward(self, x):
                return x[:, :3].square().sum(dim=1, keepdim=True)

        term = physics_informed.MakePhysicsLossTerm(Kratos.Parameters("""{
            "pde"           : "builtin:diffusion",
            "pde_arguments" : { "D" : 0.5 },
            "grad_method"   : "autodiff",
            "weight"        : 2.0
        }"""))
        model = Quadratic()
        inputs = torch.rand(30, 3)
        loss = term(model, inputs, model(inputs))
        self.assertAlmostEqual(float(loss), 2.0 * 9.0, places=4)

    def test_PhysicsTermRegularizesTraining(self):
        # fit the harmonic field u = x + 2y - z twice from the same seed:
        # the physics-informed run must end with a smaller Laplace residual
        # than the data-only run
        inputs = torch.rand(256, 3)
        targets = (inputs * torch.tensor([1.0, 2.0, -1.0])).sum(dim=1, keepdim=True)
        dataset = torch.utils.data.TensorDataset(inputs, targets)
        term = physics_informed.MakePhysicsLossTerm(Kratos.Parameters("""{
            "pde"           : "builtin:diffusion",
            "pde_arguments" : { "D" : 1.0 },
            "grad_method"   : "autodiff",
            "weight"        : 1.0
        }"""))
        training = Kratos.Parameters("""{
            "epochs"        : 40,
            "batch_size"    : 256,
            "learning_rate" : 1e-2,
            "device"        : "cpu",
            "seed"          : 0
        }""")

        residuals = {}
        for label, terms in (("data_only", None), ("physics", [term])):
            torch.manual_seed(0)
            model = torch.nn.Sequential(
                torch.nn.Linear(3, 32), torch.nn.Tanh(), torch.nn.Linear(32, 1))
            training_utils.TrainModel(model, dataset, training.Clone(), extra_loss_terms=terms)
            residuals[label] = float(term(model, inputs, model(inputs)))
        self.assertLess(residuals["physics"], residuals["data_only"])

    def test_LeastSquaresOnMeshGraph(self):
        # a linear field has zero second derivatives: the least-squares
        # residual on the real mesh graph vanishes
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.point_cloud_inference_process import (
            GatherPointCloudCoordinates)

        model = Kratos.Model()
        model_part = CreateStructuredTetModelPart(model, "Main", divisions=2)
        _, edge_index, _, _ = graph_bridge.BuildGraph(model_part)
        coordinates = GatherPointCloudCoordinates(model_part, normalize=False)

        term = physics_informed.MakePhysicsLossTerm(
            Kratos.Parameters("""{
                "pde"           : "builtin:diffusion",
                "pde_arguments" : { "D" : 1.0 },
                "grad_method"   : "least_squares"
            }"""),
            connectivity_provider=lambda: (coordinates, edge_index))

        prediction = torch.tensor(
            coordinates @ numpy.array([[1.0], [2.0], [-1.0]]), dtype=torch.float32)
        loss = term(None, None, prediction)
        self.assertLess(float(loss), 1e-6)

    def test_DiffusionPdeDimension(self):
        """dim=2 must drop the z-Laplacian entirely.

        On a planar collocation cloud (all nodes at z = 0) the 3D operator
        leaves u_zz unconstrained, and a PINN dumps curvature into z to
        cancel the source: the loss converges while the in-plane amplitude
        is wrong by half. A 2D problem needs a 2D operator.
        """
        three_d = str(physics_informed.MakeDiffusionPde(D=1.0, source=0.5).equations)
        self.assertIn("(z, 2)", three_d)
        two_d = str(physics_informed.MakeDiffusionPde(D=1.0, source=0.5, dim=2).equations)
        self.assertNotIn("z", two_d)
        self.assertIn("(x, 2)", two_d)
        self.assertIn("(y, 2)", two_d)
        with self.assertRaisesRegex(ValueError, "dim"):
            physics_informed.MakeDiffusionPde(dim=4)

    def test_BoundaryTrimMakesTheGridResidualExact(self):
        """The upstream FD stencils are wrong on the outermost shell of a
        non-periodic field. A quadratic with -D lap(u) = source EXACTLY
        (FD-exact in the interior) still averaged O(1) over the full grid;
        boundary_trim=1 drops the contaminated shell and the residual
        becomes machine zero, while a wrong field still scores."""
        D, source, spacing = 1.3, 1.0, 1.0 / 15.0
        axis = numpy.linspace(0.0, 1.0, 16)
        X, Y, Z = numpy.meshgrid(axis, axis, axis, indexing="ij")
        field = (X * (1 - X) + Y * (1 - Y) + Z * (1 - Z)) * source / (6.0 * D)
        rows = torch.as_tensor(field, dtype=torch.float64).reshape(-1, 1)

        def MakeTerm(trim):
            return physics_informed.MakePhysicsLossTerm(Kratos.Parameters("""{
                "pde"           : "builtin:diffusion",
                "pde_arguments" : { "D" : %.10f, "source" : %.10f },
                "grad_method"   : "finite_difference",
                "fd_dx"         : %.10f,
                "grid_shape"    : [16, 16, 16],
                "boundary_trim" : %d
            }""" % (D, source, spacing, trim)))

        self.assertGreater(float(MakeTerm(0)(None, None, rows)), 1.0)
        trimmed = MakeTerm(1)
        self.assertLess(float(trimmed(None, None, rows)), 1e-9)
        # a wrong field must still be graded as wrong
        self.assertAlmostEqual(float(trimmed(None, None, 2.0 * rows)), 1.0, places=6)

    def test_BoundaryTrimValidation(self):
        for trim, message in ((-1, "boundary_trim"), (4, "no interior")):
            with self.assertRaisesRegex(ValueError, message):
                physics_informed.MakePhysicsLossTerm(Kratos.Parameters("""{
                    "pde"           : "builtin:diffusion",
                    "grad_method"   : "finite_difference",
                    "grid_shape"    : [8, 8, 8],
                    "boundary_trim" : %d
                }""" % trim))

    def test_FiniteDifferenceGridResidual(self):
        # constant field on a grid: zero Laplacian
        term = physics_informed.MakePhysicsLossTerm(Kratos.Parameters("""{
            "pde"           : "builtin:diffusion",
            "pde_arguments" : { "D" : 1.0 },
            "grad_method"   : "finite_difference",
            "fd_dx"         : 0.1,
            "grid_shape"    : [8, 8, 8]
        }"""))
        prediction = torch.full((8 * 8 * 8, 1), 3.5)
        loss = term(None, None, prediction)
        self.assertAlmostEqual(float(loss), 0.0, places=10)

    def test_ConvectionDiffusionManufacturedSolution(self):
        # The fourth builtin, and the one matching ConvectionDiffusionApplication's
        # stationary residual.  u = x^2 + y^2 + z^2, c = (2, 3, 5), D = 0.1, s = 1:
        #   c . grad(u) = 2(2x) + 3(2y) + 5(2z)
        #   lap(u)      = 6
        # -> residual = 4x + 6y + 10z - 0.1*6 - 1
        informer, names, _ = physics_informed.CreatePhysicsInformer(Kratos.Parameters("""{
            "pde"           : "builtin:convection_diffusion",
            "pde_arguments" : { "D" : 0.1, "cx" : 2.0, "cy" : 3.0, "cz" : 5.0, "source" : 1.0 },
            "grad_method"   : "autodiff"
        }"""))
        self.assertEqual(names, ["convection_diffusion"])
        self.assertEqual(sorted(informer.required_inputs), ["coordinates", "u"])

        coordinates = torch.rand(1, 3, 20, requires_grad=True)  # (B, 3, N)
        residuals = informer.forward({
            "coordinates": coordinates,
            "u": coordinates.square().sum(dim=1, keepdim=True),
        })
        x, y, z = (coordinates[:, i:i + 1, :] for i in range(3))
        expected = 4.0 * x + 6.0 * y + 10.0 * z - 0.1 * 6.0 - 1.0
        self.assertTrue(torch.allclose(
            residuals["convection_diffusion"], expected, atol=1e-4))

    def test_ElasticityManufacturedSolution(self):
        # u = (x^2, y^2, z^2), lambda = mu = 1:
        # div u = 2(x+y+z), d(div u)/dx_i = 2, lap(u_i) = 2
        # -> elasticity_i = -((1+1)*2 + 1*2) = -6, constant
        informer, names, _ = physics_informed.CreatePhysicsInformer(Kratos.Parameters("""{
            "pde"           : "builtin:linear_elasticity",
            "pde_arguments" : { "lmbda" : 1.0, "mu" : 1.0 },
            "grad_method"   : "autodiff"
        }"""))
        self.assertEqual(names, ["elasticity_x", "elasticity_y", "elasticity_z"])
        # pins the component-naming contract with physicsnemo.sym
        self.assertEqual(sorted(informer.required_inputs),
                         ["coordinates", "u_x", "u_y", "u_z"])

        coordinates = torch.rand(1, 3, 20, requires_grad=True)  # (B, 3, N)
        inputs = {"coordinates": coordinates}
        for column, component in enumerate(("u_x", "u_y", "u_z")):
            inputs[component] = coordinates[:, column:column + 1, :].square()
        residuals = informer.forward(inputs)
        for name in names:
            self.assertTrue(torch.allclose(
                residuals[name], torch.full_like(residuals[name], -6.0), atol=1e-4))

    def test_NavierStokesManufacturedSolution(self):
        # v = (y^2, 0, 0), p = x, rho = mu = 1:
        # momentum_x = (v.grad)v_x - lap(v_x) + dp/dx = 0 - 2 + 1 = -1
        # momentum_y = momentum_z = 0; continuity = d(y^2)/dx = 0
        informer, names, _ = physics_informed.CreatePhysicsInformer(Kratos.Parameters("""{
            "pde"           : "builtin:incompressible_navier_stokes",
            "pde_arguments" : { "rho" : 1.0, "mu" : 1.0 },
            "grad_method"   : "autodiff"
        }"""))
        self.assertEqual(
            names, ["momentum_x", "momentum_y", "momentum_z", "continuity"])

        coordinates = torch.rand(1, 3, 20, requires_grad=True)
        # zero components as vanishing QUADRATICS: a linear-in-coords zero
        # loses its autodiff graph at the second-derivative pass
        inputs = {
            "coordinates": coordinates,
            "velocity_x": coordinates[:, 1:2, :].square(),          # y^2
            "velocity_y": 0.0 * coordinates[:, 0:1, :].square(),
            "velocity_z": 0.0 * coordinates[:, 0:1, :].square(),
            "pressure": coordinates[:, 0:1, :] + 0.0 * coordinates[:, 0:1, :].square(),  # x
        }
        residuals = informer.forward(inputs)
        self.assertTrue(torch.allclose(
            residuals["momentum_x"], torch.full_like(residuals["momentum_x"], -1.0),
            atol=1e-4))
        for name in ("momentum_y", "momentum_z", "continuity"):
            self.assertTrue(torch.allclose(
                residuals[name], torch.zeros_like(residuals[name]), atol=1e-4))

    def test_VectorFieldLossTermWithComponents(self):
        # a model computing (x^2, y^2, z^2): each elasticity residual is -6
        # -> loss = weight * sum over the 3 equations of mean(36) = 2*108
        class Quadratic3(torch.nn.Module):
            def forward(self, x):
                return x[:, :3].square()

        term = physics_informed.MakePhysicsLossTerm(Kratos.Parameters("""{
            "pde"           : "builtin:linear_elasticity",
            "pde_arguments" : { "lmbda" : 1.0, "mu" : 1.0 },
            "grad_method"   : "autodiff",
            "weight"        : 2.0,
            "fields"        : [ { "name" : "u", "width" : 3 } ]
        }"""))
        model = Quadratic3()
        inputs = torch.rand(20, 3)
        loss = term(model, inputs, model(inputs))
        self.assertAlmostEqual(float(loss), 2.0 * 3.0 * 36.0, places=2)

    def test_ComponentDefaulting(self):
        settings = Kratos.Parameters("""{
            "fields": [
                { "name" : "u", "width" : 3 },
                { "name" : "p", "width" : 1 },
                { "name" : "s", "width" : 2, "components" : ["s_a", "s_b"] }
            ]
        }""")
        specs = physics_informed._ReadFieldSpecs(settings)
        self.assertEqual(specs, [
            ("u", 3, ["u_x", "u_y", "u_z"]),
            ("p", 1, ["p"]),
            ("s", 2, ["s_a", "s_b"])])

        with self.assertRaisesRegex(ValueError, "components"):
            physics_informed._ReadFieldSpecs(Kratos.Parameters(
                '{"fields": [{"name": "q", "width": 4}]}'))
        with self.assertRaisesRegex(ValueError, "components"):
            physics_informed._ReadFieldSpecs(Kratos.Parameters(
                '{"fields": [{"name": "q", "width": 2, "components": ["only_one"]}]}'))

    def test_ElasticityEnuConversion(self):
        with self.assertRaisesRegex(ValueError, "BOTH"):
            physics_informed.MakeLinearElasticityPde(E=210.0e9)
        pde = physics_informed.MakeLinearElasticityPde(E=2.5, nu=0.25)
        # lmbda = E nu / ((1+nu)(1-2nu)) = 1.0, mu = E / (2(1+nu)) = 1.0
        self.assertEqual(len(pde.equations), 3)

    def test_ValidationErrors(self):
        with self.assertRaisesRegex(ValueError, "builtin"):
            physics_informed.CreatePde("builtin:navier_stokes", Kratos.Parameters("{}"))
        with self.assertRaisesRegex(ValueError, "dotted"):
            physics_informed.CreatePde("NotAPath", Kratos.Parameters("{}"))
        with self.assertRaisesRegex(ValueError, "grad_method"):
            physics_informed.CreatePhysicsInformer(
                Kratos.Parameters('{"grad_method": "symbolic"}'))
        term = physics_informed.MakePhysicsLossTerm(
            Kratos.Parameters('{"grad_method": "least_squares"}'))
        with self.assertRaisesRegex(ValueError, "connectivity_provider"):
            term(None, None, torch.zeros(4, 1))
        term = physics_informed.MakePhysicsLossTerm(
            Kratos.Parameters('{"grad_method": "finite_difference"}'))
        with self.assertRaisesRegex(ValueError, "grid_shape"):
            term(None, None, torch.zeros(4, 1))


if __name__ == '__main__':
    KratosUnittest.main()
