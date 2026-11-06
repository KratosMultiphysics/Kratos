import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch  # noqa: F401
    import physicsnemo.sym.eq.phy_informer  # noqa: F401
    import physicsnemo.models.mlp.fully_connected  # noqa: F401
    have_deps = True
except ImportError:
    have_deps = False


def _Harmonic(node):
    return node.X + 2.0 * node.Y - node.Z  # Laplace(u) = 0


@KratosUnittest.skipUnless(have_deps,
                           "Missing required python modules: torch, physicsnemo (.sym).")
class TestPinnSolveProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=3,
            historical_variables=(Kratos.TEMPERATURE, Kratos.NODAL_PAUX))

    def _FixBoundary(self):
        for node in self.model_part.Nodes:
            on_boundary = any(abs(v) < 1e-12 or abs(v - 1.0) < 1e-12
                              for v in (node.X, node.Y, node.Z))
            if on_boundary:
                node.Fix(Kratos.TEMPERATURE)
                node.SetSolutionStepValue(Kratos.TEMPERATURE, _Harmonic(node))

    def test_ForwardLaplaceSolve(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import pinn_solve_process
        self._FixBoundary()
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "mode"            : "forward",
                "physics"         : {
                    "pde"           : "builtin:diffusion",
                    "pde_arguments" : { "D" : 1.0 }
                },
                "solution_fields" : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "network"         : { "layer_size" : 32, "num_layers" : 3 },
                "training"        : {
                    "epochs"          : 400,
                    "learning_rate"   : 5e-3,
                    "physics_weight"  : 1.0,
                    "boundary_weight" : 10.0,
                    "seed"            : 0
                },
                "device"          : "cpu",
                "normalize_coordinates" : false
            }
        }""")
        process = pinn_solve_process.Factory(settings, self.model)
        process.ExecuteBeforeSolutionLoop()

        # the loss must have dropped substantially and the harmonic solution
        # must be reproduced everywhere (the PDE + Dirichlet data determine it)
        self.assertLess(process.loss_history[-1], process.loss_history[0] / 10.0)
        errors = [abs(node.GetSolutionStepValue(Kratos.TEMPERATURE) - _Harmonic(node))
                  for node in self.model_part.Nodes]
        self.assertLess(max(errors), 0.15)

    def test_InverseRecoversDiffusionCoefficient(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import pinn_solve_process
        # observations u = x^2 + y^2 + z^2 with residual -D lap(u) + 6:
        # zero exactly at D = 1
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(
                Kratos.NODAL_PAUX, node.X ** 2 + node.Y ** 2 + node.Z ** 2)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"    : "Main",
                "mode"               : "inverse",
                "physics"            : {
                    "pde"           : "builtin:diffusion",
                    "pde_arguments" : { "D" : null, "source" : -6.0 }
                },
                "inverse_parameters" : { "D" : 0.3 },
                "solution_fields"    : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "observation_fields" : [ { "variable_name" : "NODAL_PAUX",  "data_location" : "node_historical" } ],
                "network"            : { "layer_size" : 32, "num_layers" : 3 },
                "training"           : {
                    "epochs"         : 500,
                    "learning_rate"  : 1e-2,
                    "physics_weight" : 1.0,
                    "data_weight"    : 20.0,
                    "seed"           : 0
                },
                "device"             : "cpu",
                "normalize_coordinates" : false
            }
        }""")
        process = pinn_solve_process.Factory(settings, self.model)
        process.Solve()

        recovered = process.inverse_values["D"]
        self.assertGreater(recovered, 0.5)
        self.assertLess(recovered, 2.0)

    def test_PhysicsGradientsDoNotLeakIntoDetachedField(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import pinn_solve_process
        # inverse mode with data_weight 0: the ONLY loss is the physics term,
        # whose field gradients are detached - so the (seeded) network must
        # stay exactly at its initialization after training steps
        def make_process(epochs):
            settings = Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"    : "Main",
                    "mode"               : "inverse",
                    "physics"            : {
                        "pde"           : "builtin:diffusion",
                        "pde_arguments" : { "D" : null }
                    },
                    "inverse_parameters" : { "D" : 0.3 },
                    "observation_fields" : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_historical" } ],
                    "network"            : { "layer_size" : 8, "num_layers" : 2 },
                    "training"           : {
                        "epochs"        : %d,
                        "data_weight"   : 0.0,
                        "seed"          : 0
                    },
                    "device"             : "cpu"
                }
            }""" % epochs)
            return pinn_solve_process.Factory(settings, self.model)

        frozen = make_process(0)
        frozen.Solve()
        trained = make_process(5)
        trained.Solve()

        probe = torch.rand(10, 3)
        with torch.no_grad():
            self.assertTrue(torch.equal(frozen._network(probe), trained._network(probe)))

    def test_ValidationErrors(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import pinn_solve_process
        with self.assertRaisesRegex(ValueError, "mode"):
            pinn_solve_process.Factory(Kratos.Parameters("""{
                "Parameters": { "model_part_name" : "Main", "mode" : "adjoint" }
            }"""), self.model)
        with self.assertRaisesRegex(ValueError, "inverse_parameters"):
            pinn_solve_process.Factory(Kratos.Parameters("""{
                "Parameters": { "model_part_name" : "Main", "mode" : "inverse" }
            }"""), self.model)
        with self.assertRaisesRegex(ValueError, "autodiff"):
            pinn_solve_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name" : "Main",
                    "physics"         : { "grad_method" : "spectral" }
                }
            }"""), self.model)


def _AnisotropicHarmonic(node):
    # laplacian = 2 + 2 - 4 = 0, and unlike a linear field its SECOND
    # derivatives are non-zero, so a per-axis rescaling of the operator
    # changes the solution
    return node.X ** 2 + node.Y ** 2 - 2.0 * node.Z ** 2


@KratosUnittest.skipUnless(have_deps,
                           "Missing required python modules: torch, physicsnemo (.sym).")
class TestPinnOnANonUnitDomain(KratosUnittest.TestCase):
    """The PDE must be solved in PHYSICAL coordinates.

    normalize_coordinates rescales the network's inputs. If those same
    coordinates are also what the residual is differentiated against, the
    operator becomes sum (1/L_i^2) d2u/dx_i^2 - a different PDE, and a
    per-axis different one when the domain is anisotropic. On a unit cube
    L = 1 hides it completely, which is why every previous test passed.
    """

    _STRETCH = 4.0

    def _StretchedModelPart(self):
        model = Kratos.Model()
        model_part = CreateStructuredTetModelPart(
            model, "Main", divisions=4, historical_variables=(Kratos.TEMPERATURE,))
        for node in model_part.Nodes:
            node.Z0 *= self._STRETCH
            node.Z *= self._STRETCH
        for node in model_part.Nodes:
            on_boundary = (abs(node.X) < 1e-12 or abs(node.X - 1.0) < 1e-12
                           or abs(node.Y) < 1e-12 or abs(node.Y - 1.0) < 1e-12
                           or abs(node.Z) < 1e-12
                           or abs(node.Z - self._STRETCH) < 1e-12)
            if on_boundary:
                node.Fix(Kratos.TEMPERATURE)
                node.SetSolutionStepValue(Kratos.TEMPERATURE, _AnisotropicHarmonic(node))
        return model, model_part

    @staticmethod
    def _Settings(normalize):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "mode"            : "forward",
                "physics"         : {
                    "pde"           : "builtin:diffusion",
                    "pde_arguments" : { "D" : 1.0 }
                },
                "solution_fields" : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "network"         : { "layer_size" : 48, "num_layers" : 4 },
                "training"        : {
                    "epochs"          : 900,
                    "learning_rate"   : 5e-3,
                    "physics_weight"  : 1.0,
                    "boundary_weight" : 10.0,
                    "seed"            : 0
                },
                "device"          : "cpu",
                "normalize_coordinates" : true
            }
        }""")
        settings["Parameters"]["normalize_coordinates"].SetBool(normalize)
        return settings

    def _SolveAndMeasure(self, normalize):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import pinn_solve_process
        model, model_part = self._StretchedModelPart()
        process = pinn_solve_process.Factory(self._Settings(normalize), model)
        process.ExecuteBeforeSolutionLoop()
        interior = [node for node in model_part.Nodes
                    if not node.IsFixed(Kratos.TEMPERATURE)]
        self.assertGreater(len(interior), 0)
        return max(abs(node.GetSolutionStepValue(Kratos.TEMPERATURE)
                       - _AnisotropicHarmonic(node)) for node in interior)

    def test_DefaultNormalizationSolvesThePhysicalPde(self):
        # THE regression guard. With the coordinates handed to the residual
        # normalized, this lands near 4.3 against a field scale of 32.
        error = self._SolveAndMeasure(normalize=True)
        self.assertLess(error, 1.0)

    def test_NormalizationDoesNotChangeTheAnswer(self):
        # normalize_coordinates conditions the network's inputs; it must not
        # change which equation is solved
        normalized = self._SolveAndMeasure(normalize=True)
        raw = self._SolveAndMeasure(normalize=False)
        self.assertLess(normalized, 1.0)
        self.assertLess(raw, 1.0)
        self.assertLess(abs(normalized - raw), 1.0)


if __name__ == '__main__':
    KratosUnittest.main()
