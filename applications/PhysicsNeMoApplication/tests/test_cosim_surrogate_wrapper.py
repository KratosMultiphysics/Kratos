from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

have_cosim = KratosUtilities.CheckIfApplicationsAvailable("CoSimulationApplication")
have_mapping = KratosUtilities.CheckIfApplicationsAvailable("MappingApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

_FIXTURE_DIR = Path(__file__).parent / "cosim_cases"
_MDPA_FILE = str(_FIXTURE_DIR / "surrogate_interface")


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


def _WrapperSettings(checkpoint, time_step=0.0):
    return Kratos.Parameters("""{
        "type" : "KratosMultiphysics.PhysicsNeMoApplication.deployment.cosim_surrogate_solver_wrapper",
        "solver_wrapper_settings" : {
            "mdpa_file"      : "%s",
            "time_step"      : %f,
            "model_settings" : { "checkpoint_file" : "%s", "device" : "cpu" },
            "input_fields"   : [ { "variable_name" : "FORCE",        "data_location" : "node_historical" } ],
            "output_fields"  : [ { "variable_name" : "DISPLACEMENT", "data_location" : "node_historical" } ]
        },
        "data" : {
            "load" : { "model_part_name" : "Surrogate", "variable_name" : "FORCE",        "dimension" : 3 },
            "disp" : { "model_part_name" : "Surrogate", "variable_name" : "DISPLACEMENT", "dimension" : 3 }
        }
    }""" % (_MDPA_FILE, time_step, checkpoint))


@KratosUnittest.skipUnless(have_cosim, "Missing required application: CoSimulationApplication.")
@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestCoSimSurrogateWrapper(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_cosim_surrogate_affine.pt")
        _SaveAffineModel(self.checkpoint, 2.0, 1.0)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_FactoryResolvesDottedPath(self):
        # proves the solver-wrapper factory's PYTHONPATH fallback finds the
        # wrapper by its full module path
        from KratosMultiphysics.CoSimulationApplication.factories import solver_wrapper_factory

        wrapper = solver_wrapper_factory.CreateSolverWrapper(
            _WrapperSettings(self.checkpoint), None, "surrogate")
        self.assertEqual(type(wrapper).__name__, "CoSimSurrogateSolverWrapper")
        self.assertEqual(wrapper.model_part.NumberOfNodes(), 8)

    def test_StandaloneSolve(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import cosim_surrogate_solver_wrapper
        wrapper = cosim_surrogate_solver_wrapper.Create(
            _WrapperSettings(self.checkpoint, time_step=1.0), None, "surrogate")
        wrapper.Initialize()
        wrapper.Check()

        for node in wrapper.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.FORCE, [node.X, node.Y, node.Z])

        new_time = wrapper.AdvanceInTime(0.0)
        self.assertAlmostEqual(new_time, 1.0)
        wrapper.SolveSolutionStep()
        wrapper.Finalize()

        for node in wrapper.model_part.Nodes:
            disp = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            for value, coord in zip(disp, (node.X, node.Y, node.Z)):
                self.assertAlmostEqual(value, 2.0 * coord + 1.0, places=6)

    def test_OutputNormalizationFromTheCardIsApplied(self):
        """The wrapper writes through WriteOutputFields like every other
        deployment path, and the docs listed it as covered - but it passed
        no normalization, so a card was silently ignored here. Affine
        stand-in, so a missing shift cannot hide."""
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import cosim_surrogate_solver_wrapper
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        mean, std = [10.0, 20.0, 30.0], [2.0, 3.0, 4.0]
        model_registry.SaveModelCard(self.checkpoint, {
            "output_normalization": {"type": "mean_std", "mean": mean, "std": std}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")

        wrapper = cosim_surrogate_solver_wrapper.Create(
            _WrapperSettings(self.checkpoint, time_step=1.0), None, "surrogate")
        wrapper.Initialize()
        for node in wrapper.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.FORCE, [node.X, node.Y, node.Z])
        wrapper.AdvanceInTime(0.0)
        wrapper.SolveSolutionStep()
        wrapper.Finalize()

        for node in wrapper.model_part.Nodes:
            disp = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            for axis, coord in enumerate((node.X, node.Y, node.Z)):
                self.assertAlmostEqual(
                    disp[axis], (2.0 * coord + 1.0) * std[axis] + mean[axis], places=6)

    def test_InputNormalizationFromTheCardIsApplied(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import cosim_surrogate_solver_wrapper
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        mean, std = [1.0, 2.0, 3.0], [2.0, 4.0, 8.0]
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": mean, "std": std}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")

        wrapper = cosim_surrogate_solver_wrapper.Create(
            _WrapperSettings(self.checkpoint, time_step=1.0), None, "surrogate")
        wrapper.Initialize()
        for node in wrapper.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.FORCE, [node.X, node.Y, node.Z])
        wrapper.AdvanceInTime(0.0)
        wrapper.SolveSolutionStep()
        wrapper.Finalize()

        for node in wrapper.model_part.Nodes:
            disp = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            for axis, coord in enumerate((node.X, node.Y, node.Z)):
                self.assertAlmostEqual(
                    disp[axis], 2.0 * (coord - mean[axis]) / std[axis] + 1.0, places=6)

    def test_DrivenWrapperDoesNotOwnTime(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import cosim_surrogate_solver_wrapper
        wrapper = cosim_surrogate_solver_wrapper.Create(
            _WrapperSettings(self.checkpoint), None, "surrogate")
        self.assertEqual(wrapper.AdvanceInTime(5.0), 0.0)

    def test_UnknownInterfaceRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import cosim_surrogate_solver_wrapper
        settings = _WrapperSettings(self.checkpoint)
        settings["solver_wrapper_settings"].AddString("model_interface", "voxel")
        with self.assertRaisesRegex(ValueError, "model interface"):
            cosim_surrogate_solver_wrapper.Create(settings, None, "surrogate")


@KratosUnittest.skipUnless(have_cosim and have_mapping,
                           "Missing required applications: CoSimulationApplication, MappingApplication.")
@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestCoSimSurrogateCoupledLoop(KratosUnittest.TestCase):
    def setUp(self):
        self.affine = Path("test_cosim_loop_affine.pt")
        self.identity = Path("test_cosim_loop_identity.pt")
        self.contraction = Path("test_cosim_loop_contraction.pt")
        _SaveAffineModel(self.affine, 2.0, 1.0)
        _SaveAffineModel(self.identity, 1.0, 0.0)
        _SaveAffineModel(self.contraction, 0.5, 1.0)

    def tearDown(self):
        for path in (self.affine, self.identity, self.contraction):
            KratosUtilities.DeleteFileIfExisting(str(path))

    @staticmethod
    def _SolverBlock(checkpoint, time_step=0.0):
        # embeddable version of _WrapperSettings
        return _WrapperSettings(checkpoint, time_step)

    def test_WeakCouplingLoop(self):
        # predictor (owns time): DISP = 2 F + 1; transfer DISP -> receiver's
        # F (nearest-neighbor on the identical cube); receiver: DISP = F
        from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

        parameters = Kratos.Parameters("""{
            "problem_data" : {
                "start_time"    : 0.0,
                "end_time"      : 2.0,
                "echo_level"    : 0,
                "print_colors"  : false,
                "parallel_type" : "OpenMP"
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
        parameters["solver_settings"]["solvers"].AddValue(
            "predictor", self._SolverBlock(self.affine, time_step=1.0))
        parameters["solver_settings"]["solvers"].AddValue(
            "receiver", self._SolverBlock(self.identity))

        analysis = CoSimulationAnalysis(parameters)
        analysis.Initialize()

        wrappers = analysis._GetSolver().solver_wrappers
        for node in wrappers["predictor"].model_part.Nodes:
            node.SetSolutionStepValue(Kratos.FORCE, [node.X, 2.0 * node.Y, 3.0 * node.Z])

        analysis.RunSolutionLoop()
        analysis.Finalize()

        for node in wrappers["receiver"].model_part.Nodes:
            disp = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            expected = (2.0 * node.X + 1.0, 2.0 * 2.0 * node.Y + 1.0, 2.0 * 3.0 * node.Z + 1.0)
            for value, ref in zip(disp, expected):
                self.assertAlmostEqual(value, ref, places=6)

    def test_StrongCouplingAitken(self):
        # fixed-point loop d = 0.5 d + 1 (solution d = 2) iterated by
        # gauss_seidel_strong with an Aitken accelerator (a uniform scalar
        # contraction gives collinear residuals, which is exactly Aitken's
        # regime - MVQN's secant Jacobian would be singular here)
        from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

        parameters = Kratos.Parameters("""{
            "problem_data" : {
                "start_time"    : 0.0,
                "end_time"      : 1.0,
                "echo_level"    : 0,
                "print_colors"  : false,
                "parallel_type" : "OpenMP"
            },
            "solver_settings" : {
                "type"                    : "coupled_solvers.gauss_seidel_strong",
                "echo_level"              : 0,
                "num_coupling_iterations" : 15,
                "convergence_accelerators" : [{
                    "type"      : "aitken",
                    "solver"    : "left",
                    "data_name" : "load"
                }],
                "convergence_criteria" : [{
                    "type"          : "relative_norm_previous_residual",
                    "solver"        : "left",
                    "data_name"     : "load",
                    "abs_tolerance" : 1e-10,
                    "rel_tolerance" : 1e-10
                }],
                "data_transfer_operators" : {
                    "mapper" : {
                        "type"            : "kratos_mapping",
                        "mapper_settings" : { "mapper_type" : "nearest_neighbor" }
                    }
                },
                "coupling_sequence" : [
                    {
                        "name"            : "left",
                        "input_data_list" : [{
                            "data"                   : "load",
                            "from_solver"            : "right",
                            "from_solver_data"       : "disp",
                            "data_transfer_operator" : "mapper"
                        }],
                        "output_data_list" : []
                    },
                    {
                        "name"            : "right",
                        "input_data_list" : [{
                            "data"                   : "load",
                            "from_solver"            : "left",
                            "from_solver_data"       : "disp",
                            "data_transfer_operator" : "mapper"
                        }],
                        "output_data_list" : []
                    }
                ],
                "solvers" : {}
            }
        }""")
        parameters["solver_settings"]["solvers"].AddValue(
            "left", self._SolverBlock(self.contraction, time_step=1.0))
        parameters["solver_settings"]["solvers"].AddValue(
            "right", self._SolverBlock(self.identity))

        analysis = CoSimulationAnalysis(parameters)
        analysis.Initialize()

        # a nonzero start makes the first coupling residual nonzero (an
        # all-zero state trivially satisfies the relative criterion)
        wrappers = analysis._GetSolver().solver_wrappers
        for node in wrappers["right"].model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [5.0, 5.0, 5.0])

        analysis.RunSolutionLoop()
        analysis.Finalize()

        for node in wrappers["left"].model_part.Nodes:
            disp = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            for value in disp:
                self.assertAlmostEqual(value, 2.0, places=3)


if __name__ == '__main__':
    KratosUnittest.main()
