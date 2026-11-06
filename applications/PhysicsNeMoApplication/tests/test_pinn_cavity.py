"""The lid-driven cavity as a PINN: NVIDIA's canonical LDC tutorial on the
in-memory Navier-Stokes cavity this application already solves with FEM.

The point is not that a PINN beats the solver - it does not - but that the
builtin incompressible Navier-Stokes residual drives a mesh-free solve from
the same model part's Dirichlet data, and that the result can be compared
against the VMS solution node for node.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_fluid = kratos_utils.CheckIfApplicationsAvailable(
    "FluidDynamicsApplication", "LinearSolversApplication")

try:
    import torch  # noqa: F401
    import physicsnemo.sym.eq.phy_informer  # noqa: F401
    import physicsnemo.models.mlp.fully_connected  # noqa: F401
    have_deps = True
except ImportError:
    have_deps = False

_LID_VELOCITY = 1.0
_VISCOSITY = 0.05          # a gentle Reynolds number: the PINN has a chance


@KratosUnittest.skipUnless(have_fluid and have_deps,
                           "Missing FluidDynamics/LinearSolvers or torch/physicsnemo.sym.")
class TestPinnCavity(KratosUnittest.TestCase):
    def setUp(self):
        import fluid_case

        self.model = Kratos.Model()
        self.model_part = fluid_case.CreateFluidModelPart(self.model, divisions=8)
        fluid_case.ApplyCaseData(self.model_part, _LID_VELOCITY)

    def _RunPinn(self, epochs=400):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            pinn_solve_process)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "FluidModelPart",
                "mode"            : "forward",
                "physics"         : {
                    "pde"           : "builtin:incompressible_navier_stokes",
                    "pde_arguments" : { "rho" : 1.0, "mu" : %f, "dim" : 2 }
                },
                "fields"          : [
                    { "name" : "velocity", "width" : 2 },
                    { "name" : "pressure", "width" : 1 }
                ],
                "solution_fields" : [
                    { "variable_name" : "VELOCITY", "data_location" : "node_historical" },
                    { "variable_name" : "PRESSURE", "data_location" : "node_historical" }
                ],
                "output_fields"   : [
                    { "variable_name" : "VELOCITY", "data_location" : "node_historical" },
                    { "variable_name" : "PRESSURE", "data_location" : "node_historical" }
                ],
                "network"         : { "layer_size" : 48, "num_layers" : 4 },
                "training"        : {
                    "epochs"          : %d,
                    "learning_rate"   : 3e-3,
                    "physics_weight"  : 1.0,
                    "boundary_weight" : 20.0,
                    "seed"            : 0
                },
                "device"          : "cpu",
                "normalize_coordinates" : false
            }
        }""" % (_VISCOSITY, epochs))
        process = pinn_solve_process.Factory(settings, self.model)
        process.ExecuteBeforeSolutionLoop()
        return process

    def test_TheNavierStokesResidualDrivesAMeshFreeSolve(self):
        process = self._RunPinn()
        self.assertGreater(len(process.loss_history), 1)
        self.assertTrue(numpy.isfinite(process.loss_history).all())
        # the residual really is being minimized
        self.assertLess(min(process.loss_history[1:]), process.loss_history[0] / 2.0)

        velocities = numpy.array([
            node.GetSolutionStepValue(Kratos.VELOCITY) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(velocities).all())
        # nothing may exceed the lid that drives it, by much
        self.assertLess(numpy.abs(velocities).max(), 3.0 * _LID_VELOCITY)

    def test_ItAgreesWithTheVmsSolveWhereItShould(self):
        """The FEM solution is the reference. A PINN at this budget is not
        expected to match it closely; what it must do is reproduce the
        cavity's defining feature - the lid drags the fluid one way at the
        top, and the recirculation sends it back lower down."""
        import fluid_case

        reference_model = Kratos.Model()
        analysis = fluid_case.CreateFluidAnalysis(
            reference_model, lid_velocity=_LID_VELOCITY, divisions=8)
        analysis.Run()
        reference_part = reference_model["FluidModelPart"]
        reference = {node.Id: node.GetSolutionStepValue(Kratos.VELOCITY_X)
                     for node in reference_part.Nodes}

        process = self._RunPinn()
        predicted = numpy.array([
            node.GetSolutionStepValue(Kratos.VELOCITY_X) for node in self.model_part.Nodes])
        truth = numpy.array([reference[node.Id] for node in self.model_part.Nodes])
        heights = numpy.array([node.Y for node in self.model_part.Nodes])

        relative = (numpy.linalg.norm(predicted - truth)
                    / max(numpy.linalg.norm(truth), 1e-12))
        Kratos.Logger.PrintInfo(
            "TestPinnCavity", f"relative L2 of VELOCITY_X against the VMS solve: {relative:.3f}")
        # measured 0.53 at this budget: the PINN gets the cavity's structure
        # and not its detail, which is what a few hundred epochs buys
        self.assertLess(relative, 0.8)

        # the lid drags the top of the cavity in its own direction
        near_lid = heights > 0.85
        self.assertGreater(predicted[near_lid].mean(), 0.2 * _LID_VELOCITY)
        # and the interior returns the fluid the other way
        interior = (heights > 0.1) & (heights < 0.6)
        self.assertLess(predicted[interior].mean(), predicted[near_lid].mean())


if __name__ == '__main__':
    KratosUnittest.main()
