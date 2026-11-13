"""Integration tests driven by RomApplication's REAL basis machinery.

The neural-augmented-reduced-basis pattern end to end: thermal solves for a
conductivity sweep feed RomApplication's CalculateRomBasisOutputProcess (the
actual SVD/basis writer), rom_bridge consumes the produced basis, and a
torch surrogate maps the conductivity to the reduced coordinates, deployed
by RomSurrogateProcess at an unseen parameter. Availability-gated on the
compiled applications.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges import rom_bridge
have_rom = kratos_utils.CheckIfApplicationsAvailable(
    "RomApplication", "ConvectionDiffusionApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

_TRAIN_CONDUCTIVITIES = (0.5, 1.0, 2.0, 4.0)
_DIVISIONS = 8


def _RunRealSolve(conductivity):
    import KratosMultiphysics.ConvectionDiffusionApplication  # noqa: F401 - registers the variables
    import thermal_case
    model = Kratos.Model()
    analysis = thermal_case.CreateThermalAnalysis(
        model, conductivity=conductivity, heat_flux=1.0, divisions=_DIVISIONS)
    analysis.Run()
    return model, model["ThermalModelPart"]


def _GenerateRealBasis(output_folder):
    """Drives the real CalculateRomBasisOutputProcess over the sweep.

    A persistent model part receives each solve's TEMPERATURE (copied via
    the same VariableUtils calls the process uses); PrintOutput() records
    each snapshot and ExecuteFinalize() runs the SVD and writes the basis.
    Returns (host_model, host_part, snapshots (n_dofs, S)).
    """
    import KratosMultiphysics.ConvectionDiffusionApplication  # noqa: F401 - registers the variables
    import thermal_case
    from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import (
        CalculateRomBasisOutputProcess)

    host_model = Kratos.Model()
    host_part = thermal_case.CreateThermalModelPart(host_model, _DIVISIONS)

    settings = Kratos.Parameters("""{
        "model_part_name"          : "ThermalModelPart",
        "nodal_unknowns"           : ["TEMPERATURE"],
        "rom_basis_output_format"  : "numpy",
        "rom_basis_output_name"    : "RomParameters",
        "rom_basis_output_folder"  : "%s",
        "svd_truncation_tolerance" : 1e-8,
        "print_singular_values"    : true
    }""" % output_folder)
    process = CalculateRomBasisOutputProcess(host_model, settings)

    snapshots = []
    for conductivity in _TRAIN_CONDUCTIVITIES:
        _, solved_part = _RunRealSolve(conductivity)
        values = Kratos.VariableUtils().GetSolutionStepValuesVector(
            solved_part.Nodes, Kratos.TEMPERATURE, 0)
        Kratos.VariableUtils().SetSolutionStepValuesVector(
            host_part.Nodes, Kratos.TEMPERATURE, values, 0)
        snapshots.append(numpy.array(values))
        process.PrintOutput()
    process.ExecuteFinalize()
    return host_model, host_part, numpy.stack(snapshots, axis=1)


@KratosUnittest.skipUnless(have_rom,
                           "Missing required applications: RomApplication, "
                           "ConvectionDiffusionApplication, LinearSolversApplication.")
class TestRomApplicationBasisInterop(KratosUnittest.TestCase):
    def setUp(self):
        self.basis_folder = Path("test_rom_real_basis")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting(str(self.basis_folder))

    def test_RealBasisRoundTrip(self):
        host_model, host_part, snapshots = _GenerateRealBasis(self.basis_folder)

        basis = rom_bridge.LoadRomBasis(self.basis_folder)
        self.assertEqual(basis.nodal_unknowns, ("TEMPERATURE",))
        self.assertEqual(basis.n_nodes, host_part.NumberOfNodes())
        self.assertIsNotNone(basis.singular_values)
        # SVD basis: orthonormal columns
        numpy.testing.assert_allclose(
            basis.phi.T @ basis.phi, numpy.eye(basis.n_modes), atol=1e-10)
        # NodeIds match the mesh
        self.assertEqual(
            sorted(basis.node_ids.tolist()), sorted(node.Id for node in host_part.Nodes))
        # every training snapshot reconstructs through the basis
        for column in range(snapshots.shape[1]):
            u = snapshots[:, column]
            reconstructed = rom_bridge.ReconstructFromReducedSpace(
                basis, rom_bridge.ProjectToReducedSpace(basis, u))
            self.assertLess(
                numpy.linalg.norm(reconstructed - u) / numpy.linalg.norm(u), 1e-6)
        # gather from the host part equals the last written snapshot
        numpy.testing.assert_allclose(
            rom_bridge.GatherUnknownsVector(host_part, basis), snapshots[:, -1], atol=1e-12)


@KratosUnittest.skipUnless(have_rom and have_torch,
                           "Missing required applications/modules: RomApplication, "
                           "ConvectionDiffusionApplication, LinearSolversApplication, torch.")
class TestRomSurrogateOnRealSolves(KratosUnittest.TestCase):
    def setUp(self):
        self.basis_folder = Path("test_rom_real_surrogate_basis")
        self.checkpoint = Path("test_rom_real_surrogate.pt")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting(str(self.basis_folder))
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def test_ParameterToModesSurrogate(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import rom_surrogate_process
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
        _, _, snapshots = _GenerateRealBasis(self.basis_folder)
        basis = rom_bridge.LoadRomBasis(self.basis_folder)

        # train conductivity -> q on the real solves (TrainModel's contract fits)
        q_targets = rom_bridge.ProjectToReducedSpace(basis, snapshots).T  # (S, n_modes)
        inputs = torch.tensor([[k] for k in _TRAIN_CONDUCTIVITIES], dtype=torch.float64)
        targets = torch.tensor(q_targets, dtype=torch.float64)
        dataset = torch.utils.data.TensorDataset(inputs, targets)
        torch.manual_seed(0)  # deterministic weight init regardless of suite order
        surrogate = torch.nn.Sequential(
            torch.nn.Linear(1, 16), torch.nn.Tanh(), torch.nn.Linear(16, basis.n_modes)).double()
        history = training_utils.TrainModel(surrogate, dataset, Kratos.Parameters("""{
            "epochs"        : 2000,
            "batch_size"    : 4,
            "learning_rate" : 5e-3,
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertLess(history[-1], history[0])
        torch.jit.script(surrogate).save(str(self.checkpoint))

        # deploy at a held-out conductivity on a REAL solve's mesh
        test_conductivity = 1.5
        model, model_part = _RunRealSolve(test_conductivity)
        reference = numpy.array(Kratos.VariableUtils().GetSolutionStepValuesVector(
            model_part.Nodes, Kratos.TEMPERATURE, 0))

        carrier = Kratos.KratosGlobals.GetVariable("PROJECTED_SCALAR1")  # ConvDiff-app variable
        for node in model_part.Nodes:
            node.SetSolutionStepValue(carrier, test_conductivity)

        process = rom_surrogate_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"  : "ThermalModelPart",
                "rom_basis_folder" : "%s",
                "model_settings"   : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"     : [ { "variable_name" : "PROJECTED_SCALAR1", "data_location" : "node_historical" } ]
            }
        }""" % (self.basis_folder, self.checkpoint)), model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        predicted = numpy.array(Kratos.VariableUtils().GetSolutionStepValuesVector(
            model_part.Nodes, Kratos.TEMPERATURE, 0))
        rel_error = numpy.linalg.norm(predicted - reference) / numpy.linalg.norm(reference)
        # smooth 1D k->q map + 4 training solves: a loose but meaningful bound
        self.assertLess(rel_error, 0.05)


if __name__ == '__main__':
    KratosUnittest.main()
