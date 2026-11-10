import importlib

import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
if kratos_utilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    import KratosMultiphysics.ConvectionDiffusionApplication
if kratos_utilities.CheckIfApplicationsAvailable("ConstitutiveLawsApplication"):
    import KratosMultiphysics.ConstitutiveLawsApplication

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(hasattr(KratosMultiphysics, "Future"), "Kratos is not compiled with KRATOS_USE_FUTURE")
class TestFutureRom(KratosUnittest.TestCase):
    """Checks the ROM based on the Future schemes against the legacy ROM builder and solvers."""

    def _RunLegacyRom(self, parameters_filename):
        with open(parameters_filename,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        simulation = rom_testing_utilities.SetUpSimulationInstance(KratosMultiphysics.Model(), parameters)
        simulation.Run()
        return simulation._GetSolver().GetComputingModelPart()

    def _RunFutureRom(self, parameters_filename, projection_backend, max_iterations=None, decoder_factory=None):
        from KratosMultiphysics.RomApplication.future.future_rom_analysis import CreateFutureRomAnalysisInstance

        with open(parameters_filename,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        analysis_stage_module_name = parameters["analysis_stage"].GetString()
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_module_name.split('.')[-1].split('_'))
        analysis_stage_class = getattr(importlib.import_module(analysis_stage_module_name), analysis_stage_class_name)

        future_rom_settings = KratosMultiphysics.Parameters("""{}""")
        future_rom_settings.AddString("projection_backend", projection_backend)
        if max_iterations is not None: # Otherwise taken from the solver settings of the case
            future_rom_settings.AddInt("max_iterations", max_iterations)
        simulation = CreateFutureRomAnalysisInstance(analysis_stage_class, KratosMultiphysics.Model(), parameters, future_rom_settings, decoder_factory=decoder_factory)
        simulation.Run()
        return simulation._GetSolver().GetComputingModelPart()

    def _CheckResults(self, obtained_output, expected_output, relative_tolerance):
        obtained_output = np.array(obtained_output)
        expected_output = np.array(expected_output)
        error = np.linalg.norm(obtained_output - expected_output) / np.linalg.norm(expected_output)
        self.assertLess(error, relative_tolerance)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticRom2D(self):
        with KratosUnittest.WorkFolderScope("structural_static_test_files/ROM/", __file__):
            legacy_model_part = self._RunLegacyRom("../ProjectParameters.json")
            expected_output = rom_testing_utilities.GetVectorNodalResults(legacy_model_part, KratosMultiphysics.DISPLACEMENT)
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    future_model_part = self._RunFutureRom("../ProjectParameters.json", projection_backend)
                    obtained_output = rom_testing_utilities.GetVectorNodalResults(future_model_part, KratosMultiphysics.DISPLACEMENT)
                    # Linear analysis: one iteration per step on the moved mesh, as the legacy linear strategy
                    self._CheckResults(obtained_output, expected_output, 1.0e-12)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticTotalLagrangianRom2D(self):
        # Geometrically nonlinear cantilever (total Lagrangian elements, tip displacement 27% of the length). These elements
        # use the current coordinates, so the scheme must move the mesh as the legacy solver does ("move_mesh_flag")
        with KratosUnittest.WorkFolderScope("future_rom_test_files/TotalLagrangianROM/", __file__):
            legacy_model_part = self._RunLegacyRom("ProjectParameters.json")
            expected_output = rom_testing_utilities.GetVectorNodalResults(legacy_model_part, KratosMultiphysics.DISPLACEMENT)
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    future_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, max_iterations=20)
                    obtained_output = rom_testing_utilities.GetVectorNodalResults(future_model_part, KratosMultiphysics.DISPLACEMENT)
                    self._CheckResults(obtained_output, expected_output, 1.0e-8)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    @KratosUnittest.skipUnless(have_torch, "Missing required python module: PyTorch.")
    def testStructuralStaticTotalLagrangianTorchDecoder2D(self):
        # The linear decoder written as a PyTorch function: its Jacobian, obtained by automatic differentiation in every
        # iteration, must give the solution of the LinearDecoder
        with KratosUnittest.WorkFolderScope("future_rom_test_files/TotalLagrangianROM/", __file__):
            def TorchLinearDecoderFactory(rom_data_folder, dof_set):
                from KratosMultiphysics.RomApplication.future.rom_decoders import LinearDecoder, TorchDecoder
                linear_decoder = LinearDecoder.FromRomData(rom_data_folder, dof_set)
                phi, u0 = torch.as_tensor(linear_decoder.phi), torch.as_tensor(linear_decoder.u0)
                return TorchDecoder(lambda q: u0 + phi @ q, phi.shape[1])
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    linear_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend)
                    expected_output = rom_testing_utilities.GetVectorNodalResults(linear_model_part, KratosMultiphysics.DISPLACEMENT)
                    torch_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, decoder_factory=TorchLinearDecoderFactory)
                    obtained_output = rom_testing_utilities.GetVectorNodalResults(torch_model_part, KratosMultiphysics.DISPLACEMENT)
                    self._CheckResults(obtained_output, expected_output, 1.0e-12)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "ConstitutiveLawsApplication")
    def testStructuralStaticJ2PlasticityRom2D(self):
        # Material nonlinearity: J2 plasticity with linear hardening, half of the elements yield. The internal variables of the
        # constitutive law are updated at the end of each step by the Future scheme, as by the legacy one
        with KratosUnittest.WorkFolderScope("future_rom_test_files/J2PlasticityROM/", __file__):
            legacy_model_part = self._RunLegacyRom("ProjectParameters.json")
            expected_output = rom_testing_utilities.GetVectorNodalResults(legacy_model_part, KratosMultiphysics.DISPLACEMENT)
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    future_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend)
                    obtained_output = rom_testing_utilities.GetVectorNodalResults(future_model_part, KratosMultiphysics.DISPLACEMENT)
                    self._CheckResults(obtained_output, expected_output, 1.0e-8)

    @KratosUnittest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
    def testConvDiffStationaryRom2D(self):
        with KratosUnittest.WorkFolderScope("thermal_static_test_files/ROM/", __file__):
            legacy_model_part = self._RunLegacyRom("../ProjectParameters.json")
            expected_output = rom_testing_utilities.GetScalarNodalResults(legacy_model_part, KratosMultiphysics.TEMPERATURE)
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    # The legacy ROM stops after the first iteration of this non-linear case, as its residual criterion gets a null residual
                    future_model_part = self._RunFutureRom("../ProjectParameters.json", projection_backend, max_iterations=1)
                    obtained_output = rom_testing_utilities.GetScalarNodalResults(future_model_part, KratosMultiphysics.TEMPERATURE)
                    self._CheckResults(obtained_output, expected_output, 1.0e-10)


    @KratosUnittest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
    def testConvDiffStationaryRomFixedDofsDecoding(self):
        # With "fixed_dofs_decoding": "prescribed" the decoder gives the prescribed values on the fixed DOFs (the imposed temperature
        # of the left wall) and the ROM solution on the free ones. The solution is the one of the default decoding ("basis"), as the
        # projection does not use the rows of the fixed DOFs. Note that the basis of this test includes the fixed DOFs
        from KratosMultiphysics.RomApplication.future.future_rom_analysis import CreateFutureRomAnalysisInstance
        from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

        class FixedDofsCheck(ConvectionDiffusionAnalysis):
            def OutputSolutionStep(self):
                super().OutputSolutionStep()
                solver = self.future_rom_solver
                if not hasattr(solver.decoder, "fixed_dofs"): # Only the FixedDofsDecoder ("prescribed")
                    return
                decoded = solver.decoder.Decode(solver.GetRomSolution())
                # The fixed DOFs of the step, read by the decoder once the conditions were applied (the processes free them when finalizing)
                self.decoded = {(dof.Id(), bool(solver.decoder.fixed_dofs[i])): decoded[i] for i, dof in enumerate(solver.GetEffectiveDofSet())}

        def Run(fixed_dofs_decoding):
            with open("../ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            settings = KratosMultiphysics.Parameters('{"max_iterations" : 1}')
            settings.AddString("fixed_dofs_decoding", fixed_dofs_decoding)
            simulation = CreateFutureRomAnalysisInstance(FixedDofsCheck, KratosMultiphysics.Model(), parameters, settings)
            simulation.Run()
            return simulation

        with KratosUnittest.WorkFolderScope("thermal_static_test_files/ROM/", __file__):
            simulation = Run("prescribed")
            model_part = simulation._GetSolver().GetComputingModelPart()
            fixed = {node_id: value for (node_id, is_fixed), value in simulation.decoded.items() if is_fixed}
            free = {node_id: value for (node_id, is_fixed), value in simulation.decoded.items() if not is_fixed}
            self.assertGreater(len(fixed), 0)
            for node_id, value in fixed.items():
                self.assertEqual(value, 303.15)
            for node_id, value in free.items():
                self.assertAlmostEqual(value, model_part.GetNode(node_id).GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), delta=1.0e-10 * abs(value))

            default_model_part = Run("basis")._GetSolver().GetComputingModelPart()
            obtained_output = rom_testing_utilities.GetScalarNodalResults(model_part, KratosMultiphysics.TEMPERATURE)
            expected_output = rom_testing_utilities.GetScalarNodalResults(default_model_part, KratosMultiphysics.TEMPERATURE)
            self._CheckResults(obtained_output, expected_output, 1.0e-12)

    def testCreateDecoder(self):
        # The decoder is built from the ROM folder alone (no analysis stage, no solver): its rows follow the given DOFs,
        # whatever their order, and its type is the one of the RomParameters.json
        import json
        import tempfile
        from KratosMultiphysics.RomApplication.future.rom_decoders import CreateDecoder, GetTypeOfDecoder, LinearDecoder

        class Dof:
            def __init__(self, node_id, variable_name):
                self.node_id, self.variable_name = node_id, variable_name
            def Id(self):
                return self.node_id
            def GetVariable(self):
                return type("Variable", (), {"Name": lambda _: self.variable_name})()

        class DofSet(list):
            def GetValues(self):
                return [10.0 * dof.Id() + (dof.variable_name == "TEMPERATURE") for dof in self]

        node_ids = [3, 7, 11]
        nodal_unknowns = ["DISPLACEMENT_X", "TEMPERATURE"]
        basis = np.random.default_rng(0).standard_normal((len(node_ids) * len(nodal_unknowns), 3))
        dof_set = DofSet([Dof(11, "TEMPERATURE"), Dof(3, "DISPLACEMENT_X"), Dof(7, "TEMPERATURE"), Dof(3, "TEMPERATURE")])
        rows = [5, 0, 3, 1]
        with tempfile.TemporaryDirectory() as rom_data_folder:
            np.save(f"{rom_data_folder}/NodeIds.npy", np.array(node_ids))
            np.save(f"{rom_data_folder}/RightBasisMatrix.npy", basis)
            rom_parameters = {"rom_settings": {"nodal_unknowns": nodal_unknowns, "number_of_rom_dofs": 2}}
            with open(f"{rom_data_folder}/RomParameters.json", "w") as parameter_file:
                json.dump(rom_parameters, parameter_file)
            decoder = CreateDecoder(rom_data_folder, dof_set)
            self.assertIsInstance(decoder, LinearDecoder)
            self.assertTrue(np.array_equal(decoder.Jacobian(None), basis[rows, :2]))
            q = np.array([0.5, -2.0])
            self.assertTrue(np.allclose(decoder.Decode(q), np.array(dof_set.GetValues()) + basis[rows, :2] @ q, rtol=0.0, atol=1.0e-14))

            rom_parameters["type_of_decoder"] = "unknown"
            with open(f"{rom_data_folder}/RomParameters.json", "w") as parameter_file:
                json.dump(rom_parameters, parameter_file)
            with self.assertRaisesRegex(Exception, "Unknown 'type_of_decoder' 'unknown'"):
                CreateDecoder(rom_data_folder, dof_set)

        # The type of the RomParameters.json, else the one of the legacy projection strategy, else linear
        self.assertEqual(GetTypeOfDecoder({}), "linear")
        self.assertEqual(GetTypeOfDecoder({"projection_strategy": "galerkin"}), "linear")
        self.assertEqual(GetTypeOfDecoder({"projection_strategy": "galerkin_ann"}), "ann_enhanced")
        self.assertEqual(GetTypeOfDecoder({"projection_strategy": "lspg_rbf"}), "rbf_enhanced")
        self.assertEqual(GetTypeOfDecoder({"projection_strategy": "galerkin_ann", "type_of_decoder": "linear"}), "linear")


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
