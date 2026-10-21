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
            def TorchLinearDecoderFactory(rows, u0):
                from KratosMultiphysics.RomApplication.future.rom_decoders import TorchDecoder
                phi = torch.as_tensor(np.load("rom_data/RightBasisMatrix.npy")[rows, :])
                u0 = torch.as_tensor(u0)
                return TorchDecoder(lambda q: u0 + phi @ q, phi.shape[1])
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    linear_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend)
                    expected_output = rom_testing_utilities.GetVectorNodalResults(linear_model_part, KratosMultiphysics.DISPLACEMENT)
                    torch_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, decoder_factory=TorchLinearDecoderFactory)
                    obtained_output = rom_testing_utilities.GetVectorNodalResults(torch_model_part, KratosMultiphysics.DISPLACEMENT)
                    self._CheckResults(obtained_output, expected_output, 1.0e-12)

    def testQuadraticDecoder(self):
        from KratosMultiphysics.RomApplication.future.rom_decoders import QuadraticDecoder, QuadraticFeatures, FitQuadraticManifold
        random_generator = np.random.default_rng(0)
        phi = np.linalg.qr(random_generator.standard_normal((20, 3)))[0]
        h = random_generator.standard_normal((20, 6))
        h -= phi @ (phi.T @ h) # out of span(Phi), as the fitted ones
        decoder = QuadraticDecoder(phi, h, random_generator.standard_normal(20))

        # Jacobian against central finite differences
        q, delta = random_generator.standard_normal(3), 1.0e-6
        finite_differences = np.array([(decoder.Decode(q + delta * e) - decoder.Decode(q - delta * e)) / (2.0 * delta) for e in np.eye(3)]).T
        self.assertLess(np.abs(decoder.Jacobian(q) - finite_differences).max(), 1.0e-8)

        # The fit recovers H from snapshots on the quadratic manifold
        reduced_coordinates = random_generator.standard_normal((3, 30))
        snapshots = phi @ reduced_coordinates + h @ QuadraticFeatures(reduced_coordinates)
        self.assertLess(np.abs(FitQuadraticManifold(snapshots, phi, regularization=0.0) - h).max(), 1.0e-10)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticQuadraticDecoder2D(self):
        # Beam clamped at both ends with two independent loads: the membrane stretching makes the axial displacement quadratic in
        # the deflections. With 2 modes the quadratic decoder is about 18 times more accurate than the linear one
        from KratosMultiphysics.RomApplication.future.rom_decoders import QuadraticDecoder
        with KratosUnittest.WorkFolderScope("future_rom_test_files/QuadraticClampedBeamROM/", __file__):
            def QuadraticDecoderFactory(rows, u0):
                return QuadraticDecoder(np.load("rom_data/RightBasisMatrix.npy")[rows, :], np.load("rom_data/QuadraticMatrix.npy")[rows, :], u0)
            def RelativeError(model_part):
                obtained_output = np.array([node.GetSolutionStepValue(variable) for node in model_part.Nodes for variable in [KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y]])
                expected_output = np.load("ExpectedOutputFOM.npy")
                return np.linalg.norm(obtained_output - expected_output) / np.linalg.norm(expected_output)
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    linear_error = RelativeError(self._RunFutureRom("ProjectParameters.json", projection_backend))
                    quadratic_error = RelativeError(self._RunFutureRom("ProjectParameters.json", projection_backend, decoder_factory=QuadraticDecoderFactory))
                    self.assertLess(quadratic_error, 3.0e-2)
                    self.assertLess(quadratic_error, 0.1 * linear_error)

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


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
