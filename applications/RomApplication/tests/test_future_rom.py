import importlib
import shutil
import tempfile
from pathlib import Path

import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
if kratos_utilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    import KratosMultiphysics.ConvectionDiffusionApplication


@KratosUnittest.skipUnless(hasattr(KratosMultiphysics, "Future"), "Kratos is not compiled with KRATOS_USE_FUTURE")
class TestFutureRom(KratosUnittest.TestCase):
    """Checks the ROM based on the Future schemes against the legacy ROM builder and solvers."""

    def _RunLegacyRom(self, parameters_filename, parameters=None):
        if parameters is None:
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
        simulation = rom_testing_utilities.SetUpSimulationInstance(KratosMultiphysics.Model(), parameters)
        simulation.Run()
        return simulation._GetSolver().GetComputingModelPart()

    def _RunFutureRom(self, parameters_filename, projection_backend, max_iterations=None, run_hrom=False, parameters=None):
        simulation = self._CreateFutureRomSimulation(parameters_filename, projection_backend, max_iterations, run_hrom, parameters)
        simulation.Run()
        return simulation._GetSolver().GetComputingModelPart()

    def _CreateFutureRomSimulation(self, parameters_filename, projection_backend, max_iterations=None, run_hrom=False, parameters=None, train_hrom=False):
        from KratosMultiphysics.RomApplication.future.future_rom_analysis import CreateFutureRomAnalysisInstance

        if parameters is None:
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
        analysis_stage_module_name = parameters["analysis_stage"].GetString()
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_module_name.split('.')[-1].split('_'))
        analysis_stage_class = getattr(importlib.import_module(analysis_stage_module_name), analysis_stage_class_name)

        future_rom_settings = KratosMultiphysics.Parameters("""{}""")
        future_rom_settings.AddString("projection_backend", projection_backend)
        if max_iterations is not None: # Otherwise taken from the solver settings of the case
            future_rom_settings.AddInt("max_iterations", max_iterations)
        future_rom_settings.AddBool("run_hrom", run_hrom)
        future_rom_settings.AddBool("train_hrom", train_hrom)
        return CreateFutureRomAnalysisInstance(analysis_stage_class, KratosMultiphysics.Model(), parameters, future_rom_settings)

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
    def testStructuralStaticHRom2D(self):
        # Same HROM (Galerkin, the weights of the legacy HROM test) with the Future HROM scheme and with the legacy builder and solver.
        # Only the first step is compared: in the second one the load does not change, and the legacy HROM adds a spurious increment
        # of relative size 8e-8 (as the legacy ROM does, see testStructuralStaticRom2D), while the Future one keeps its solution
        with KratosUnittest.WorkFolderScope("structural_static_test_files/HROM/", __file__):
            def GetParameters():
                with open("../ProjectParameters.json",'r') as parameter_file:
                    parameters = KratosMultiphysics.Parameters(parameter_file.read())
                parameters["problem_data"]["end_time"].SetDouble(parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble())
                return parameters
            legacy_model_part = self._RunLegacyRom(None, parameters=GetParameters())
            expected_output = rom_testing_utilities.GetVectorNodalResults(legacy_model_part, KratosMultiphysics.DISPLACEMENT)
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    future_model_part = self._RunFutureRom(None, projection_backend, run_hrom=True, parameters=GetParameters())
                    obtained_output = rom_testing_utilities.GetVectorNodalResults(future_model_part, KratosMultiphysics.DISPLACEMENT)
                    self._CheckResults(obtained_output, expected_output, 1.0e-9)

    def _GetGalerkinHRomParameters(self, input_filename=None):
        """The problem of the RomManager test, loaded with 0.1 and without output, to be run in future_rom_test_files/GalerkinHROM"""
        with open("../../test_rom_manager/ProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        solver_settings = parameters["solver_settings"]
        for settings, key in [(solver_settings["model_import_settings"], "input_filename"), (solver_settings["material_import_settings"], "materials_filename")]:
            settings[key].SetString("../../test_rom_manager/" + settings[key].GetString())
        if input_filename is not None:
            solver_settings["model_import_settings"]["input_filename"].SetString(input_filename)
        for load_process in parameters["processes"]["loads_process_list"].values():
            load_process["Parameters"]["value"].SetString("(0.1)")
        parameters.RemoveValue("output_processes")
        return parameters

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticHRomLoadOutsideSelectedElements(self):
        # Galerkin HROM whose selected condition (a load) shares no node with its selected element: the rows of the DOFs not
        # touched by the selected element are empty and hold the load, so they must not be modified. The HROM integrates the
        # residual of the ROM exactly, so both must give the same solution
        with KratosUnittest.WorkFolderScope("future_rom_test_files/GalerkinHROM/", __file__):
            rom_model_part = self._RunFutureRom(None, "cpp", parameters=self._GetGalerkinHRomParameters())
            expected_output = rom_testing_utilities.GetVectorNodalResults(rom_model_part, KratosMultiphysics.DISPLACEMENT)
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    hrom_model_part = self._RunFutureRom(None, projection_backend, run_hrom=True, parameters=self._GetGalerkinHRomParameters())
                    obtained_output = rom_testing_utilities.GetVectorNodalResults(hrom_model_part, KratosMultiphysics.DISPLACEMENT)
                    self._CheckResults(obtained_output, expected_output, 1.0e-9)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticHRomModelPart(self):
        # The HROM run on the HROM model part (written as the HROM training does) must give the solution of the HROM run on the
        # full one, on the nodes of the HROM model part
        with KratosUnittest.WorkFolderScope("future_rom_test_files/GalerkinHROM/", __file__):
            full_model_part = self._RunFutureRom(None, "cpp", run_hrom=True, parameters=self._GetGalerkinHRomParameters())
            for projection_backend in ["numpy", "cpp"]:
                with self.subTest(projection_backend=projection_backend):
                    hrom_model_part = self._RunFutureRom(None, projection_backend, run_hrom=True,
                        parameters=self._GetGalerkinHRomParameters("simple_structure_multiple_parametersHROM"))
                    self.assertLess(hrom_model_part.NumberOfElements(), full_model_part.NumberOfElements())
                    expected_output = [full_model_part.GetNode(node.Id).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) for node in hrom_model_part.Nodes]
                    obtained_output = [node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) for node in hrom_model_part.Nodes]
                    self._CheckResults(obtained_output, expected_output, 1.0e-12)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticTotalLagrangianHRomAllWeightsOne(self):
        # HROM with all the elements and conditions and all the weights equal to 1: it assembles the system of the ROM, so it must
        # give the solution of the ROM (geometrically nonlinear cantilever). Run in a copy of the case, as the HROM weights are written
        with tempfile.TemporaryDirectory() as case_folder:
            shutil.copytree(Path(__file__).parent / "future_rom_test_files/TotalLagrangianROM", case_folder, dirs_exist_ok=True)
            with KratosUnittest.WorkFolderScope(case_folder, __file__):
                for projection_backend in ["numpy", "cpp"]:
                    with self.subTest(projection_backend=projection_backend):
                        rom_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, max_iterations=20)
                        expected_output = rom_testing_utilities.GetVectorNodalResults(rom_model_part, KratosMultiphysics.DISPLACEMENT)

                        root_model_part = rom_model_part.GetRootModelPart()
                        for name, entities in [("Element", root_model_part.Elements), ("Condition", root_model_part.Conditions)]:
                            entity_ids = [entity.Id for entity in entities]
                            np.save(f"rom_data/HROM_{name}Ids.npy", np.array(entity_ids))
                            np.save(f"rom_data/HROM_{name}Weights.npy", np.ones(len(entity_ids)))
                        hrom_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, max_iterations=20, run_hrom=True)
                        obtained_output = rom_testing_utilities.GetVectorNodalResults(hrom_model_part, KratosMultiphysics.DISPLACEMENT)
                        self._CheckResults(obtained_output, expected_output, 1.0e-12)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticTotalLagrangianHRomTraining(self):
        # HROM training of the Future ROM on the geometrically nonlinear cantilever: the residual snapshots must be the ones of the
        # legacy training (ProjectedResidualsOutputProcess) and the trained HROM must give the solution of the ROM. Run in a copy of
        # the case, as the training writes the HROM weights and the HROM model part
        with tempfile.TemporaryDirectory() as case_folder:
            shutil.copytree(Path(__file__).parent / "future_rom_test_files/TotalLagrangianROM", case_folder, dirs_exist_ok=True)
            with KratosUnittest.WorkFolderScope(case_folder, __file__):
                with open("ProjectParameters.json",'r') as parameter_file:
                    parameters = KratosMultiphysics.Parameters(parameter_file.read())
                parameters["output_processes"].AddEmptyArray("rom_output")
                parameters["output_processes"]["rom_output"].Append(KratosMultiphysics.Parameters("""{
                    "python_module" : "projected_residuals_output_process",
                    "kratos_module" : "KratosMultiphysics.RomApplication",
                    "process_name"  : "ProjectedResidualsOutputProcess",
                    "Parameters"    : {
                        "model_part_name" : "Structure",
                        "output_path"     : "rom_data/Residuals"
                    }
                }"""))
                self._RunLegacyRom(None, parameters=parameters)
                end_time = parameters["problem_data"]["end_time"].GetDouble()
                time_step = parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
                legacy_snapshots = [np.load(f"rom_data/Residuals/Residual_{step}.npy") for step in range(1, round(end_time / time_step) + 1)]

                for projection_backend in ["numpy", "cpp"]:
                    with self.subTest(projection_backend=projection_backend):
                        simulation = self._CreateFutureRomSimulation("ProjectParameters.json", projection_backend, max_iterations=20, train_hrom=True)
                        simulation.Run()
                        self.assertEqual(len(simulation.hrom_residual_snapshots), len(legacy_snapshots))
                        for future_snapshot, legacy_snapshot in zip(simulation.hrom_residual_snapshots, legacy_snapshots):
                            self._CheckResults(future_snapshot, legacy_snapshot, 1.0e-10)
                        rom_model_part = simulation._GetSolver().GetComputingModelPart()
                        expected_output = rom_testing_utilities.GetVectorNodalResults(rom_model_part, KratosMultiphysics.DISPLACEMENT)

                        self.assertLess(len(np.load("rom_data/HROM_ElementIds.npy")), rom_model_part.NumberOfElements())
                        self.assertTrue(Path("cantileverHROM.mdpa").exists())
                        hrom_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, max_iterations=20, run_hrom=True)
                        obtained_output = rom_testing_utilities.GetVectorNodalResults(hrom_model_part, KratosMultiphysics.DISPLACEMENT)
                        self._CheckResults(obtained_output, expected_output, 1.0e-6)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticTotalLagrangianHRomTrainingFromSnapshots(self):
        # HROM training from full order snapshots (rows of the basis, one column per step), without solving the ROM:
        #  - the snapshots of the ROM itself: the same reduced coordinates and residual snapshots as the training of the solved ROM
        #  - the FOM snapshots: their projection on the basis, and the trained HROM gives the ROM solution
        from unittest import mock
        from KratosMultiphysics.RomApplication.future.future_rom_analysis import TrainHRomFromSnapshots
        from KratosMultiphysics.RomApplication.future.rom_future_solver import FutureRomSolver

        with tempfile.TemporaryDirectory() as case_folder:
            shutil.copytree(Path(__file__).parent / "future_rom_test_files/TotalLagrangianROM", case_folder, dirs_exist_ok=True)
            with KratosUnittest.WorkFolderScope(case_folder, __file__):
                with open("ProjectParameters.json", "r") as parameter_file:
                    parameters = KratosMultiphysics.Parameters(parameter_file.read())
                analysis_stage_module_name = parameters["analysis_stage"].GetString()
                analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_module_name.split('.')[-1].split('_'))
                analysis_stage_class = getattr(importlib.import_module(analysis_stage_module_name), analysis_stage_class_name)
                node_ids = np.load("rom_data/NodeIds.npy")

                class FomAnalysis(analysis_stage_class):
                    """Saves the snapshots in the rows of the basis"""
                    def FinalizeSolutionStep(self):
                        super().FinalizeSolutionStep()
                        model_part = self._GetSolver().GetComputingModelPart()
                        variables = (KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y)
                        self.snapshots.append([model_part.GetNode(int(i)).GetSolutionStepValue(variable) for i in node_ids for variable in variables])
                fom = FomAnalysis(KratosMultiphysics.Model(), parameters.Clone())
                fom.snapshots = []
                fom.Run()
                np.save("fom_snapshots.npy", np.array(fom.snapshots).T)

                for projection_backend in ["numpy", "cpp"]:
                    with self.subTest(projection_backend=projection_backend):
                        future_rom_settings = KratosMultiphysics.Parameters("""{}""")
                        future_rom_settings.AddString("projection_backend", projection_backend)
                        future_rom_settings.AddInt("max_iterations", 20)
                        rom = self._CreateFutureRomSimulation("ProjectParameters.json", projection_backend, max_iterations=20, train_hrom=True)
                        rom.Run()
                        expected_output = rom_testing_utilities.GetVectorNodalResults(rom._GetSolver().GetComputingModelPart(), KratosMultiphysics.DISPLACEMENT)
                        # The ROM states in the rows of the basis (the initial state is zero)
                        basis = np.load("rom_data/RightBasisMatrix.npy")[:, :rom.future_rom_solver.decoder.NumberOfRomDofs()]
                        rom_snapshots = basis @ np.array(rom.rom_solutions).T

                        for snapshots in [rom_snapshots, "fom_snapshots.npy"]:
                            with mock.patch.object(FutureRomSolver, "SolveSolutionStep", side_effect=Exception("The ROM must not be solved.")):
                                training = TrainHRomFromSnapshots(analysis_stage_class, KratosMultiphysics.Model(), parameters.Clone(), snapshots, future_rom_settings)
                            if isinstance(snapshots, str):
                                for obtained, fom_snapshot in zip(training.rom_solutions, fom.snapshots):
                                    self._CheckResults(obtained, np.linalg.lstsq(basis, fom_snapshot, rcond=None)[0], 1.0e-10)
                            else:
                                for obtained, expected in zip(training.rom_solutions, rom.rom_solutions):
                                    self._CheckResults(obtained, expected, 1.0e-12)
                                for obtained, expected in zip(training.hrom_residual_snapshots, rom.hrom_residual_snapshots):
                                    self._CheckResults(obtained, expected, 1.0e-10)

                            hrom_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, max_iterations=20, run_hrom=True)
                            obtained_output = rom_testing_utilities.GetVectorNodalResults(hrom_model_part, KratosMultiphysics.DISPLACEMENT)
                            self._CheckResults(obtained_output, expected_output, 1.0e-6)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticTotalLagrangianHRomNoConditions(self):
        # HROM with all the elements and no condition selected, with the empty weights saved as 1D and as 2D (one column per set) arrays:
        # the load of the cantilever is applied by its conditions only, so the HROM gives a null displacement
        with tempfile.TemporaryDirectory() as case_folder:
            shutil.copytree(Path(__file__).parent / "future_rom_test_files/TotalLagrangianROM", case_folder, dirs_exist_ok=True)
            with KratosUnittest.WorkFolderScope(case_folder, __file__):
                model_part = KratosMultiphysics.Model().CreateModelPart("cantilever")
                KratosMultiphysics.ModelPartIO("cantilever", KratosMultiphysics.IO.READ | KratosMultiphysics.IO.SKIP_TIMER).ReadModelPart(model_part)
                element_ids = np.array([element.Id for element in model_part.Elements])
                np.save("rom_data/HROM_ElementIds.npy", element_ids)
                np.save("rom_data/HROM_ConditionIds.npy", np.array([], dtype=int))
                for empty_weights in [np.ones(0), np.ones((0, 1))]:
                    with self.subTest(empty_weights_shape=empty_weights.shape):
                        np.save("rom_data/HROM_ElementWeights.npy", np.ones((len(element_ids),) + empty_weights.shape[1:]))
                        np.save("rom_data/HROM_ConditionWeights.npy", empty_weights)
                        hrom_model_part = self._RunFutureRom("ProjectParameters.json", "cpp", max_iterations=20, run_hrom=True)
                        obtained_output = rom_testing_utilities.GetVectorNodalResults(hrom_model_part, KratosMultiphysics.DISPLACEMENT)
                        self.assertLess(np.abs(obtained_output).max(), 1.0e-12)

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

    def testDecoderEncode(self):
        # Encode inverts Decode: the linear decoder by least squares on the basis, any other one by Gauss-Newton (base class, checked
        # with a quadratic decoder). For a vector out of the decoder range the linear one gives the least squares fit
        from KratosMultiphysics.RomApplication.future.rom_decoders import LinearDecoder, RomDecoder

        class QuadraticDecoder(RomDecoder):
            """u = u0 + Phi q + Psi (q * q)"""
            def __init__(self, phi, psi, u0):
                self.phi, self.psi, self.u0 = phi, psi, u0
            def NumberOfRomDofs(self):
                return self.phi.shape[1]
            def Decode(self, q):
                return self.u0 + self.phi @ q + self.psi @ (q * q)
            def Jacobian(self, q):
                return self.phi + 2.0 * self.psi * q

        rng = np.random.default_rng(0)
        phi, psi, u0 = rng.standard_normal((20, 3)), 0.1 * rng.standard_normal((20, 3)), rng.standard_normal(20)
        q = np.array([0.5, -1.0, 2.0])

        linear = LinearDecoder(phi, u0)
        self.assertTrue(np.allclose(linear.Encode(linear.Decode(q)), q, rtol=0.0, atol=1.0e-12))
        u = rng.standard_normal(20)
        self.assertLess(np.abs(phi.T @ (linear.Decode(linear.Encode(u)) - u)).max(), 1.0e-12)

        quadratic = QuadraticDecoder(phi, psi, u0)
        self.assertTrue(np.allclose(quadratic.Encode(quadratic.Decode(q), q + 0.1), q, rtol=0.0, atol=1.0e-12))

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
