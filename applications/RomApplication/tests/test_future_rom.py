import importlib
import json
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

    def _RunLegacyRom(self, parameters_filename):
        with open(parameters_filename,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        simulation = rom_testing_utilities.SetUpSimulationInstance(KratosMultiphysics.Model(), parameters)
        simulation.Run()
        return simulation._GetSolver().GetComputingModelPart()

    def _RunFutureRom(self, parameters_filename, projection_backend, max_iterations=None):
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
        simulation = CreateFutureRomAnalysisInstance(analysis_stage_class, KratosMultiphysics.Model(), parameters, future_rom_settings)
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


    def testDecoderEncode(self):
        # Encode inverts Decode on the manifold of the decoder: exactly for the linear decoder (any full rank basis, with origin),
        # by Gauss-Newton for the others (base class)
        from KratosMultiphysics.RomApplication.future.rom_decoders import LinearDecoder, RomDecoder
        rng = np.random.default_rng(0)
        phi, u0, q = rng.standard_normal((40, 3)), rng.standard_normal(40), rng.standard_normal(3)
        decoder = LinearDecoder(phi, u0)
        self.assertLess(np.linalg.norm(decoder.Encode(decoder.Decode(q)) - q), 1.0e-12 * np.linalg.norm(q))

        class QuadraticTermDecoder(RomDecoder):
            def NumberOfRomDofs(self):
                return 3
            def Decode(self, q):
                return phi @ q + 0.1 * u0 * (q @ q)
            def Jacobian(self, q):
                return phi + 0.2 * np.outer(u0, q)
        decoder = QuadraticTermDecoder()
        self.assertLess(np.linalg.norm(decoder.Encode(decoder.Decode(q)) - q), 1.0e-10 * np.linalg.norm(q))

    def _RunTotalLagrangian(self, rom_folder=None, time_step=None):
        # The cantilever of testStructuralStaticTotalLagrangianRom2D: the FOM (rom_folder None) or the Future ROM of the given folder.
        # Returns the displacements of each step (one column per step, rows node by node) and the model of each step
        from KratosMultiphysics.RomApplication.future.future_rom_analysis import CreateFutureRomAnalysisInstance
        from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

        def Parameters():
            with open("ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            if time_step is not None:
                parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(time_step)
            return parameters

        def Recording(base_class):
            class RecordingAnalysis(base_class):
                def Initialize(self):
                    super().Initialize()
                    self.snapshots, self.models = [], []
                def FinalizeSolutionStep(self):
                    super().FinalizeSolutionStep()
                    self.snapshots.append(np.array(rom_testing_utilities.GetVectorNodalResults(self._GetSolver().GetComputingModelPart(), KratosMultiphysics.DISPLACEMENT)).ravel())
                    decoder = self.future_rom_solver.decoder if hasattr(self, "future_rom_solver") else None
                    self.models.append(decoder.GetActiveModel() if hasattr(decoder, "GetActiveModel") else 0)
            return RecordingAnalysis

        if rom_folder is None:
            simulation = Recording(StructuralMechanicsAnalysis)(KratosMultiphysics.Model(), Parameters())
        else:
            future_rom_settings = KratosMultiphysics.Parameters('{"max_iterations" : 20}')
            future_rom_class = type(CreateFutureRomAnalysisInstance(StructuralMechanicsAnalysis, KratosMultiphysics.Model(), Parameters(), future_rom_settings, rom_basis_output_folder=str(rom_folder)))
            simulation = Recording(future_rom_class)(KratosMultiphysics.Model(), Parameters())
        simulation.Run()
        return np.array(simulation.snapshots).T, simulation.models

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticTotalLagrangianOneLocalRom(self):
        # A LocalRomDecoder with one model gives the solution of that model alone: the local ROMs reduce to the single decoder
        with KratosUnittest.WorkFolderScope("future_rom_test_files/TotalLagrangianROM/", __file__), tempfile.TemporaryDirectory() as temp_dir:
            single, _ = self._RunTotalLagrangian(Path("rom_data"))
            local_folder = Path(temp_dir)
            shutil.copytree("rom_data", local_folder / "local_0")
            shutil.copy(Path("rom_data") / "NodeIds.npy", local_folder)
            with open(Path("rom_data") / "RomParameters.json", 'r') as parameter_file:
                rom_parameters = json.load(parameter_file)
            rom_parameters["local_roms"] = {"models": ["local_0"], "selector": {"type": "process_info", "variable": "TIME", "thresholds": []}}
            with open(local_folder / "RomParameters.json", 'w') as parameter_file:
                json.dump(rom_parameters, parameter_file)
            local, models = self._RunTotalLagrangian(local_folder)
            self.assertEqual(models, [0] * len(models))
            self._CheckResults(local, single, 1.0e-14)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticTotalLagrangianLocalRoms(self):
        # Two local affine models (origin = mean of their snapshots, 2 modes each), fitted to the first and second half of the load
        # steps (2 overlapping steps), switched by TIME or by the nearest centre, against the global basis of the same size (2 modes)
        with KratosUnittest.WorkFolderScope("future_rom_test_files/TotalLagrangianROM/", __file__), tempfile.TemporaryDirectory() as temp_dir:
            time_step, number_of_steps, number_of_modes = 0.1, 30, 2
            fom, _ = self._RunTotalLagrangian(time_step=time_step)
            self.assertEqual(fom.shape[1], number_of_steps)
            folder = Path(temp_dir)
            node_ids = np.array([node.Id for node in self._ModelPartNodes()])

            def WriteRomParameters(path, local_roms=None):
                rom_parameters = {"rom_settings": {"nodal_unknowns": ["DISPLACEMENT_X", "DISPLACEMENT_Y"], "number_of_rom_dofs": number_of_modes}}
                if local_roms is not None:
                    rom_parameters["local_roms"] = local_roms
                path.mkdir(parents=True, exist_ok=True)
                np.save(path / "NodeIds.npy", node_ids)
                with open(path / "RomParameters.json", 'w') as parameter_file:
                    json.dump(rom_parameters, parameter_file)

            WriteRomParameters(folder / "global")
            np.save(folder / "global" / "RightBasisMatrix.npy", np.linalg.svd(fom, full_matrices=False)[0][:, :number_of_modes])
            global_rom, _ = self._RunTotalLagrangian(folder / "global", time_step)

            half = number_of_steps // 2
            centres = []
            for k, steps in enumerate([np.arange(0, half + 2), np.arange(half - 2, number_of_steps)]):
                local_folder = folder / "local" / f"local_{k}"
                local_folder.mkdir(parents=True)
                centre = fom[:, steps].mean(axis=1)
                np.save(local_folder / "ReferenceSnapshot.npy", centre)
                np.save(local_folder / "RightBasisMatrix.npy", np.linalg.svd(fom[:, steps] - centre[:, None], full_matrices=False)[0][:, :number_of_modes])
                centres.append(centre)
            np.save(folder / "local" / "Centres.npy", np.array(centres).T)

            error = lambda solution: np.linalg.norm(solution - fom) / np.linalg.norm(fom)
            selectors = {
                "process_info": {"type": "process_info", "variable": "TIME", "thresholds": [(half + 0.5) * time_step]},
                "nearest_centre": {"type": "nearest_centre", "centres": "Centres.npy", "margin": 0.0}}
            for name, selector in selectors.items():
                with self.subTest(selector=name):
                    WriteRomParameters(folder / "local", {"models": ["local_0", "local_1"], "selector": selector})
                    local_rom, models = self._RunTotalLagrangian(folder / "local", time_step)
                    # One switch, from the first model to the second, at the half of the load (the nearest centre within 2 steps)
                    switch_steps = list(np.flatnonzero(np.diff(models)) + 1)
                    self.assertEqual(len(switch_steps), 1)
                    self.assertEqual(models[0], 0)
                    self.assertEqual(models[-1], 1)
                    self.assertLessEqual(abs(switch_steps[0] - half), 2)
                    self.assertLess(error(local_rom), 0.25 * error(global_rom))
                    self.assertLess(error(local_rom), 1.0e-4)

    def _ModelPartNodes(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Structure")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        KratosMultiphysics.ModelPartIO("cantilever").ReadModelPart(model_part)
        return list(model_part.Nodes)


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
