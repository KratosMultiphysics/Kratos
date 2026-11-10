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

    def _CreateFutureRomSimulation(self, parameters_filename, projection_backend, max_iterations=None, run_hrom=False, parameters=None, train_hrom=False, extra_settings=None):
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
        for name, value in (extra_settings or {}).items():
            future_rom_settings.AddBool(name, value)
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
    def testStructuralStaticTotalLagrangianHRomTrainingReplay(self):
        # HROM training by replay: the ROM run saves its reduced coordinates; the replay imposes them in each step without solving and
        # gives the residual snapshots of the run, and its trained HROM gives the solution of the ROM
        with tempfile.TemporaryDirectory() as case_folder:
            shutil.copytree(Path(__file__).parent / "future_rom_test_files/TotalLagrangianROM", case_folder, dirs_exist_ok=True)
            with KratosUnittest.WorkFolderScope(case_folder, __file__):
                for projection_backend in ["numpy", "cpp"]:
                    with self.subTest(projection_backend=projection_backend):
                        run = self._CreateFutureRomSimulation("ProjectParameters.json", projection_backend, max_iterations=20, train_hrom=True, extra_settings={"save_rom_solutions": True})
                        run.Run()
                        rom_model_part = run._GetSolver().GetComputingModelPart()
                        expected_output = rom_testing_utilities.GetVectorNodalResults(rom_model_part, KratosMultiphysics.DISPLACEMENT)
                        self.assertEqual(np.load("rom_data/RomSolutions.npy").shape, (len(run.rom_solutions), len(run.rom_solutions[0])))

                        replay = self._CreateFutureRomSimulation("ProjectParameters.json", projection_backend, train_hrom=True, extra_settings={"replay_rom_solutions": True})
                        replay.Initialize()
                        def NoSolve():
                            raise Exception("The replay must not solve.")
                        replay.future_rom_solver.SolveSolutionStep = NoSolve
                        replay.RunSolutionLoop()
                        replay.Finalize()
                        self.assertEqual(len(replay.hrom_residual_snapshots), len(run.hrom_residual_snapshots))
                        for replayed, solved in zip(replay.hrom_residual_snapshots, run.hrom_residual_snapshots):
                            self._CheckResults(replayed, solved, 1.0e-12)
                        replayed_output = rom_testing_utilities.GetVectorNodalResults(replay._GetSolver().GetComputingModelPart(), KratosMultiphysics.DISPLACEMENT)
                        self._CheckResults(replayed_output, expected_output, 1.0e-12)

                        hrom_model_part = self._RunFutureRom("ProjectParameters.json", projection_backend, max_iterations=20, run_hrom=True)
                        obtained_output = rom_testing_utilities.GetVectorNodalResults(hrom_model_part, KratosMultiphysics.DISPLACEMENT)
                        self._CheckResults(obtained_output, expected_output, 1.0e-6)

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
