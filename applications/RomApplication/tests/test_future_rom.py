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
