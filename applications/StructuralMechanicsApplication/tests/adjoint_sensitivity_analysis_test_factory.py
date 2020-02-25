from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics import IsDistributedRun
from structural_mechanics_test_factory import SelectAndVerifyLinearSolver

has_hdf5_application = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")

class AdjointSensitivityAnalysisTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            # Reading the ProjectParameters
            with open(self.primal_file_name,'r') as parameter_file:
                primal_parameters = KratosMultiphysics.Parameters(parameter_file.read())
            with open(self.adjoint_file_name,'r') as parameter_file:
                self.adjoint_parameters = KratosMultiphysics.Parameters(parameter_file.read())
            self.problem_name = primal_parameters["problem_data"]["problem_name"].GetString()
            self.model_part_name = primal_parameters["solver_settings"]["model_part_name"].GetString()

            # To avoid many prints
            if (primal_parameters["problem_data"]["echo_level"].GetInt() == 0 or self.adjoint_parameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            SelectAndVerifyLinearSolver(primal_parameters, self.skipTest)
            SelectAndVerifyLinearSolver(self.adjoint_parameters, self.skipTest)

            # solve primal problem
            model_primal = KratosMultiphysics.Model()
            primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, primal_parameters)
            primal_analysis.Run()
            # create adjoint analysis
            model_adjoint = KratosMultiphysics.Model()
            self.adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, self.adjoint_parameters)
            self.adjoint_analysis.Initialize()

    def test_execution(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.adjoint_analysis.RunSolutionLoop()
            self.perform_additional_checks()

    def perform_additional_checks(self):
        pass

    def tearDown(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.adjoint_analysis.Finalize()
            kratos_utilities.DeleteFileIfExisting(self.problem_name + ".time")
            kratos_utilities.DeleteFileIfExisting(self.model_part_name + ".h5")
            kratos_utilities.DeleteFileIfExisting(self.model_part_name + "-1.0000.h5")
            kratos_utilities.DeleteFileIfExisting(self.model_part_name + "-1.1000.h5")


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisBeamStructureLocalStress(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_local_stress_adjoint_parameters.json"

    def perform_additional_checks(self):
        reference_values = [-87.62277093392399, 9.497391494932984, 38.125186783868, 0.6250049974719261, 0.15624887499699122]
        sensitivities_to_check = []
        element_list = [1,2,3,4,5,6,10]
        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_model_part.Elements[element_id].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))
        sensitivities_to_check.append(adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 3)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[4], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[5], reference_values[2], 3)
        self.assertAlmostEqual(sensitivities_to_check[6], reference_values[3], 3)
        self.assertAlmostEqual(sensitivities_to_check[7], reference_values[4], 5)


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisBeamStructureNodalDisplacement(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_nodal_disp_adjoint_parameters.json"

    def perform_additional_checks(self):
        reference_values = [-0.45410279537614157, -0.37821875982596204, -0.006200296058668847, 0.0004340210813670321]
        sensitivities_to_check = []
        element_list = [1,6,10]
        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_model_part.Elements[element_id].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))
        sensitivities_to_check.append(adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 4)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 4)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 4)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[3], 5)


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisBeamStructureStrainEnergy(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_strain_energy_adjoint_parameters.json"

    def perform_additional_checks(self):
        reference_values = [-9.082055907522943, -7.5643751965193164, -0.12400592117339182, 0.017360843254681547]
        sensitivities_to_check = []
        element_list = [1,6,10]
        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_model_part.Elements[element_id].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))
        sensitivities_to_check.append(adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 4)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 4)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 4)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[3], 5)


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisBeamStructureNodalReaction(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_nodal_reaction_adjoint_parameters.json"

    def perform_additional_checks(self):
        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        reference_values = [-0.31249774999397384, -1.249]
        sensitivities_to_check = []
        sensitivities_to_check.append(adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])
        sensitivities_to_check.append(adjoint_model_part.Elements[10].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 4)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 2)


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisShellStructureLocalStress(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_local_stress_adjoint_parameters.json"

    def perform_additional_checks(self):
        reference_values = [1.7135092490964121, -6.860092387341681, 0.14749301178647778, -0.0823339298948347]
        sensitivities_to_check = []
        element_list = [1,2,8]
        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_model_part.Elements[element_id].GetValue(StructuralMechanicsApplication.THICKNESS_SENSITIVITY))
        sensitivities_to_check.append(adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[3], 5)


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisShellStructureNodalDisplacement(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_nodal_disp_adjoint_parameters.json"

    def perform_additional_checks(self):
        reference_values = [-0.09916013365433643, -0.23348175177098657, -0.04942512089147077, 0.012125502238309537]
        sensitivities_to_check = []
        element_list = [1,2,8]
        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_model_part.Elements[element_id].GetValue(StructuralMechanicsApplication.THICKNESS_SENSITIVITY))
        sensitivities_to_check.append(adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[3], 5)


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisShellStructureStrainEnergy(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_strain_energy_adjoint_parameters.json"

    def perform_additional_checks(self):
        reference_values = [-0.4958006682716821, -1.1674087588549331, -0.2471256044520311, 0.12125502238309535]
        sensitivities_to_check = []
        element_list = [1,2,8]
        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_model_part.Elements[element_id].GetValue(StructuralMechanicsApplication.THICKNESS_SENSITIVITY))
        sensitivities_to_check.append(adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[3], 5)


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisSpringDamperElement(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_spring_damper_element_3d2n/ProjectParameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_spring_damper_element_3d2n/AdjointParameters.json"


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisLinearTrussStructure(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/linear_truss_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/linear_truss_test_local_stress_adjoint_parameters.json"


@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisNonLinearTrussStructure(AdjointSensitivityAnalysisTestFactory):
    primal_file_name = "adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/nonlinear_truss_test_parameters.json"
    adjoint_file_name = "adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/nonlinear_truss_test_local_stress_adjoint_parameters.json"

if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureLocalStress]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureNodalDisplacement]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureStrainEnergy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureNodalReaction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisShellStructureLocalStress]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisShellStructureNodalDisplacement]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisShellStructureStrainEnergy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisSpringDamperElement]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisLinearTrussStructure]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisNonLinearTrussStructure]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
