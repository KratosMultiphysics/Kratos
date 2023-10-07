import os

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics import IsDistributedRun
from structural_mechanics_test_factory import SelectAndVerifyLinearSolver
import numpy as np

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
        reference_value = -0.31249774999397384
        sensitivity_to_check = adjoint_model_part.Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2]
        self.assertAlmostEqual(sensitivity_to_check, reference_value)


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

def testAdjointStrain():
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("Structure")
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
    model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
    sub_model_part = model_part.CreateSubModelPart("Parts_AREAS")

    n1: KratosMultiphysics.Node = model_part.CreateNewNode(1, 0, 0, 0)
    n2: KratosMultiphysics.Node = model_part.CreateNewNode(2, 1, 0, 0)
    n3: KratosMultiphysics.Node = model_part.CreateNewNode(3, 1, 1, 0)

    n1.AddDof(KratosMultiphysics.DISPLACEMENT_X)
    n1.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
    n1.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
    n2.AddDof(KratosMultiphysics.DISPLACEMENT_X)
    n2.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
    n2.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
    n3.AddDof(KratosMultiphysics.DISPLACEMENT_X)
    n3.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
    n3.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

    n1.AddDof(KratosMultiphysics.ROTATION_X)
    n1.AddDof(KratosMultiphysics.ROTATION_Y)
    n1.AddDof(KratosMultiphysics.ROTATION_Z)
    n2.AddDof(KratosMultiphysics.ROTATION_X)
    n2.AddDof(KratosMultiphysics.ROTATION_Y)
    n2.AddDof(KratosMultiphysics.ROTATION_Z)
    n3.AddDof(KratosMultiphysics.ROTATION_X)
    n3.AddDof(KratosMultiphysics.ROTATION_Y)
    n3.AddDof(KratosMultiphysics.ROTATION_Z)

    n1.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)

    n1.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z)

    p1 = model_part.CreateNewProperties(1)
    material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
    material_settings["Parameters"]["materials_filename"].SetString("/software/kratos/digital_twin/applications/StructuralMechanicsApplication/tests/adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/StructuralMaterials.json")
    KratosMultiphysics.ReadMaterialsUtility(material_settings, model)

    elem1: KratosMultiphysics.Element = model_part.CreateNewElement("AdjointFiniteDifferencingShellThinElement3D3N", 1, [1, 2, 3], p1)
    sub_model_part.AddElements([1])

    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, node.Id)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, node.Id + 1)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, node.Id + 2)
        node.SetSolutionStepValue(KratosMultiphysics.ROTATION_X, node.Id)
        node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Y, node.Id + 1)
        node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Z, node.Id + 2)

    elem1.Initialize(model_part.ProcessInfo)
    residuals = KratosMultiphysics.Matrix()
    lhs = KratosMultiphysics.Matrix()
    rhs = KratosMultiphysics.Vector()
    elem1.CalculateLocalSystem(lhs, rhs, model_part.ProcessInfo)

    ref_gp_values = elem1.CalculateOnIntegrationPoints(StructuralMechanicsApplication.SHELL_STRAIN, model_part.ProcessInfo)
    adjoint_response = StructuralMechanicsApplication.AdjointElementStrainResponseFunction(model_part, KratosMultiphysics.Parameters("""{"weight": 2, "step_size":1e-9, "gradient_mode": "semi_analytic", "strain_type": "y", "element_id": 1}"""))
    adjoint_response.Initialize()

    ad_sensitivities = KratosMultiphysics.Vector()
    adjoint_response.CalculateGradient(elem1, lhs, ad_sensitivities, model_part.ProcessInfo)
    # print(ad_sensitivities)

    ref_value = ref_gp_values[0][1, 1] + ref_gp_values[1][1, 1] + ref_gp_values[2][1, 1]
    print(adjoint_response.CalculateValue(model_part), ref_value / 6)
    ref_value = adjoint_response.CalculateValue(model_part)

    dof_order = [
        KratosMultiphysics.DISPLACEMENT_X,
        KratosMultiphysics.DISPLACEMENT_Y,
        KratosMultiphysics.DISPLACEMENT_Z,
        KratosMultiphysics.ROTATION_X,
        KratosMultiphysics.ROTATION_Y,
        KratosMultiphysics.ROTATION_Z
    ]

    delta = 1e-8
    fd_sensitiivties = KratosMultiphysics.Vector(18, 0.0)
    local_index = 0
    for node in elem1.GetGeometry():
        for dof in dof_order:
            node.SetSolutionStepValue(dof, node.GetSolutionStepValue(dof) + delta)
            value = adjoint_response.CalculateValue(model_part)
            node.SetSolutionStepValue(dof, node.GetSolutionStepValue(dof) - delta)
            fd_sensitiivties[local_index] = (value - ref_value) / delta
            local_index += 1

    print(fd_sensitiivties)
    print(ad_sensitivities)
    print(np.linalg.norm(np.array(ad_sensitivities - fd_sensitiivties)))

def testAdjointDisp():
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("Structure")
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
    model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
    sub_model_part = model_part.CreateSubModelPart("Parts_AREAS")

    n1: KratosMultiphysics.Node = model_part.CreateNewNode(1, 0, 0, 0)
    n2: KratosMultiphysics.Node = model_part.CreateNewNode(2, 1, 0, 0)
    n3: KratosMultiphysics.Node = model_part.CreateNewNode(3, 1, 1, 0)

    n1.AddDof(KratosMultiphysics.DISPLACEMENT_X)
    n1.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
    n1.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
    n2.AddDof(KratosMultiphysics.DISPLACEMENT_X)
    n2.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
    n2.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
    n3.AddDof(KratosMultiphysics.DISPLACEMENT_X)
    n3.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
    n3.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

    n1.AddDof(KratosMultiphysics.ROTATION_X)
    n1.AddDof(KratosMultiphysics.ROTATION_Y)
    n1.AddDof(KratosMultiphysics.ROTATION_Z)
    n2.AddDof(KratosMultiphysics.ROTATION_X)
    n2.AddDof(KratosMultiphysics.ROTATION_Y)
    n2.AddDof(KratosMultiphysics.ROTATION_Z)
    n3.AddDof(KratosMultiphysics.ROTATION_X)
    n3.AddDof(KratosMultiphysics.ROTATION_Y)
    n3.AddDof(KratosMultiphysics.ROTATION_Z)

    n1.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)

    n1.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y)
    n1.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y)
    n2.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y)
    n3.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z)

    p1 = model_part.CreateNewProperties(1)
    material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
    material_settings["Parameters"]["materials_filename"].SetString("/software/kratos/digital_twin/applications/StructuralMechanicsApplication/tests/adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/StructuralMaterials.json")
    KratosMultiphysics.ReadMaterialsUtility(material_settings, model)

    elem1: KratosMultiphysics.Element = model_part.CreateNewElement("AdjointFiniteDifferencingShellThinElement3D3N", 1, [1, 2, 3], p1)
    sub_model_part.AddElements([1])

    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, node.Id)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, node.Id + 1)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, node.Id + 2)
        node.SetSolutionStepValue(KratosMultiphysics.ROTATION_X, node.Id)
        node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Y, node.Id + 1)
        node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Z, node.Id + 2)

    elem1.Initialize(model_part.ProcessInfo)
    residuals = KratosMultiphysics.Matrix()
    lhs = KratosMultiphysics.Matrix()
    rhs = KratosMultiphysics.Vector()
    elem1.CalculateLocalSystem(lhs, rhs, model_part.ProcessInfo)

    ref_gp_values = elem1.CalculateOnIntegrationPoints(StructuralMechanicsApplication.SHELL_STRAIN, model_part.ProcessInfo)
    adjoint_response = StructuralMechanicsApplication.AdjointNodalDispResponseFunction(model_part, KratosMultiphysics.Parameters("""{"weight": 2, "step_size":1e-9, "gradient_mode": "semi_analytic", "direction": [0, -1, 0], "node_id": 2}"""))
    adjoint_response.Initialize()

    ref_value = -n2.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)[1]

    ad_sensitivities = KratosMultiphysics.Vector()
    adjoint_response.CalculateGradient(elem1, lhs, ad_sensitivities, model_part.ProcessInfo)

    print(adjoint_response.CalculateValue(model_part), ref_value / 2)
    ref_value = adjoint_response.CalculateValue(model_part)

    dof_order = [
        KratosMultiphysics.DISPLACEMENT_X,
        KratosMultiphysics.DISPLACEMENT_Y,
        KratosMultiphysics.DISPLACEMENT_Z,
        KratosMultiphysics.ROTATION_X,
        KratosMultiphysics.ROTATION_Y,
        KratosMultiphysics.ROTATION_Z
    ]

    delta = 1e-8
    fd_sensitiivties = KratosMultiphysics.Vector(18, 0.0)
    local_index = 0
    for node in elem1.GetGeometry():
        for dof in dof_order:
            node.SetSolutionStepValue(dof, node.GetSolutionStepValue(dof) + delta)
            value = adjoint_response.CalculateValue(model_part)
            node.SetSolutionStepValue(dof, node.GetSolutionStepValue(dof) - delta)
            fd_sensitiivties[local_index] = (value - ref_value) / delta
            local_index += 1

    print(fd_sensitiivties)
    print(ad_sensitivities)
    print(np.linalg.norm(np.array(ad_sensitivities - fd_sensitiivties)))

if __name__ == '__main__':
    # testAdjointStrain()
    testAdjointDisp()
    # suites = KratosUnittest.KratosSuites
    # smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureLocalStress]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureNodalDisplacement]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureStrainEnergy]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisBeamStructureNodalReaction]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisShellStructureLocalStress]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisShellStructureNodalDisplacement]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisShellStructureStrainEnergy]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisSpringDamperElement]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisLinearTrussStructure]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityAnalysisNonLinearTrussStructure]))
    # allSuite = suites['all']
    # allSuite.addTests(smallSuite)
    # KratosUnittest.runTests(suites)
