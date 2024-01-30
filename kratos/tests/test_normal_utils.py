import os
import math

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def PostProcess(model_part):
    gid_output = GiDOutputProcess(model_part,
                                "gid_output",
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "nodal_results" : ["NORMAL"]
                                        }
                                    }
                                    """)
                                )

    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

def CalculateAnalyticalNormal(node):
    norm = math.sqrt(node.X**2+node.Y**2+node.Z**2)
    normal = KratosMultiphysics.Array3([node.X/norm, node.Y/norm, node.Z/norm])
    return normal

def CalculateNorm(array_3d_value):
    return math.sqrt(array_3d_value[0]**2+array_3d_value[1]**2+array_3d_value[2]**2)

def RemoveFiles(mdpa_name):
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

def FiniteDifferenceNormalShapeSensitivityTest(UnitTestObject, model_part, check_condition_ids, delta, tolerance):
    dimensionality = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
    normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()

    ## calculate analytical shape derivatives for all conditions
    normal_calculation_utils.CalculateNormalShapeDerivativesOnSimplex(model_part.Conditions, dimensionality)
    def get_matrix_row(matrix, row):
        v = KratosMultiphysics.Array3(0.0)
        for i in range(matrix.Size2()):
            v[i] = matrix[row, i]
        return v

    ## calculate finite difference shape derivatives
    normal_calculation_utils.CalculateOnSimplex(model_part.Conditions, dimensionality)
    reference_normals = {}
    for condition in model_part.Conditions:
        if (condition.Id in check_condition_ids):
            reference_normals[condition.Id] = condition.GetValue(KratosMultiphysics.NORMAL)

    def check_normal_sensitivities(condition, row_index, analytical_normal_shape_sensitivities):
        normal_calculation_utils.CalculateOnSimplex(
            model_part.Conditions, dimensionality)
        perturbed_normal = condition.GetValue(
            KratosMultiphysics.NORMAL)
        UnitTestObject.assertVectorAlmostEqual((perturbed_normal - reference_normals[condition.Id]) / delta, get_matrix_row(
            analytical_normal_shape_sensitivities, row_index), tolerance)

    for condition in model_part.Conditions:
        if (condition.Id in check_condition_ids):
            analytical_normal_shape_sensitivities = condition.GetValue(
                KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE)
            for i, node in enumerate(condition.GetGeometry()):
                node.X += delta
                check_normal_sensitivities(
                    condition, i * dimensionality, analytical_normal_shape_sensitivities)
                node.X -= delta

                node.Y += delta
                check_normal_sensitivities(
                    condition, i * dimensionality + 1, analytical_normal_shape_sensitivities)
                node.Y -= delta

                if dimensionality == 3:
                    node.Z += delta
                    check_normal_sensitivities(
                        condition, i * dimensionality + 2, analytical_normal_shape_sensitivities)
                    node.Z -= delta


class TestNormalUtilsCoarseSphere(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCES_VECTOR)
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LOCAL_AXES_MATRIX)
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)

    def test_ComputeSimplexNormalModelPart(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            solution_normal /= CalculateNorm(solution_normal)

            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

    @KratosUnittest.skipIf(KratosMultiphysics.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_ComputeSimplexNormalModelPartWithLineCondition(self):
        #Adding one line, to make sure it is getting ignored
        self.model_part.CreateNewCondition("LineCondition3D2N", 1000, [1,2], self.model_part.GetProperties()[1])
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part)
        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            solution_normal /= CalculateNorm(solution_normal)

            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

    def test_ComputeSimplexNormalModelPartWithCustomVariable(self):
        ## Calculate the normals using NODAL_VAUX as custom variable
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(
            self.model_part,
            KratosMultiphysics.NODAL_VAUX)

        ## DEBUG
        #PostProcess(self.model_part)

        ## Check results
        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NODAL_VAUX)
            solution_normal /= CalculateNorm(solution_normal)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

    def test_ComputeSimplexNormalModelPartNonHistorical(self):
        ## Calculate the normals using NODAL_VAUX as custom variable
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(self.model_part)

        ## DEBUG
        #PostProcess(self.model_part)

        ## Check results
        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetValue(KratosMultiphysics.NORMAL)
            solution_normal /= CalculateNorm(solution_normal)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

    def test_ComputeUnitNormalModelPart(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(self.model_part)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

    def test_ComputeUnitNormalModelPartWithCustomVariable(self):
        ## Calculate the unit normals using NODAL_VAUX as custom variable
        enforce_generic_algorithm = False
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(
            self.model_part,
            enforce_generic_algorithm,
            KratosMultiphysics.NODAL_VAUX)

        ## DEBUG
        #PostProcess(self.model_part)

        ## Check results
        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NODAL_VAUX)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

    def test_ComputeUnitNormalModelPartNonHistorical(self):
        ## Calculate the unit normals using NODAL_VAUX as custom variable
        enforce_generic_algorithm = False
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormalsNonHistorical(
            self.model_part,
            enforce_generic_algorithm)

        ## DEBUG
        #PostProcess(self.model_part)

        ## Check results
        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

    def test_ComputeNodesMeanNormalModelPart(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(self.model_part, True)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.1)

    def test_ComputeNodesMeanNormalModelPartWithCustomVariable(self):
        ## Calculate the unit normals using NODAL_VAUX as custom variable
        enforce_generic_algorithm = True
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(
            self.model_part,
            enforce_generic_algorithm,
            KratosMultiphysics.NODAL_VAUX)

        ## DEBUG
        #PostProcess(self.model_part)

        ## Check results
        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NODAL_VAUX)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.1)

    def test_ComputeNodesMeanNormalModelPartNonHistorical(self):
        ## Calculate the unit normals using NODAL_VAUX as custom variable
        enforce_generic_algorithm = True
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormalsNonHistorical(
            self.model_part,
            enforce_generic_algorithm)

        ## DEBUG
        #PostProcess(self.model_part)

        ## Check results
        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.1)

    def test_InvertNormal(self):
        KratosMultiphysics.MortarUtilities.InvertNormal(self.model_part.Conditions)
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(self.model_part, True)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL) * -1.0
            self.assertLess(CalculateNorm(normal - solution_normal), 0.1)

    def test_CalculateShapeDerivativesOnSimplexConditions3D(self):
        FiniteDifferenceNormalShapeSensitivityTest(self, self.model_part, [1, 2, 6, 10, 15, 20, 25, 30, 31, 32, 50, 80, 100, 120], 1e-8, 8)

class TestNormalUtilsQuadSphere(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/quad_sphere")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)

    def test_ComputeUnitNormalQuadModelPart(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(self.model_part)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNorm(normal - solution_normal), 0.15)

class TestNormalUtils2DSymmetricalSquare(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)

    def test_CalculateShapeDerivativesOnSimplexConditions2D(self):
        FiniteDifferenceNormalShapeSensitivityTest(self, self.model_part, [1, 2, 3, 4, 5, 10, 11, 12, 13, 19, 20], 1e-7, 8)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
