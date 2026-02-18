# Modules
import os
import math
import pytest

# Kratos
from KratosMultiphysics import Model, Parameters, Array3, Logger                        # Classes
from KratosMultiphysics import VariableUtils, NormalCalculationUtils, MortarUtilities   # Utilities
from KratosMultiphysics import IsDistributedRun                                         # Parell Interface
from KratosMultiphysics import (
    NORMAL, NODAL_VAUX, DOMAIN_SIZE, NORMAL_SHAPE_DERIVATIVE,
    DOMAIN_SIZE, BULK_MODULUS, EXTERNAL_FORCES_VECTOR, LOCAL_AXES_MATRIX
)                                                                                       # Variables

# Kratos Modules
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def PostProcess(model_part):
    gid_output = GiDOutputProcess(model_part,
        "gid_output",
        Parameters("""
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
    normal = Array3([node.X/norm, node.Y/norm, node.Z/norm])

    return normal

def CalculateNorm(array_3d_value):
    return math.sqrt(array_3d_value[0]**2+array_3d_value[1]**2+array_3d_value[2]**2)

def RemoveFiles(mdpa_name):
    DeleteFileIfExisting(mdpa_name + ".time")

def FiniteDifferenceNormalShapeSensitivity(model_part, check_condition_ids, delta, tolerance):
    dimensionality = model_part.ProcessInfo[DOMAIN_SIZE]
    normal_calculation_utils = NormalCalculationUtils()

    ## calculate analytical shape derivatives for all conditions
    normal_calculation_utils.CalculateNormalShapeDerivativesOnSimplex(model_part.Conditions, dimensionality)
    def get_matrix_row(matrix, row):
        v = Array3(0.0)
        for i in range(matrix.Size2()):
            v[i] = matrix[row, i]
        return v

    ## calculate finite difference shape derivatives
    normal_calculation_utils.CalculateOnSimplex(model_part.Conditions, dimensionality)
    reference_normals = {}
    for condition in model_part.Conditions:
        if (condition.Id in check_condition_ids):
            reference_normals[condition.Id] = condition.GetValue(NORMAL)

    def check_normal_sensitivities(condition, row_index, analytical_normal_shape_sensitivities):
        normal_calculation_utils.CalculateOnSimplex(model_part.Conditions, dimensionality)
        perturbed_normal = condition.GetValue(NORMAL)
        
        dif = perturbed_normal - reference_normals[condition.Id]
        ref = get_matrix_row(analytical_normal_shape_sensitivities, row_index)
        
        for i in range(3):
            assert (dif[i]/delta) == pytest.approx(ref[i], abs=tolerance)

        # UnitTestObject.assertVectorAlmostEqual((perturbed_normal - reference_normals[condition.Id]) / delta, get_matrix_row(
        #     analytical_normal_shape_sensitivities, row_index), tolerance)

    for condition in model_part.Conditions:
        if (condition.Id in check_condition_ids):
            analytical_normal_shape_sensitivities = condition.GetValue(
                NORMAL_SHAPE_DERIVATIVE)
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

@pytest.fixture(scope="class", autouse=True)
def test_model():
    model = Model()

    return model

class TestNormalUtilsEmptyModelPart():
    @pytest.fixture(scope="class", autouse=True)
    def test_modelpart(self, test_model):
        modelpart_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions")
        modelpart      = test_model.CreateModelPart("Main")

        modelpart.ProcessInfo[DOMAIN_SIZE] = 3
        modelpart.AddNodalSolutionStepVariable(NORMAL)
        modelpart.AddNodalSolutionStepVariable(BULK_MODULUS)
        modelpart.AddNodalSolutionStepVariable(NODAL_VAUX)
        modelpart.AddNodalSolutionStepVariable(EXTERNAL_FORCES_VECTOR)
        modelpart.AddNodalSolutionStepVariable(LOCAL_AXES_MATRIX)
        
        ReadModelPart(modelpart_name, modelpart)
        VariableUtils().SetHistoricalVariableToZero(NORMAL, modelpart.Nodes)

        yield modelpart
        
        RemoveFiles(modelpart_name)

    @pytest.fixture(autouse=True)
    def test_clean_modelpart(self, test_modelpart):
        VariableUtils().SetHistoricalVariableToZero(NORMAL, test_modelpart.Nodes)

        return test_modelpart

    def test_ComputeSimplexNormalModelPart(self, test_clean_modelpart):
        NormalCalculationUtils().CalculateOnSimplex(test_clean_modelpart)

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(NORMAL)
            solution_normal /= CalculateNorm(solution_normal)

            assert CalculateNorm(normal - solution_normal) < 0.15

    @pytest.mark.skipif(IsDistributedRun(), "This test is designed for serial runs only.")
    def test_ComputeSimplexNormalModelPartWithLineCondition(self, test_clean_modelpart):
        # Adding one line, to make sure it is getting ignored
        
        # CHARLIE: This is potentially dangerous both with the new pytest and the unittest implementation as both
        # treat the modelpart as a class variable, meaning that this line will be added to the model for all
        # furhter tests. 
        test_clean_modelpart.CreateNewCondition("LineCondition3D2N", 1000, [1,2], test_clean_modelpart.GetProperties()[1])

        NormalCalculationUtils().CalculateOnSimplex(test_clean_modelpart)
        
        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(NORMAL)
            solution_normal /= CalculateNorm(solution_normal)

            assert CalculateNorm(normal - solution_normal) < 0.15

    def test_ComputeSimplexNormalModelPartWithCustomVariable(self, test_clean_modelpart):
        # Calculate the normals using NODAL_VAUX as custom variable
        NormalCalculationUtils().CalculateOnSimplex( test_clean_modelpart, NODAL_VAUX)

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(NODAL_VAUX)
            solution_normal /= CalculateNorm(solution_normal)

            assert CalculateNorm(normal - solution_normal) < 0.15

    def test_ComputeSimplexNormalModelPartNonHistorical(self, test_clean_modelpart):
        NormalCalculationUtils().CalculateOnSimplexNonHistorical(test_clean_modelpart)

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetValue(NORMAL)
            solution_normal /= CalculateNorm(solution_normal)

            assert CalculateNorm(normal - solution_normal) < 0.15

    def test_ComputeUnitNormalModelPart(self, test_clean_modelpart):
        NormalCalculationUtils().CalculateUnitNormals(test_clean_modelpart)

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            
            solution_normal = node.GetSolutionStepValue(NORMAL)
            
            assert CalculateNorm(normal - solution_normal) < 0.15

    def test_ComputeUnitNormalModelPartWithCustomVariable(self, test_clean_modelpart):
        # Calculate the unit normals using NODAL_VAUX as custom variable
        enforce_generic_algorithm = False
        NormalCalculationUtils().CalculateUnitNormals(
            test_clean_modelpart,
            enforce_generic_algorithm,
            NODAL_VAUX
        )

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(NODAL_VAUX)
            assert CalculateNorm(normal - solution_normal) < 0.15

    def test_ComputeUnitNormalModelPartNonHistorical(self, test_clean_modelpart):
        # Calculate the unit normals using NODAL_VAUX as custom variable
        enforce_generic_algorithm = False
        NormalCalculationUtils().CalculateUnitNormalsNonHistorical(
            test_clean_modelpart,
            enforce_generic_algorithm
        )

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetValue(NORMAL)

            assert CalculateNorm(normal - solution_normal) < 0.15

    def test_ComputeNodesMeanNormalModelPart(self, test_clean_modelpart):
        NormalCalculationUtils().CalculateUnitNormals(test_clean_modelpart, True)

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            
            solution_normal = node.GetSolutionStepValue(NORMAL)
            
            assert CalculateNorm(normal - solution_normal) < 0.1

    def test_ComputeNodesMeanNormalModelPartWithCustomVariable(self, test_clean_modelpart):
        # Calculate the unit normals using NODAL_VAUX as custom variable
        enforce_generic_algorithm = True
        NormalCalculationUtils().CalculateUnitNormals(
            test_clean_modelpart,
            enforce_generic_algorithm,
            NODAL_VAUX
        )

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(NODAL_VAUX)

            assert CalculateNorm(normal - solution_normal) < 0.1

    def test_ComputeNodesMeanNormalModelPartNonHistorical(self, test_clean_modelpart):
        enforce_generic_algorithm = True
        NormalCalculationUtils().CalculateUnitNormalsNonHistorical(
            test_clean_modelpart,
            enforce_generic_algorithm
        )

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetValue(NORMAL)

            assert CalculateNorm(normal - solution_normal) < 0.1

    def test_InvertNormal(self, test_clean_modelpart):
        MortarUtilities.InvertNormal(test_clean_modelpart.Conditions)
        NormalCalculationUtils().CalculateUnitNormals(test_clean_modelpart, True)

        for node in test_clean_modelpart.GetSubModelPart("Skin_Part").Nodes:

            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(NORMAL) * -1.0

            assert CalculateNorm(normal - solution_normal) < 0.1

    def test_CalculateShapeDerivativesOnSimplexConditions3D(self, test_clean_modelpart):
        FiniteDifferenceNormalShapeSensitivity(test_clean_modelpart, [1, 2, 6, 10, 15, 20, 25, 30, 31, 32, 50, 80, 100, 120], 1e-8, 8)

class TestNormalUtilsQuadSphere():
    @pytest.fixture(scope="class", autouse=True)
    def test_modelpart(self, test_model):
        modelpart_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/quad_sphere")
        modelpart      = test_model.CreateModelPart("Main")

        modelpart.ProcessInfo[DOMAIN_SIZE] = 3
        modelpart.AddNodalSolutionStepVariable(NORMAL)
        
        ReadModelPart(modelpart_name, modelpart)
        VariableUtils().SetHistoricalVariableToZero(NORMAL, modelpart.Nodes)

        yield modelpart
        
        RemoveFiles(modelpart_name)

    @pytest.fixture(autouse=True)
    def test_clean_modelpart(self, test_modelpart):
        VariableUtils().SetHistoricalVariableToZero(NORMAL, test_modelpart.Nodes)

        return test_modelpart

    def test_ComputeUnitNormalQuadModelPart(self, test_clean_modelpart):
        NormalCalculationUtils().CalculateUnitNormals(test_clean_modelpart)

        for node in test_clean_modelpart.Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(NORMAL)

            assert CalculateNorm(normal - solution_normal) < 0.15

class TestNormalUtils2DSymmetricalSquare():
    @pytest.fixture(scope="class", autouse=True)
    def test_modelpart(test_model):
        modelpart_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square")
        modelpart      = test_model.CreateModelPart("Main")

        modelpart.ProcessInfo[DOMAIN_SIZE] = 2
        modelpart.AddNodalSolutionStepVariable(NORMAL)
        
        ReadModelPart(modelpart_name, modelpart)
        VariableUtils().SetHistoricalVariableToZero(NORMAL, modelpart.Nodes)

        yield modelpart
        
        RemoveFiles(modelpart_name)

    @pytest.fixture(autouse=True)
    def test_clean_modelpart(test_modelpart):
        VariableUtils().SetHistoricalVariableToZero(NORMAL, test_modelpart.Nodes)

        return test_modelpart

    def test_CalculateShapeDerivativesOnSimplexConditions2D(self, test_clean_modelpart):
        FiniteDifferenceNormalShapeSensitivity(test_clean_modelpart, [1, 2, 3, 4, 5, 10, 11, 12, 13, 19, 20], 1e-7, 8)
