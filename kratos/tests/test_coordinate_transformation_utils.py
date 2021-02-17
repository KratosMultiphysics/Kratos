import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def SetupModelPart(cls, mdpa_name, domain_size):
    cls.current_model = KratosMultiphysics.Model()
    cls.model_part = cls.current_model.CreateModelPart("Main")
    cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = domain_size
    cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
    cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
    cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
    cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCES_VECTOR)
    cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LOCAL_AXES_MATRIX)
    cls.mdpa_name = GetFilePath(mdpa_name)
    ReadModelPart(cls.mdpa_name, cls.model_part)

    KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, False, cls.model_part.Nodes)

    for condition in cls.model_part.Conditions:
        condition.Set(KratosMultiphysics.SLIP, True)
        for node in condition.GetGeometry():
            node.Set(KratosMultiphysics.SLIP, True)

    cls.model_part.GetCommunicator().SynchronizeOrNodalFlags(KratosMultiphysics.SLIP)

def RemoveFiles(mdpa_name):
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

def FiniteDifferenceNormalShapeSensitivityTest(UnitTestObject, model_part, check_node_ids, delta, tolerance):
    domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

    # find nodal neighbour id map
    process = KratosMultiphysics.FindGlobalNodalNeighboursForConditionsProcess(model_part.GetCommunicator().GetDataCommunicator(), model_part)
    process.Execute()
    neighbour_node_id_map = process.GetNeighbourIds(model_part.Nodes)

    # calculate nodal normal shape sensitivities
    KratosMultiphysics.NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(model_part.Conditions, domain_size)
    KratosMultiphysics.SensitivityUtilities.AssignConditionDerivativesToNodes(
        model_part,
        domain_size,
        KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE,
        neighbour_node_id_map,
        1.0/domain_size,
        KratosMultiphysics.SLIP,
        True)

    ref_values = {}
    coordinate_transformation_utils = KratosMultiphysics.CoordinateTransformationUtils(domain_size, domain_size, KratosMultiphysics.SLIP)
    # calculate nodal normals
    KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(model_part)

    # calculate nodal rotation matrices
    for node in model_part.Nodes:
        if (node.Is(KratosMultiphysics.SLIP) and node.Id in check_node_ids):
            ref_values[node.Id] = KratosMultiphysics.Matrix()
            coordinate_transformation_utils.CalculateRotationOperatorPure(ref_values[node.Id], node)

    def CheckSensitivities(node, derivative_node_index, derivative_direction_index, ref_value):
        # calculate analytical shape sensitivities
        analytical_sensitivity = KratosMultiphysics.Matrix()
        coordinate_transformation_utils.CalculateRotationOperatorPureShapeSensitivities(analytical_sensitivity, derivative_node_index, derivative_direction_index, node)

        # calculate fd shape sensitivities
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(model_part)
        perturbed_rotation_matrix = KratosMultiphysics.Matrix()
        coordinate_transformation_utils.CalculateRotationOperatorPure(perturbed_rotation_matrix, node)
        fd_shape_sensitivity = (perturbed_rotation_matrix - ref_value) / delta

        UnitTestObject.assertMatrixAlmostEqual(fd_shape_sensitivity, analytical_sensitivity, tolerance)

    for ref_id, ref_value in ref_values.items():
        base_node = model_part.GetNode(ref_id)
        derivative_ids = [base_node.Id]
        derivative_ids.extend(neighbour_node_id_map[base_node.Id])

        for i, derivative_node_id in enumerate(derivative_ids):
            derivative_node = model_part.GetNode(derivative_node_id)

            derivative_node.X += delta
            CheckSensitivities(base_node, i, 0, ref_value)
            derivative_node.X -= delta

            derivative_node.Y += delta
            CheckSensitivities(base_node, i, 1, ref_value)
            derivative_node.Y -= delta

            if domain_size == 3:
                derivative_node.Z += delta
                CheckSensitivities(base_node, i, 2, ref_value)
                derivative_node.Z -= delta

class TestCoordinateTransformationUtilitiesCoarseSphere(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        SetupModelPart(cls, "auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions", 3)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)

    def test_CalculateRotationOperatorPureShapeSensitivities(self):
        check_node_ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 50, 70, 80, 85]
        FiniteDifferenceNormalShapeSensitivityTest(self, self.model_part, check_node_ids, 1e-8, 6)


class TestCoordinateTransformationUtilities2DSymmetricalSquare(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        SetupModelPart(cls, "auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square", 2)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)

    def test_CalculateRotationOperatorPureShapeSensitivities(self):
        check_node_ids = [1, 2, 3, 4, 8, 9, 30, 60, 61]
        FiniteDifferenceNormalShapeSensitivityTest(self, self.model_part, check_node_ids, 1e-8, 6)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

