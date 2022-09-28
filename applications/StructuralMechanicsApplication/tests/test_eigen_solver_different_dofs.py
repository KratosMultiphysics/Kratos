import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as KSM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

'''
Test description:
This test does an eigenvalue analysis on a simple structure
that contains elements with different Dofs,
Beams with DISPLACEMENT and ROTATION dofs and Trusses with only DISPLACEMENT dofs
'''

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestEigenSolverWithDifferentDofs(KratosUnittest.TestCase):
    # muting the output
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

    def test_eigen_with_different_dofs_block_builder(self):
        self.execute_test_eigen_with_different_dofs(use_block_builder=True)

    def test_eigen_with_different_dofs_elimination_builder(self):
        self.execute_test_eigen_with_different_dofs(use_block_builder=False)

    def execute_test_eigen_with_different_dofs(self, use_block_builder):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open(GetFilePath("eigen_test/Eigen_different_dofs_parameters.json"),'r') as parameter_file:
                parameters = KM.Parameters(parameter_file.read())

            parameters["solver_settings"]["block_builder"].SetBool(use_block_builder)

            model = KM.Model()
            StructuralMechanicsAnalysis(model, parameters).Run()
            self.__CheckEigenSolution(model["Structure"])

    def __CheckEigenSolution(self, model_part):
        exp_eigen_values = KM.Vector([2.3779548, 2.465771262])

        self.assertVectorAlmostEqual(exp_eigen_values, model_part.ProcessInfo[KSM.EIGENVALUE_VECTOR])

        eigen_vec_res_bottom_nodes = KM.Matrix(2,6)
        eigen_vec_res_bottom_nodes.fill(0.0)

        eigen_vec_res_top_node_raw = [
            [0.0, -0.00579401, 0.0,-0.9999832,0.0, 0.0],
            [-0.00716752, 0.0, 0.0, 0.0, 0.99997431, 0.0]
        ]

        eigen_vec_res_top_node = KM.Matrix(2,6)
        for i in range(eigen_vec_res_top_node.Size1()):
            for j in range(eigen_vec_res_top_node.Size2()):
                eigen_vec_res_top_node[i,j] = eigen_vec_res_top_node_raw[i][j]

        for node in model_part.Nodes:
            if node.Id == 1: # this is the node at the top
                self.__CompareMatrixAbs(eigen_vec_res_top_node, node[KSM.EIGENVECTOR_MATRIX])
            else:
                # here we make sure that the fixed dofs have a zero eigenvector
                self.assertMatrixAlmostEqual(eigen_vec_res_bottom_nodes, node[KSM.EIGENVECTOR_MATRIX])

    def __CompareMatrixAbs(self, mat_1, mat_2, tol=7):
        self.assertEqual(mat_1.Size1(), mat_2.Size1())
        self.assertEqual(mat_1.Size2(), mat_2.Size2())

        for i in range(mat_1.Size1()):
            for j in range(mat_1.Size2()):
                self.assertAlmostEqual(abs(mat_1[i,j]), abs(mat_2[i,j]), tol)


if __name__ == '__main__':
    KratosUnittest.main()
