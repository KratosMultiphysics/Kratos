import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import eigen_solver_factory
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestConditionNumber(KratosUnittest.TestCase):

    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_condition_number(self):
        # Read the matrices
        K = KratosMultiphysics.CompressedMatrix()
        KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath("auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm"),K)

        # Construct the solver
        settings_max = KratosMultiphysics.Parameters("""
        {
            "solver_type"             : "power_iteration_highest_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "LinearSolversApplication.sparse_lu",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        """)
        eigen_solver_max = eigen_solver_factory.ConstructSolver(settings_max)
        settings_min = KratosMultiphysics.Parameters("""
        {
            "solver_type"             : "power_iteration_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "LinearSolversApplication.sparse_lu",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        """)
        eigen_solver_min = eigen_solver_factory.ConstructSolver(settings_min)

        # Solve
        condition_number_utility = KratosMultiphysics.ConditionNumberUtility()
        condition_number = condition_number_utility.GetConditionNumber(K, eigen_solver_max, eigen_solver_min)

        self.assertLessEqual(abs(condition_number-194.5739)/194.5739, 1e-3)

if __name__ == '__main__':
    KratosUnittest.main()
