import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import eigen_solver_factory
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestEigenSolvers(KratosUnittest.TestCase):

    def _RunParametrized(self, my_params_string, eigen_value_estimated = "lowest" ):
        all_settings = KratosMultiphysics.Parameters( my_params_string )

        for i in range(all_settings["test_list"].size()):
            settings = all_settings["test_list"][i]
            self._auxiliary_test_function(settings, "auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm", eigen_value_estimated)

    def _auxiliary_test_function(self, settings, matrix_name="auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm", eigen_value_estimated = "lowest"):
        space = KratosMultiphysics.UblasSparseSpace()

        # Read the matrices
        K = KratosMultiphysics.CompressedMatrix()
        KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath(matrix_name),K)

        n = K.Size1()
        M = KratosMultiphysics.CompressedMatrix(n, n)

        for i in range(n):
            for j in range(n):
                if (i == j):
                    M[i, j] = 1.0

        # create result containers (they will be resized inside the solver)
        eigenvalues = KratosMultiphysics.Vector(n)
        eigenvectors = KratosMultiphysics.Matrix(n, 1)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigen_solver.Solve(K, M, eigenvalues, eigenvectors)
        eigenvalue = eigenvalues[0]

        if (eigen_value_estimated == "lowest"):
            self.assertLessEqual(abs(eigenvalue - 0.061463)/0.061463, 5.0e-3)
        else:
            self.assertLessEqual(abs(eigenvalue - 11.959)/11.959, 5.0e-3)

    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_lowest_power_in_core(self):
        import platform

        if platform.system() != "Windows":
            self._RunParametrized("""
                {
                    "test_list" : [
                        {
                            "solver_type"             : "power_iteration_eigenvalue_solver",
                            "max_iteration"           : 10000,
                            "tolerance"               : 1e-8,
                            "required_eigen_number"   : 1,
                            "shifting_convergence"    : 0.25,
                            "verbosity"               : 0,
                            "linear_solver_settings"      : {
                                "solver_type"             : "LinearSolversApplication.sparse_lu",
                                "max_iteration"           : 500,
                                "tolerance"               : 1e-9,
                                "scaling"                 : false,
                                "verbosity"               : 0
                            }
                        }
                    ]
                }
                """)
        else:
            self.skipTest("Test temporaly disabled in Windows")

    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_highest_power_in_core(self):
        import platform

        if platform.system() != "Windows":
            self._RunParametrized("""
                {
                    "test_list" : [
                        {
                            "solver_type"             : "power_iteration_highest_eigenvalue_solver",
                            "max_iteration"           : 10000,
                            "tolerance"               : 1e-8,
                            "required_eigen_number"   : 1,
                            "shifting_convergence"    : 0.25,
                            "verbosity"               : 0,
                            "linear_solver_settings"      : {
                            "solver_type"             : "LinearSolversApplication.sparse_lu",
                                "max_iteration"           : 500,
                                "tolerance"               : 1e-9,
                                "scaling"                 : false,
                                "verbosity"               : 0
                            }
                        }
                    ]
                }
                """, "highest")
        else:
            self.skipTest("Test temporaly disabled in Windows")

    def test_rayleigh_in_core(self):
        import platform

        if platform.system() != "Windows":
            self._RunParametrized("""
                {
                    "test_list" : [
                        {
                            "solver_type"             : "rayleigh_quotient_iteration_eigenvalue_solver",
                            "max_iteration"           : 10000,
                            "tolerance"               : 1e-9,
                            "required_eigen_number"   : 1,
                            "shifting_convergence"    : 0.25,
                            "verbosity"               : 0,
                            "linear_solver_settings"      : {
                                "solver_type"             : "skyline_lu_factorization",
                                "max_iteration"           : 500,
                                "tolerance"               : 1e-9,
                                "scaling"                 : false,
                                "verbosity"               : 0
                            }
                        }
                    ]
                }
                """)
        else:
            self.skipTest("Test temporaly disabled in Windows")

    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_eigen_eigensystem_solver(self):
        import platform

        if platform.system() != "Windows":
            self._RunParametrized("""
                {
                    "test_list" : [
                        {
                            "solver_type": "eigen_eigensystem",
                            "number_of_eigenvalues": 3,
                            "max_iteration": 1000,
                            "tolerance": 1e-8,
                            "echo_level": 1
                        }
                    ]
                }
                """)
        else:
            self.skipTest("Test temporaly disabled in Windows")

    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_FEAST_with_eigen_solver(self):
        import platform

        if platform.system() != "Windows":
            from KratosMultiphysics import LinearSolversApplication
            if not LinearSolversApplication.HasFEAST():
                self.skipTest("FEAST is not available")
            self._RunParametrized("""
                {
                    "test_list" : [
                        {
                            "solver_type": "feast",
                            "symmetric": true,
                            "e_min": 0.01,
                            "e_max": 0.20,
                            "subspace_size": 5
                        }
                    ]
                }
                """)
        else:
            self.skipTest("Test temporaly disabled in Windows")

if __name__ == '__main__':
    KratosUnittest.main()
