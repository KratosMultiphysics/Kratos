from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName


class TestEigenSolvers(KratosUnittest.TestCase):

    def _RunParametrized(self, my_params_string, eigen_value_estimated = "lowest" ):
        all_settings = KratosMultiphysics.Parameters( my_params_string )

        for i in range(all_settings["test_list"].size()):
            settings = all_settings["test_list"][i]
            self._auxiliary_test_function(settings, "A.mm", eigen_value_estimated)

    def _auxiliary_test_function(self, settings, matrix_name="A.mm", eigen_value_estimated = "lowest"):
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
                else:
                    M[i, j] = 0.0

        # Construct the solver
        import eigen_solver_factory
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigenvalue = eigen_solver.GetEigenValue(K, M)

        if (eigen_value_estimated == "lowest"):
            self.assertLessEqual(abs(eigenvalue - 0.061463)/0.061463, 5.0e-3)
        else:
            self.assertLessEqual(abs(eigenvalue - 11.959)/11.959, 5.0e-3)

    def test_lowest_power_in_core(self):
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
                            "solver_type"             : "SuperLUSolver",
                            "max_iteration"           : 500,
                            "tolerance"               : 1e-9,
                            "scaling"                 : false,
                            "verbosity"               : 0
                        }
                    }
                ]
            }
            """)

    def test_highest_power_in_core(self):
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
                        "solver_type"             : "SuperLUSolver",
                            "max_iteration"           : 500,
                            "tolerance"               : 1e-9,
                            "scaling"                 : false,
                            "verbosity"               : 0
                        }
                    }
                ]
            }
            """, "highest")

    def test_rayleigh_in_core(self):
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
                            "solver_type"             : "SkylineLUFactorizationSolver",
                            "max_iteration"           : 500,
                            "tolerance"               : 1e-9,
                            "scaling"                 : false,
                            "verbosity"               : 0
                        }
                    }
                ]
            }
            """)

    def test_eigen_sparse_eigensystem_solver(self):
        try:
            import KratosMultiphysics.EigenSolversApplication
        except:
            self.skipTest("KratosMultiphysics.EigenSolversApplication is not available")

            self._RunParametrized("""
                {
                    "test_list" : [
                        {
                            "solver_type": "eigen_sparse_eigensystem",
                            "number_of_eigenvalues": 1,
                            "max_iteration": 100,
                            "tolerance": 1e-8,
                            "echo_level": 1
                        }
                    ]
                }
                """)

if __name__ == '__main__':
    KratosUnittest.main()
