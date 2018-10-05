from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestLinearSolvers(KratosUnittest.TestCase):

    def _RunParametrized(self, my_params_string ):
        all_settings = KratosMultiphysics.Parameters( my_params_string )

        for i in range(all_settings["test_list"].size()):
            settings = all_settings["test_list"][i]
            self._auxiliary_test_function(settings)

    def _auxiliary_test_function(self, settings, matrix_name="A.mm"):
        space = KratosMultiphysics.UblasSparseSpace()

        #read the matrices
        A = KratosMultiphysics.CompressedMatrix()
        KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath(matrix_name),A)

        Aoriginal = KratosMultiphysics.CompressedMatrix(A) #create a copy of A

        n = A.Size1()
        b = KratosMultiphysics.Vector(n)
        space.SetToZeroVector(b)

        for i in range(len(b)):
            b[i] = i/len(b)

        x = KratosMultiphysics.Vector(n)
        #KratosMultiphysics.ReadMatrixMarketVector("b.mm",b)

        boriginal = KratosMultiphysics.Vector(b) #create a copy of b

        space.SetToZeroVector(x)
        #space.SetToZeroVector(boriginal)
        #space.UnaliasedAdd(boriginal, 1.0, b) #boriginal=1*bs

        #construct the solver
        import new_linear_solver_factory
        linear_solver = new_linear_solver_factory.ConstructSolver(settings)

        #solve
        linear_solver.Solve(A,x,b)

        #test the results
        tmp = KratosMultiphysics.Vector(n)
        tmp *= 0.0
        space.Mult(Aoriginal,x,tmp)

        check = KratosMultiphysics.Vector(n)
        check = boriginal - tmp

        achieved_norm = space.TwoNorm(check)

        tolerance = 1e-9
        if(settings.Has("tolerance")):
            tolerance = settings["tolerance"].GetDouble()

        target_norm = tolerance*space.TwoNorm(boriginal)

        if(not (achieved_norm <= target_norm)):
            KratosMultiphysics.Logger.PrintInfo("Test linear solvers: ", "Echo of settings for failing test:\n", settings.PrettyPrintJsonString())
            KratosMultiphysics.Logger.PrintInfo("Test linear solvers: ", "Achieved_norm",achieved_norm, "\n", "Target_norm", target_norm)
        self.assertTrue(achieved_norm <= target_norm)


    def test_tfqmr_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "TFQMRSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "ILU0Preconditioner"
                    },
                    {
                        "solver_type" : "TFQMRSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "DiagonalPreconditioner"
                    },
                    {
                        "solver_type" : "TFQMRSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 1000,
                        "preconditioner_type" : "None"
                    }
                ]
            }
            """)

    def test_cg_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "CGSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "DiagonalPreconditioner"
                    },
                    {
                        "solver_type" : "CGSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 1000,
                        "preconditioner_type" : "None"
                    }
                ]
            }
            """)

    def test_deflated_cg_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "DeflatedCGSolver"
                    }
                ]
            }
            """)

    def test_bicgstab_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "BICGSTABSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "ILU0Preconditioner",
                        "scaling": false
                    },
                    {
                        "solver_type" : "BICGSTABSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "ILU0Preconditioner",
                        "scaling": false
                    },
                    {
                        "solver_type" : "BICGSTABSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "DiagonalPreconditioner"
                    },
                    {
                        "solver_type" : "BICGSTABSolver",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "None"
                    }
                ]
            }
            """)

    def test_skyline_lu(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "SkylineLUFactorizationSolver",
                        "scaling": false
                    }
                ]
            }
            """)

    def test_superlu(self):
        try:
            import KratosMultiphysics.ExternalSolversApplication
        except:
            self.skipTest("KratosMultiphysics.ExternalSolversApplication is not available")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "SuperLUSolver",
                        "scaling": false
                    },
                    {
                        "solver_type" : "SuperLUIterativeSolver",
                        "scaling": false
                    },
                    {
                        "solver_type" : "SuperLUIterativeSolver",
                        "scaling": true
                    }
                ]
            }
            """)

    ##@KratosUnittest.skipUnless(hasattr(KratosMultiphysics,  "PastixSolver"), "Pastix solver is not included in the compilation of the External Solvers Application")
    def test_pastix(self):
        try:
            import KratosMultiphysics.ExternalSolversApplication
        except:
            self.skipTest("ExternalSolversApplication is not available")

        if( not hasattr(KratosMultiphysics.ExternalSolversApplication,  "PastixSolver") ):
            self.skipTest("Pastix solver is not included in the compilation of the External Solvers Application")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "PastixSolver",
                        "solution_method": "Direct",
                            "tolerance":1e-6,
                            "max_iteration":100,
                            "gmres_krylov_space_dimension":100,
                            "ilu_level_of_fill":1,
                            "is_symmetric":false,
                            "verbosity":0,
                            "scaling": false,
                            "block_size": 1,
                        "use_block_matrices_if_possible" : true
                    }
                ]
            }
            """)

    def test_bicgstab_iluk(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {

                        "solver_type" : "AMGCL",
                        "smoother_type":"iluk",
                        "krylov_type": "bicgstab",
                        "coarsening_type": "aggregation",
                        "max_iteration": 100,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 2,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                    }
                ]
            }
            """)

    def test_lgmres_iluk(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {

                        "solver_type" : "AMGCL",
                        "smoother_type":"iluk",
                        "krylov_type": "lgmres",
                        "coarsening_type": "aggregation",
                        "max_iteration": 300,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 1 ,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                    }
                ]
            }
            """)

    def test_amgcl_unpreconditioned(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "bicgstab",
                        "preconditioner_type": "dummy",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "gmres",
                        "preconditioner_type": "dummy",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "preconditioner_type": "dummy",
                        "verbosity" : 1
                    }
                ]
            }
            """)

    def test_amgcl_no_amg_only_preconditioner(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "smoother_type":"ilu0",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "smoother_type":"ilu0",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1,
                        "block_size" : 2
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "smoother_type":"iluk",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "smoother_type":"spai0",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "smoother_type":"damped_jacobi",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "smoother_type":"gauss_seidel",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "AMGCL",
                        "krylov_type": "lgmres",
                        "smoother_type":"chebyshev",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    }
                ]
            }
            """)

    def test_amgcl_bicgstab_ilu0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "AMGCL",
                        "smoother_type":"ilu0",
                        "krylov_type": "bicgstab",
                        "coarsening_type": "aggregation",
                        "max_iteration": 100,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 0,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                    }
                ]
            }
            """)

    def test_amgcl_idr_ilu0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "AMGCL",
                        "smoother_type":"ilu0",
                        "krylov_type": "idrs",
                        "coarsening_type": "aggregation",
                        "max_iteration": 100,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 0,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                    }
                ]
            }
            """)

    def test_amgcl_bicgstab_spai0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "AMGCL",
                        "smoother_type":"spai0",
                        "krylov_type": "bicgstab",
                        "coarsening_type": "aggregation",
                        "max_iteration": 100,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 1,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                    }
                ]
            }
            """)

    def test_cg_spai0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {

                        "solver_type" : "AMGCL",
                        "smoother_type":"spai0",
                        "krylov_type": "cg",
                        "coarsening_type": "ruge_stuben",
                        "max_iteration": 100,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 1,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                    }
                ]
            }
            """)


    def test_amgcl_bicgstabl(self):
        self._RunParametrized("""
            {
                "test_list" : [
                {

                        "solver_type" : "AMGCL",
                        "smoother_type":"iluk",
                        "krylov_type": "bicgstabl",
                        "coarsening_type": "aggregation",
                        "max_iteration": 100,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 0,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                }]
            }
            """)

if __name__ == '__main__':
    KratosUnittest.main()
