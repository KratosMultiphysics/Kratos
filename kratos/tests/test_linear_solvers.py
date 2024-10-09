import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import platform

psystem = platform.system()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestLinearSolvers(KratosUnittest.TestCase):

    def _RunParametrized(self, my_params_string ):
        all_settings = KratosMultiphysics.Parameters( my_params_string )

        for i in range(all_settings["test_list"].size()):
            settings = all_settings["test_list"][i]
            self._auxiliary_test_function(settings)

    def _auxiliary_test_function(self, settings, matrix_name="test_files/sparse_matrix_files/A.mm"):
        space = KratosMultiphysics.UblasSparseSpace()

        #read the matrices
        A = KratosMultiphysics.CompressedMatrix()
        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath(matrix_name),A)
        self.assertTrue(file_read, msg="The MatrixFile could not be read")

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
        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(settings)

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

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_tfqmr_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "tfqmr",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "ilu0"
                    },
                    {
                        "solver_type" : "tfqmr",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "diagonal"
                    },
                    {
                        "solver_type" : "tfqmr",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 1000,
                        "preconditioner_type" : "none"
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_cg_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "cg",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "diagonal"
                    },
                    {
                        "solver_type" : "cg",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 1000,
                        "preconditioner_type" : "none"
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_deflated_cg_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "deflated_cg"
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_bicgstab_in_core(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "bicgstab",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "ilu0",
                        "scaling": false
                    },
                    {
                        "solver_type" : "bicgstab",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "ilu0",
                        "scaling": false
                    },
                    {
                        "solver_type" : "bicgstab",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "diagonal"
                    },
                    {
                        "solver_type" : "bicgstab",
                        "tolerance" : 1.0e-6,
                        "max_iteration" : 500,
                        "preconditioner_type" : "none"
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_skyline_lu(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "skyline_lu_factorization",
                        "scaling": false
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_bicgstab_iluk(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {

                        "solver_type" : "amgcl",
                        "smoother_type":"iluk",
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

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_lgmres_iluk(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {

                        "solver_type" : "amgcl",
                        "smoother_type":"iluk",
                        "krylov_type": "lgmres",
                        "coarsening_type": "aggregation",
                        "max_iteration": 300,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 100,
                        "verbosity" : 0 ,
                        "tolerance": 1e-6,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 100
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_amgcl_unpreconditioned(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "bicgstab",
                        "preconditioner_type": "dummy",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "gmres",
                        "preconditioner_type": "dummy",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "preconditioner_type": "dummy",
                        "verbosity" : 1
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_amgcl_no_amg_only_preconditioner(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "smoother_type":"ilu0",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "smoother_type":"ilu0",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1,
                        "block_size" : 2
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "smoother_type":"iluk",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "smoother_type":"spai0",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "smoother_type":"damped_jacobi",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "smoother_type":"gauss_seidel",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    },
                    {
                        "solver_type" : "amgcl",
                        "krylov_type": "lgmres",
                        "smoother_type":"chebyshev",
                        "preconditioner_type": "relaxation",
                        "verbosity" : 1
                    }
                ]
            }
            """)

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_amgcl_bicgstab_ilu0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amgcl",
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

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_amgcl_idr_ilu0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amgcl",
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

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_amgcl_bicgstab_spai0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amgcl",
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

    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_cg_spai0(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {

                        "solver_type" : "amgcl",
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


    @KratosUnittest.skipIf(psystem == 'Windows', "This test is temporaly disabled for windows.")
    def test_amgcl_bicgstabl(self):
        self._RunParametrized("""
            {
                "test_list" : [
                {

                        "solver_type" : "amgcl",
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
