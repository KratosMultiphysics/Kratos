from __future__ import print_function, absolute_import, division
import KratosMultiphysics
from KratosMultiphysics.mpi import *
#import KratosMultiphysics.mpi
import KratosMultiphysics.TrilinosApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
import os

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName


class TestLinearSolvers(KratosUnittest.TestCase):

    def _RunParametrized(self, my_params_string ):
        all_settings = KratosMultiphysics.Parameters( my_params_string )

        for i in range(all_settings["test_list"].size()):
            settings = all_settings["test_list"][i]
            self._auxiliary_test_function(settings)

    def _auxiliary_test_function(self, settings, matrix_name="A.mm", absolute_norm=False):
        comm = KratosMultiphysics.TrilinosApplication.CreateCommunicator()
        space = KratosMultiphysics.TrilinosApplication.TrilinosSparseSpace()

        #read the matrices
        pA = space.ReadMatrixMarketMatrix(GetFilePath(matrix_name),comm)
        n = space.Size1(pA.GetReference())

        pAoriginal = space.ReadMatrixMarketMatrix(GetFilePath(matrix_name),comm)
        pb  = space.CreateEmptyVectorPointer(comm)
        space.ResizeVector(pb,n)
        space.SetToZeroVector(pb.GetReference())
        space.SetValue(pb.GetReference(), 0, 1.0)  #pb[1] = 1.0

        px = space.CreateEmptyVectorPointer(comm)
        space.ResizeVector(px,n)
        pboriginal = space.CreateEmptyVectorPointer(comm) #create a copy of b
        space.ResizeVector(pboriginal,n)
        space.SetToZeroVector(pboriginal.GetReference())
        space.UnaliasedAdd(pboriginal.GetReference(), 1.0,pb.GetReference())

        space.SetToZeroVector(px.GetReference())
        #space.SetToZeroVector(boriginal)
        #space.UnaliasedAdd(boriginal, 1.0, b) #boriginal=1*bs


        #construct the solver
        linear_solver = trilinos_linear_solver_factory.ConstructSolver(settings)

        #solve
        linear_solver.Solve(pA.GetReference(),px.GetReference(),pb.GetReference())

        #test the results
        ptmp = space.CreateEmptyVectorPointer(comm)
        space.ResizeVector(ptmp,n)
        space.SetToZeroVector(ptmp.GetReference())
        space.Mult(pAoriginal.GetReference(),px.GetReference(),ptmp.GetReference())


        pcheck = space.CreateEmptyVectorPointer(comm)
        space.ResizeVector(pcheck,n)
        space.SetToZeroVector(pcheck.GetReference())
        space.UnaliasedAdd(pcheck.GetReference(), 1.0,pboriginal.GetReference())
        space.UnaliasedAdd(pcheck.GetReference(), -1.0,ptmp.GetReference())


        achieved_norm = space.TwoNorm(pcheck.GetReference())

        tolerance = 1e-9
        if(settings.Has("tolerance")):
            tolerance = settings["tolerance"].GetDouble()

        nproc = KratosMultiphysics.DataCommunicator.GetDefault().Size()
        target_norm = tolerance*space.TwoNorm(pboriginal.GetReference())*nproc #multiplying by nproc the target tolerance to give some slack. Not really nice :-(

        if(achieved_norm > target_norm):
            print("target_norm   :",target_norm)
            print("achieved_norm :",achieved_norm)
        self.assertTrue(achieved_norm <= target_norm)

        #destroy the preconditioner - this is needed since  the solver should be destroyed before the destructor of the system matrix is called
        del linear_solver

    def test_amesos_superludist(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Superludist") ):
            self.skipTest("Amesos_Superludist is not among the available Amesos Solvers")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos",
                        "amesos_solver_type" : "Amesos_Superludist"
                    }
                ]
            }
            """)

    def test_amesos_mumps(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Mumps") ):
            self.skipTest("Amesos_Mumps is not among the available Amesos Solvers")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos",
                        "amesos_solver_type" : "Amesos_Mumps"
                    }
                ]
            }
            """)


    def test_amesos_klu(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Klu") ):
            self.skipTest("Amesos_Klu is not among the available Amesos Solvers")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos",
                        "amesos_solver_type" : "Amesos_Klu"
                    }
                ]
            }
            """)

    def test_amesos_superludist_2(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Superludist") ):
            self.skipTest("Amesos_Superludist is not among the available Amesos Solvers")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "super_lu_dist"
                    }
                ]
            }
            """)

    def test_amesos_mumps_2(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Mumps") ):
            self.skipTest("Amesos_Mumps is not among the available Amesos Solvers")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "mumps"
                    }
                ]
            }
            """)


    def test_amesos_klu_2(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Klu") ):
            self.skipTest("Amesos_Klu is not among the available Amesos Solvers")

        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type" : "klu"
                    }
                ]
            }
            """)

    def test_aztec_cg(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type": "cg",
                        "tolerance" : 1.0e-9,
                        "max_iteration" : 200,
                        "preconditioner_type" : "None",
                        "overlap_level":0,
                        "scaling":false,
                        "verbosity":0,
                        "trilinos_aztec_parameter_list": {},
                        "trilinos_preconditioner_parameter_list": {}
                    },
                    {
                        "solver_type": "cg",
                        "tolerance" : 1.0e-9,
                        "max_iteration" : 200,
                        "preconditioner_type" : "None",
                        "overlap_level":0,
                        "verbosity":0,
                        "scaling":true,
                        "trilinos_aztec_parameter_list": {},
                        "trilinos_preconditioner_parameter_list": {}
                    },
                    {
                        "solver_type": "cg",
                        "tolerance" : 1.0e-9,
                        "max_iteration" : 200,
                        "preconditioner_type" : "None",
                        "overlap_level":0,
                        "verbosity":0,
                        "scaling":true,
                        "trilinos_aztec_parameter_list": {},
                        "trilinos_preconditioner_parameter_list": {}
                    }
                ]
            }
            """)

    def test_aztec_gmres(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                        "solver_type": "gmres",
                        "tolerance" : 1.0e-9,
                        "max_iteration" : 200,
                        "preconditioner_type" : "ILU0",
                        "overlap_level":1,
                        "verbosity":0,
                        "scaling":true,
                        "trilinos_aztec_parameter_list": {},
                        "trilinos_preconditioner_parameter_list": {}
                    },
                    {
                        "solver_type": "gmres",
                        "tolerance" : 1.0e-9,
                        "max_iteration" : 200,
                        "preconditioner_type" : "ILUT",
                        "overlap_level":1,
                        "verbosity":0,
                        "scaling":true,
                        "trilinos_aztec_parameter_list": {},
                        "trilinos_preconditioner_parameter_list": {}
                    },
                    {
                        "solver_type": "gmres",
                        "tolerance" : 1.0e-9,
                        "max_iteration" : 200,
                        "preconditioner_type" : "ILUT",
                        "overlap_level":0,
                        "verbosity":0,
                        "scaling":true,
                        "trilinos_aztec_parameter_list": {},
                        "trilinos_preconditioner_parameter_list": {}
                    }
                ]
            }
            """)

    def test_ml_symmetric(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "multi_level",
                    "tolerance" : 1.0e-6,
                    "max_iteration" : 200,
                    "max_levels" : 3,
                    "scaling":false,
                    "reform_preconditioner_at_each_step":true,
                    "symmetric":true,
                    "verbosity":0,
                    "trilinos_aztec_parameter_list": {},
                    "trilinos_ml_parameter_list": {}
                    }
                ]
            }
            """)

    def test_ml_symmetric_scaling(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "multi_level",
                    "tolerance" : 1.0e-6,
                    "max_iteration" : 200,
                    "max_levels" : 3,
                    "scaling":true,
                    "reform_preconditioner_at_each_step":true,
                    "symmetric":true,
                    "verbosity":0,
                    "trilinos_aztec_parameter_list": {},
                    "trilinos_ml_parameter_list": {}
                    }
                ]
            }
            """)

    def test_ml_nonsymmetric(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "multi_level",
                    "tolerance" : 1.0e-6,
                    "max_iteration" : 200,
                    "max_levels" : 3,
                    "scaling":false,
                    "reform_preconditioner_at_each_step":true,
                    "symmetric":false,
                    "verbosity":0,
                    "trilinos_aztec_parameter_list": {},
                    "trilinos_ml_parameter_list": {}
                    }
                ]
            }
            """)

    def test_ml_nonsymmetric_scaling(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "multi_level",
                    "tolerance" : 1.0e-6,
                    "max_iteration" : 200,
                    "max_levels" : 3,
                    "scaling":true,
                    "reform_preconditioner_at_each_step":true,
                    "symmetric":false,
                    "verbosity":0,
                    "trilinos_aztec_parameter_list": {},
                    "trilinos_ml_parameter_list": {}
                    }
                ]
            }
            """)

    def test_amgcl_mpi_solver_cg(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "amgcl",
                    "tolerance":1.0e-9,
                    "scaling":false,
                    "krylov_type":"cg",
                    "verbosity":0
                    }
                ]
            }
            """)

    def test_amgcl_mpi_solver_bicgstab(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "amgcl",
                    "scaling":false,
                    "tolerance":1.0e-9,
                    "krylov_type":"bicgstab",
                    "verbosity":0
                    }
                ]
            }
            """)

    def test_amgcl_mpi_solver_bicgstabl(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "amgcl",
                    "tolerance":1.0e-9,
                    "scaling":false,
                    "krylov_type":"bicgstabl",
                    "verbosity":0
                    }
                ]
            }
            """)

    def test_amgcl_mpi_solver_gmres(self):
        self._RunParametrized("""
            {
                "test_list" : [
                    {
                    "solver_type": "amgcl",
                    "tolerance":1.0e-9,
                    "scaling":false,
                    "krylov_type":"gmres",
                    "verbosity":0
                    }
                ]
            }
            """)


if __name__ == '__main__':
    KratosUnittest.main()
