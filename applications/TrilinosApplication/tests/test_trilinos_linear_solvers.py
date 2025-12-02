import KratosMultiphysics
import KratosMultiphysics.TrilinosApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
from KratosMultiphysics.mpi import DataCommunicatorFactory

import pathlib

def GetFilePath(fileName):
    return str(pathlib.Path(__file__).absolute().parent / fileName)

class TestLinearSolvers(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        default_data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        cls.data_comm_name = "AllExceptFirst"
        ranks_sub_comm = [i for i in range(1,default_data_comm.Size())]
        cls.sub_comm = DataCommunicatorFactory.CreateFromRanksAndRegister(default_data_comm, ranks_sub_comm, cls.data_comm_name)

    @classmethod
    def tearDownClass(cls):
        KratosMultiphysics.ParallelEnvironment.UnregisterDataCommunicator(cls.data_comm_name)

    def _RunParametrized(self, my_params_string ):
        all_settings = KratosMultiphysics.Parameters( my_params_string )

        for i in range(all_settings["test_list"].size()):
            settings = all_settings["test_list"][i]
            self._auxiliary_test_function(settings)

    def _RunParametrizedWithSubComm(self, my_params_string ):
        all_settings = KratosMultiphysics.Parameters( my_params_string )

        # create a sub communicator, containing all ranks except rank 0
        default_data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        size = default_data_comm.Size()
        if size < 3:
            self.skipTest("This test needs at least 3 mpi processes")

        if default_data_comm.Rank() == 0:
            self.assertFalse(self.sub_comm.IsDefinedOnThisRank())
        else:
            self.assertTrue(self.sub_comm.IsDefinedOnThisRank())

        if not self.sub_comm.IsDefinedOnThisRank():
            # this rank does not participate
            return

        for i in range(all_settings["test_list"].size()):
            settings = all_settings["test_list"][i]
            self._auxiliary_test_function(settings, self.sub_comm)

    def _auxiliary_test_function(self, settings, data_comm=KratosMultiphysics.Testing.GetDefaultDataCommunicator(), matrix_name="A.mm", absolute_norm=False):
        comm = KratosMultiphysics.TrilinosApplication.CreateEpetraCommunicator(data_comm)
        space = KratosMultiphysics.TrilinosApplication.TrilinosSparseSpace()

        #read the matrices
        pA = space.ReadMatrixMarketMatrix(GetFilePath( "auxiliary_files/matrix_market_files/" + matrix_name),comm)
        n = space.Size1(pA.GetReference())

        pAoriginal = space.ReadMatrixMarketMatrix(GetFilePath( "auxiliary_files/matrix_market_files/" + matrix_name),comm)
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

        nproc = data_comm.Size()
        target_norm = tolerance*space.TwoNorm(pboriginal.GetReference())*nproc #multiplying by nproc the target tolerance to give some slack. Not really nice :-(

        self.assertLessEqual(achieved_norm, target_norm)

        #destroy the preconditioner - this is needed since  the solver should be destroyed before the destructor of the system matrix is called
        del linear_solver

@KratosUnittest.skipUnless(hasattr(KratosMultiphysics.TrilinosApplication, 'AmesosSolver'), "AmesosSolver was explicitly disabled.")
class TestAmesosLinearSolvers(TestLinearSolvers):
    def test_amesos_superludist(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Superludist") ):
            self.skipTest("Amesos_Superludist is not among the available Amesos Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos",
                        "amesos_solver_type" : "Amesos_Superludist"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos_mumps(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Mumps") ):
            self.skipTest("Amesos_Mumps is not among the available Amesos Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos",
                        "amesos_solver_type" : "Amesos_Mumps"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos_klu(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Klu") ):
            self.skipTest("Amesos_Klu is not among the available Amesos Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos",
                        "amesos_solver_type" : "Amesos_Klu"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos_superludist_2(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Superludist") ):
            self.skipTest("Amesos_Superludist is not among the available Amesos Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "super_lu_dist"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos_mumps_2(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Mumps") ):
            self.skipTest("Amesos_Mumps is not among the available Amesos Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "mumps"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos_klu_2(self):
        if( not KratosMultiphysics.TrilinosApplication.AmesosSolver.HasSolver("Amesos_Klu") ):
            self.skipTest("Amesos_Klu is not among the available Amesos Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "klu"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

@KratosUnittest.skipUnless(hasattr(KratosMultiphysics.TrilinosApplication, 'Amesos2Solver'), "Amesos2Solver is not included in the compilation.")
class TestAmesos2LinearSolvers(TestLinearSolvers):
    def test_amesos2_superludist(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("amesos2_superludist") ):
            self.skipTest("amesos2_superludist is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos2",
                        "amesos2_solver_type" : "amesos2_superludist"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        # NOTE: SuperLUDist apparently is not compatible with subcommunicators
        # with self.subTest('SubComm'):
        #     self._RunParametrizedWithSubComm(params_string)

    def test_amesos2_mumps(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("amesos2_mumps") ):
            self.skipTest("amesos2_mumps is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos2",
                        "amesos2_solver_type" : "amesos2_mumps"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos2_klu(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("amesos2_klu2") ):
            self.skipTest("amesos2_klu2 is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos2",
                        "amesos2_solver_type" : "amesos2_klu2"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos2_basker(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("basker") ):
            self.skipTest("basker is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "amesos2",
                        "amesos2_solver_type" : "basker"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos2_superludist_2(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("amesos2_superludist") ):
            self.skipTest("amesos2_superludist is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "super_lu_dist2"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos2_mumps_2(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("amesos2_mumps") ):
            self.skipTest("amesos2_mumps is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "mumps2"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos2_klu_2(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("amesos2_klu2") ):
            self.skipTest("amesos2_klu2 is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "klu2"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amesos2_basker_2(self):
        if( not KratosMultiphysics.TrilinosApplication.Amesos2Solver.HasSolver("basker") ):
            self.skipTest("basker is not among the available Amesos2 Solvers")

        params_string = """
            {
                "test_list" : [
                    {
                        "solver_type" : "basker"
                    }
                ]
            }
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

@KratosUnittest.skipUnless(hasattr(KratosMultiphysics.TrilinosApplication, 'AztecSolver'), "AztecSolver was explicitly disabled.")
class TestAztecLinearSolvers(TestLinearSolvers):
    def test_aztec_cg(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_aztec_gmres(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)


@KratosUnittest.skipUnless(hasattr(KratosMultiphysics.TrilinosApplication, 'MultiLevelSolver'), "MultiLevelSolver was explicitly disabled.")
class TestMLLinearSolvers(TestLinearSolvers):
    def test_ml_symmetric(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_ml_symmetric_scaling(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_ml_nonsymmetric(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_ml_nonsymmetric_scaling(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

class TestAMGCLMPILinearSolvers(TestLinearSolvers):
    def test_amgcl_mpi_solver_cg(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amgcl_mpi_solver_bicgstab(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amgcl_mpi_solver_bicgstabl(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)

    def test_amgcl_mpi_solver_gmres(self):
        params_string = """
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
            """

        with self.subTest('All ranks (MPI_COMM_WORLD)'):
            self._RunParametrized(params_string)

        with self.subTest('SubComm'):
            self._RunParametrizedWithSubComm(params_string)


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
