
from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from new_linear_solver_factory import ConstructSolver

class TestEigenDirectSolver(KratosUnittest.TestCase):
    def _test_eigen_direct_solver(self, class_name, solver_type):
        # check if application is available
        try:
            import KratosMultiphysics.EigenSolversApplication
        except:
            self.skipTest("EigenSolversApplication is not available")
        
        # check if solver is available
        if (not hasattr(KratosMultiphysics.EigenSolversApplication, class_name)):
            self.skipTest(class_name + " is not included in the compilation of the EigenSolversApplication")


        space = KratosMultiphysics.UblasSparseSpace()


        settings = KratosMultiphysics.Parameters('{ "solver_type" : "' + solver_type + '" }')
        
        solver = ConstructSolver(settings)


        a = KratosMultiphysics.CompressedMatrix()
        
        KratosMultiphysics.ReadMatrixMarketMatrix('../../../kratos/tests/A.mm', a) # symmetric test matrix
        
        dimension = a.Size1()

        b_exp = KratosMultiphysics.Vector(dimension) # [1, 2, ..., dimension-1, dimension]

        for i in range(dimension):
            b_exp[i] = i + 1

        x = KratosMultiphysics.Vector(dimension)
        

        solver.Solve(a, x, b_exp)


        b_act = KratosMultiphysics.Vector(dimension)
        space.Mult(a, x, b_act)

        for i in range(dimension):
            self.assertLess(abs(b_act[i] - b_exp[i]), 1e-7)

    def test_eigen_sparse_lu(self):
        self._test_eigen_direct_solver('SparseLUSolver', 'eigen_sparse_lu')
    
    def test_eigen_pardiso_lu(self):
        self._test_eigen_direct_solver('PardisoLUSolver', 'eigen_pardiso_lu')

    def test_eigen_pardiso_ldlt(self):
        self._test_eigen_direct_solver('PardisoLDLTSolver', 'eigen_pardiso_ldlt')

    def test_eigen_pardiso_llt(self):
        self._test_eigen_direct_solver('PardisoLLTSolver', 'eigen_pardiso_llt')

if __name__ == '__main__':
    KratosUnittest.main()