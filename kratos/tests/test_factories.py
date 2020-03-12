from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestFactories(KratosUnittest.TestCase):

    def _auxiliary_test_function_BuilderAndSolver(self, settings, name):
        #builder_and_solver = KM.BuilderAndSolverFactory().Create(settings)
        self.assertTrue(KM.BuilderAndSolverFactory().Has(settings["name"].GetString()))
        #self.assertEqual(builder_and_solver.Info(), name)
        
    def test_ResidualBasedEliminationBuilderAndSolver(self):
        settings = KM.Parameters("""
        {
            "name" : "elimination_builder_and_solver"
        }
        """)
        self._auxiliary_test_function_BuilderAndSolver(settings, "ResidualBasedEliminationBuilderAndSolver")

    def test_ResidualBasedEliminationBuilderAndSolverWithConstraints(self):
        settings = KM.Parameters("""
        {
            "name" : "elimination_builder_and_solver_with_constraints"
        }
        """)
        self._auxiliary_test_function_BuilderAndSolver(settings, "ResidualBasedEliminationBuilderAndSolverWithConstraints")
        
    def test_ResidualBasedBlockBuilderAndSolver(self):
        settings = KM.Parameters("""
        {
            "name" : "block_builder_and_solver"
        }
        """)
        self._auxiliary_test_function_BuilderAndSolver(settings, "ResidualBasedBlockBuilderAndSolver")

if __name__ == '__main__':
    KratosUnittest.main()
