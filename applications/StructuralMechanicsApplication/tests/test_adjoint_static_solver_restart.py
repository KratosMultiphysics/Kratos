import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_adjoint_static_solver import StructuralMechanicsAdjointStaticSolver

class TestAdjointStaticSolverRestart(KratosUnittest.TestCase):
    """PrepareModelPart() of the adjoint static solver on restarted model parts.

    A restarted adjoint model part (IS_RESTARTED and IS_ADJOINT) is already replaced and must be
    kept as is; a restarted primal model part (IS_RESTARTED only) must still be replaced.
    """
    def _CreateSolver(self, is_adjoint: bool) -> StructuralMechanicsAdjointStaticSolver:
        self.model = KratosMultiphysics.Model()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name"       : "Structure",
            "domain_size"           : 3,
            "model_import_settings" : { "input_type": "use_input_model_part" }
        }""")
        solver = StructuralMechanicsAdjointStaticSolver(self.model, settings)
        model_part = solver.GetComputingModelPart()
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        model_part.CreateNewElement("SmallDisplacementElement3D4N", 1, [1, 2, 3, 4], model_part.GetProperties()[1])
        # the solver constructor resets IS_RESTARTED, so mark the restart afterwards
        model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        if is_adjoint:
            model_part.ProcessInfo[StructuralMechanicsApplication.IS_ADJOINT] = True
        return solver

    @staticmethod
    def _GetElementName(solver: StructuralMechanicsAdjointStaticSolver) -> str:
        return KratosMultiphysics.CompareElementsAndConditionsUtility.GetRegisteredName(solver.GetComputingModelPart().GetElement(1))

    def test_restarted_adjoint_model_part_is_kept(self):
        solver = self._CreateSolver(is_adjoint=True)
        solver.PrepareModelPart()
        self.assertEqual(self._GetElementName(solver), "SmallDisplacementElement3D4N")

    def test_restarted_primal_model_part_is_replaced(self):
        solver = self._CreateSolver(is_adjoint=False)
        solver.PrepareModelPart()
        model_part = solver.GetComputingModelPart()
        self.assertEqual(self._GetElementName(solver), "AdjointFiniteDifferencingSmallDisplacementElement3D4N")
        self.assertTrue(model_part.ProcessInfo[StructuralMechanicsApplication.IS_ADJOINT])

if __name__ == '__main__':
    KratosUnittest.main()
