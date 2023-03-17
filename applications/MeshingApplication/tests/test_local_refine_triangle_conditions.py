# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestLocalRefineTriangleMeshConditions(KratosUnittest.TestCase):

    def _ComputeArea(self,Conditions):
        total_area = 0.0
        for cond in Conditions:
            cond_area = cond.GetGeometry().Area()
            self.assertTrue(cond_area>0.0)#check no inverted elems
            total_area+=cond_area
        return total_area

    def test_refine_condition_mesh(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        # We add a dummy var
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/cube_with_5_faces").ReadModelPart(main_model_part)

        #IMPORTANT! must compute the neighbours first!!
        KratosMultiphysics.FindGlobalNodalNeighboursForConditionsProcess(main_model_part).Execute()

        #Flaging all elems to be refined (but only those touching BCs will be refined by the process)
        for cond in main_model_part.Conditions:
            cond.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)

        #refining
        refiner = MeshingApplication.LocalRefineTriangleMeshConditions(main_model_part)
        refine_on_reference = False
        interpolate_internal_variables = True
        refiner.LocalRefineMesh( refine_on_reference, interpolate_internal_variables)

        #checking
        refined_area = self._ComputeArea(main_model_part.Conditions)
        self.assertAlmostEqual(refined_area, 5.0, 12) #5 faces of a 1x1x1 cube
        self.assertEqual(main_model_part.NumberOfElements(), 162)
        self.assertEqual(main_model_part.NumberOfNodes(), 205)
        self.assertEqual(main_model_part.NumberOfConditions(), 360)


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
