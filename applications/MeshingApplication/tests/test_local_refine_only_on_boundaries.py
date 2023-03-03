﻿# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestLocalRefineOnlyOnBoundaries(KratosUnittest.TestCase):

    def _ComputeVolume(self,Elements):
        vol = 0.0
        for elem in Elements:
            elem_vol = elem.GetGeometry().Volume()
            self.assertTrue(elem_vol>0.0)#check no inverted elems
            vol+=elem_vol
        return vol

    def test_refine_on_boundary_edges(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("main_model_part")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        # We add a dummy var
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/cube_with_5_faces").ReadModelPart(main_model_part)

        #IMPORTANT! must compute the neighbours first!!
        KratosMultiphysics.FindNodalNeighboursProcess(main_model_part).Execute()

        #Flaging all elems to be refined (but only those touching BCs will be refined by the process)
        for elem in main_model_part.Elements:
            elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)

        #refining
        refiner = MeshingApplication.LocalRefineTetrahedraMeshOnlyOnBoundaries(main_model_part)
        refine_on_reference = False
        interpolate_internal_variables = True
        refiner.LocalRefineMesh( refine_on_reference, interpolate_internal_variables)

        #checking
        refined_vol = self._ComputeVolume(main_model_part.Elements)
        self.assertAlmostEqual(refined_vol, 1.0, 12)
        self.assertEqual(main_model_part.NumberOfElements(), 604)
        self.assertEqual(current_model["main_model_part.MainPart.VolumeParts.Body1"].NumberOfConditions(), 0)
        self.assertEqual(current_model["main_model_part.MainPart.VolumeParts.Body1"].NumberOfElements(), 604)
        self.assertEqual(main_model_part.NumberOfNodes(), 224)
        self.assertEqual(main_model_part.NumberOfConditions(), 360)
        self.assertEqual(current_model["main_model_part.MainPart.Wall_BC.Wall1"].NumberOfConditions(), 360)
        self.assertEqual(current_model["main_model_part.MainPart.Wall_BC.Wall1"].NumberOfElements(), 0)

if __name__ == '__main__':
    KratosUnittest.main()