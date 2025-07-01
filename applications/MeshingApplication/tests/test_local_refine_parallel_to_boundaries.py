# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestLocalRefineParallelToBoundaries(KratosUnittest.TestCase):

    def _ComputeVolume(self,Elements):
        vol = 0.0
        for elem in Elements:
            elem_vol = elem.GetGeometry().Volume()
            self.assertTrue(elem_vol>0.0)#check no inverted elems
            vol+=elem_vol
        return vol

    def test_refine_boundary_elems(self):
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
        KratosMultiphysics.FindNodalNeighboursProcess(main_model_part).Execute()

        #Flaging all elems to be refined (but only those touching BCs will be refined by the process)
        for elem in main_model_part.Elements:
            elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)

        #refining
        refiner = MeshingApplication.LocalRefineTetrahedraMeshParallelToBoundaries(main_model_part)
        refine_on_reference = False
        interpolate_internal_variables = True
        refiner.LocalRefineMesh( refine_on_reference, interpolate_internal_variables)

        #checking
        refined_vol = self._ComputeVolume(main_model_part.Elements)
        self.assertAlmostEqual(refined_vol, 1.0, 12)
        self.assertEqual(main_model_part.NumberOfElements(), 701)
        self.assertEqual(main_model_part.NumberOfNodes(), 169)
        self.assertEqual(main_model_part.NumberOfConditions(), 90)


        # export to gid commented, uncomment in case you want to see how the refined mesh looks like
        #the initial geometry was a cube with 3X3x3 tets, identical size
        #one of the sides of the cube does not have BCs, so that side will be refined.
        #in the remaining 5 faces, where triangles are presents, their edges are not cut, so the original number of BCs is kept

        # gid_output = GiDOutputProcess(main_model_part,
        #                             "gid_output",
        #                             KratosMultiphysics.Parameters("""
        #                                 {
        #                                     "result_file_configuration" : {
        #                                         "gidpost_flags": {
        #                                             "GiDPostMode": "GiD_PostBinary",
        #                                             "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                             "WriteConditionsFlag": "WriteConditions",
        #                                             "MultiFileFlag": "SingleFile"
        #                                         },
        #                                         "nodal_results"       : ["DISTANCE"]
        #                                     }
        #                                 }
        #                                 """)
        #                             )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()



if __name__ == '__main__':
    KratosUnittest.main()
