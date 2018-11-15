from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestSkinDetectionProcess(KratosUnittest.TestCase):

    def test_SkinDetectionProcess(self):
        current_model = KratosMultiphysics.Model()

        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        model_part = current_model.CreateModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/coarse_sphere"))
        model_part_io.ReadModelPart(model_part)

        # We set a flag in the already knon node in the skin
        for node in model_part.GetSubModelPart("Skin_Part").Nodes:
            node.Set(KratosMultiphysics.ACTIVE, True)

        detect_skin = KratosMultiphysics.SkinDetectionProcess3D(model_part)
        detect_skin.Execute()

        ## DEBUG
        #self._post_process(model_part)

        for node in model_part.Nodes:
            self.assertEqual(node.Is(KratosMultiphysics.INTERFACE), node.Is(KratosMultiphysics.ACTIVE))

    def test_SkinDetectionProcessWithAssign(self):
        current_model = KratosMultiphysics.Model()
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        model_part = current_model.CreateModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/coarse_sphere"))
        model_part_io.ReadModelPart(model_part)

        skin_detection_parameters = KratosMultiphysics.Parameters("""
        {
            "list_model_parts_to_assign_conditions" : ["Skin_Part"]
        }
        """)

        detect_skin = KratosMultiphysics.SkinDetectionProcess3D(model_part, skin_detection_parameters)
        detect_skin.Execute()

        ## DEBUG
        #self._post_process(model_part)

        self.assertEqual(model_part.GetSubModelPart("Skin_Part").NumberOfConditions(), model_part.NumberOfConditions())

    def _post_process(self, model_part):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_flags_results" : ["INTERFACE","ACTIVE"]
                                            }
                                        }
                                        """)
                                    )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()
