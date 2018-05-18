from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMortarUtilities(KratosUnittest.TestCase):

    def test_ComputeNodesMeanNormalModelPart(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        model = KratosMultiphysics.Model()
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("coarse_sphere"))
        model_part_io.ReadModelPart(model_part)
        model.AddModelPart(model_part)

        detect_skin = KratosMultiphysics.SkinDetectionProcess3D(model_part)
        detect_skin.Execute()

        normal_compute = KratosMultiphysics.MortarUtilities
        normal_compute.ComputeNodesMeanNormalModelPart(model_part.GetSubModelPart("SkinModelPart"))

        ## DEBUG
        #self._post_process(model_part)

        import from_json_check_result_process

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : ["NORMAL"],
            "input_file_name"      : "mortar_mapper_python_tests/normal_results.json",
            "model_part_name"      : "Main",
            "sub_model_part_name"  : "Skin_Part"
        }
        """)

        check = from_json_check_result_process.FromJsonCheckResultProcess(model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        ## The following is used to create the solution database
        #import json_output_process

        #out_parameters = KratosMultiphysics.Parameters("""
        #{
            #"output_variables"     : ["NORMAL"],
            #"output_file_name"     : "mortar_mapper_python_tests/normal_results.json",
            #"model_part_name"      : "Main",
            #"sub_model_part_name"  : "Skin_Part"
        #}
        #""")

        #out = json_output_process.JsonOutputProcess(model, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

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
                                                "nodal_results" : ["NORMAL"]
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

