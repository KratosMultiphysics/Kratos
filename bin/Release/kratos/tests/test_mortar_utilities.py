from __future__ import print_function, absolute_import, division

import os
import math

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMortarUtilities(KratosUnittest.TestCase):

    def test_ComputeNodesMeanNormalModelPart(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere"))
        model_part_io.ReadModelPart(model_part)

        detect_skin = KratosMultiphysics.SkinDetectionProcess3D(model_part)
        detect_skin.Execute()

        KratosMultiphysics.MortarUtilities.ComputeNodesMeanNormalModelPart(model_part, True)

        ## DEBUG
        #self._post_process(model_part)

        for node in model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = []
            norm = math.sqrt(node.X**2+node.Y**2+node.Z**2)
            normal.append(node.X/norm)
            normal.append(node.Y/norm)
            normal.append(node.Z/norm)

            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)

            residual = math.sqrt((solution_normal[0]-normal[0])**2+(solution_normal[1]-normal[1])**2+(solution_normal[2]-normal[2])**2)
            self.assertLess(residual, 0.1)

    def test_InvertNormal(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere"))
        model_part_io.ReadModelPart(model_part)

        detect_skin = KratosMultiphysics.SkinDetectionProcess3D(model_part)
        detect_skin.Execute()

        KratosMultiphysics.MortarUtilities.InvertNormal(model_part.Conditions)
        KratosMultiphysics.MortarUtilities.ComputeNodesMeanNormalModelPart(model_part, True)

        ## DEBUG
        #self._post_process(model_part)

        for node in model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = []
            norm = math.sqrt(node.X**2+node.Y**2+node.Z**2)
            normal.append(-node.X/norm)
            normal.append(-node.Y/norm)
            normal.append(-node.Z/norm)

            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)

            residual = math.sqrt((solution_normal[0]-normal[0])**2+(solution_normal[1]-normal[1])**2+(solution_normal[2]-normal[2])**2)
            self.assertLess(residual, 0.1)

    def _post_process(self, model_part):
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

