from __future__ import print_function, absolute_import, division

import os
import math

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as KratosUtils
dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("MetisApplication")
if dependencies_are_available:
    import KratosMultiphysics.MetisApplication as KratosMetis

# Importing testing utilities
import KratosMultiphysics.testing.utilities as testing_utilities

# DEBUG
from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestNormalUtilsMPI(KratosUnittest.TestCase):

    @KratosUnittest.skipUnless(dependencies_are_available, "MetisApplication is not available")
    def test_ComputeSimplexNormalModelPartMPI(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        testing_utilities.ReadModelPart(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere"), model_part)

        detect_skin = KratosMultiphysics.SkinDetectionProcess3D(model_part)
        detect_skin.Execute()

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(model_part.Conditions, 3)

        ## DEBUG
        self._post_process(model_part)

        for node in model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = []
            norm = math.sqrt(node.X**2+node.Y**2+node.Z**2)
            normal.append(node.X/norm)
            normal.append(node.Y/norm)
            normal.append(node.Z/norm)

            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            solution_normal_norm = math.sqrt(solution_normal[0]**2+solution_normal[1]**2+solution_normal[2]**2)
            solution_normal /= solution_normal_norm

            residual = math.sqrt((solution_normal[0]-normal[0])**2+(solution_normal[1]-normal[1])**2+(solution_normal[2]-normal[2])**2)
            self.assertLess(residual, 0.15)

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

