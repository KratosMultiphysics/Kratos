from __future__ import print_function, absolute_import, division

import os
import math
import time

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as KratosUtils
dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("MetisApplication")
if dependencies_are_available:
    import KratosMultiphysics.MetisApplication as KratosMetis

# Importing testing utilities
import KratosMultiphysics.testing.utilities as testing_utilities

## DEBUG
#from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestNormalUtilsMPI(KratosUnittest.TestCase):

    def tearDown(self):
        # Clean up temporary files
        KratosUtils.DeleteFileIfExisting(GetFilePath("aux_model_part_with_skin.out.mdpa"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("aux_model_part_with_skin.out.time"))

    @KratosUnittest.skipUnless(dependencies_are_available, "MetisApplication is not available")
    def test_ComputeSimplexNormalModelPartMPI(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        domain_size = 3
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = domain_size

        comm = KratosMultiphysics.DataCommunicator.GetDefault()
        if comm.Rank() == 0:
            self._generate_mdpa_skin(current_model, GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere"))
        while not os.path.exists(GetFilePath("aux_model_part_with_skin.out.mdpa")):
            time.sleep(0.1)
        testing_utilities.ReadModelPart(GetFilePath("aux_model_part_with_skin.out"), model_part)

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(model_part, domain_size)

        #### DEBUG
        ##self._post_process(model_part, comm.Rank())

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

    def _generate_mdpa_skin(self, current_model, fileName):
        model_part = current_model.CreateModelPart("AuxMain")
        model_part_io = KratosMultiphysics.ModelPartIO(fileName)
        model_part_io.ReadModelPart(model_part)
        detect_skin = KratosMultiphysics.SkinDetectionProcess3D(model_part)
        detect_skin.Execute()
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("aux_model_part_with_skin.out"), KratosMultiphysics.IO.WRITE)
        model_part_io.WriteModelPart(model_part)

    def _post_process(self, model_part, rank = 0):
        gid_output = GiDOutputProcess(model_part,
                                    "gid_output_" + str(rank),
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results" : ["NORMAL","PARTITION_INDEX"]
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

