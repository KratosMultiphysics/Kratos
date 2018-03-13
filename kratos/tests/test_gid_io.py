from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from gid_output_process import GiDOutputProcess
import os

def GetFilePath(fileName):

    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestGidIO(KratosUnittest.TestCase):

    def __WriteOutput(self, model_part, output_file):

        gid_output = GiDOutputProcess(model_part,
                                    output_file,
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration": {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostAscii",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "file_label": "time",
                                                "output_control_type": "step",
                                                "output_frequency": 1.0,
                                                "body_output": true,
                                                "node_output": false,
                                                "skin_output": false,
                                                "plane_output": [],
                                                "nodal_results": ["DISPLACEMENT","VISCOSITY"],
                                                "nodal_nonhistorical_results": [],
                                                "nodal_flags_results": [],
                                                "gauss_point_results": [],
                                                "additional_list_files": []
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

    def __InitialRead(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)
        return model_part

    def __Check(self,output_file,reference_file):
        import compare_two_files_check_process

        ## Settings string in json format
        params = KratosMultiphysics.Parameters("""
            {
                "file_name_1"            : "",
                "file_name_2"            : "",
                "deterministic"          : true
            }
        """)

        params["file_name_1"].SetString(GetFilePath(output_file))
        params["file_name_2"].SetString(GetFilePath(reference_file))

        cmp_process = compare_two_files_check_process.CompareTwoFilesCheckProcess(KratosMultiphysics.ModelPart(), params)

        cmp_process.ExecuteInitialize()
        cmp_process.ExecuteBeforeSolutionLoop()
        cmp_process.ExecuteInitializeSolutionStep()
        cmp_process.ExecuteFinalizeSolutionStep()
        cmp_process.ExecuteBeforeOutputStep()
        cmp_process.ExecuteAfterOutputStep()
        cmp_process.ExecuteFinalize()

    def test_gid_io_all(self):
        model_part = self.__InitialRead()

        self.__WriteOutput(model_part,"all_active_out")

        self.__Check("all_active_out_0.post.msh","all_active_ref.ref")

    def test_gid_io_deactivation(self):
        model_part = self.__InitialRead()

        model_part.Elements[3].Set(KratosMultiphysics.ACTIVE,False)
        model_part.Elements[1796].Set(KratosMultiphysics.ACTIVE,False)

        model_part.Conditions[1947].Set(KratosMultiphysics.ACTIVE,False)
        model_part.Conditions[1948].Set(KratosMultiphysics.ACTIVE,False)

        self.__WriteOutput(model_part,"deactivated_out")

        self.__Check("deactivated_out_0.post.msh","deactivated_ref.ref")

    def test_DoubleFreeError(self):

        output_file_1 = "outFile"
        output_file_2 = "otherFile"

        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostAscii
        multifile = KratosMultiphysics.MultiFileFlag.MultipleFiles
        deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

        gid_io_1 = KratosMultiphysics.GidIO(output_file_1, gid_mode, multifile,
                         deformed_mesh_flag, write_conditions)

        gid_io_2 = KratosMultiphysics.GidIO(output_file_2, gid_mode, multifile,
                         deformed_mesh_flag, write_conditions)

        gid_io_1 = None
        gid_io_2 = None



if __name__ == '__main__':
    KratosUnittest.main()
