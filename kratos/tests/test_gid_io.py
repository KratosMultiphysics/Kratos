from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics import compare_two_files_check_process

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

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
                                                "nodal_nonhistorical_results": ["TEMPERATURE","INERTIA","ELEMENTAL_DISTANCES"],
                                                "nodal_flags_results": ["ISOLATED"],
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

    def __InitialRead(self, current_model):
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)
        return model_part

    def __Check(self,output_file,reference_file):

        ## Settings string in json format
        params = KratosMultiphysics.Parameters("""
            {
                "reference_file_name"   : "",
                "output_file_name"      : ""
            }
        """)

        params["reference_file_name"].SetString(GetFilePath(reference_file))
        params["output_file_name"].SetString(output_file)

        cmp_process = compare_two_files_check_process.CompareTwoFilesCheckProcess(params)

        cmp_process.ExecuteInitialize()
        cmp_process.ExecuteBeforeSolutionLoop()
        cmp_process.ExecuteInitializeSolutionStep()
        cmp_process.ExecuteFinalizeSolutionStep()
        cmp_process.ExecuteBeforeOutputStep()
        cmp_process.ExecuteAfterOutputStep()
        cmp_process.ExecuteFinalize()

    def test_gid_io_all(self):
        current_model = KratosMultiphysics.Model()

        model_part = self.__InitialRead(current_model)

        self.__WriteOutput(model_part,"all_active_out")

        self.__Check("all_active_out_0.post.msh","auxiliar_files_for_python_unittest/reference_files/all_active_ref.ref")

    def test_gid_io_deactivation(self):
        current_model = KratosMultiphysics.Model()

        model_part = self.__InitialRead(current_model)

        model_part.Elements[3].Set(KratosMultiphysics.ACTIVE,False)
        model_part.Elements[1796].Set(KratosMultiphysics.ACTIVE,False)

        model_part.Conditions[1947].Set(KratosMultiphysics.ACTIVE,False)
        model_part.Conditions[1948].Set(KratosMultiphysics.ACTIVE,False)

        self.__WriteOutput(model_part,"deactivated_out")

        self.__Check("deactivated_out_0.post.msh","auxiliar_files_for_python_unittest/reference_files/deactivated_ref.ref")

    def test_gid_io_results(self):
        current_model = KratosMultiphysics.Model()

        model_part = self.__InitialRead(current_model)

        # Initialize flag
        for node in model_part.Nodes:
            node.Set(KratosMultiphysics.ISOLATED, False)

        model_part.Nodes[1].Set(KratosMultiphysics.ISOLATED, True)
        model_part.Nodes[2].Set(KratosMultiphysics.ISOLATED, True)
        model_part.Nodes[973].Set(KratosMultiphysics.ISOLATED, True)
        model_part.Nodes[974].Set(KratosMultiphysics.ISOLATED, True)

        # Initialize value
        for node in model_part.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, 0.0)

        model_part.Nodes[1].SetValue(KratosMultiphysics.TEMPERATURE, 100.0)
        model_part.Nodes[2].SetValue(KratosMultiphysics.TEMPERATURE, 200.0)
        model_part.Nodes[973].SetValue(KratosMultiphysics.TEMPERATURE, 300.0)
        model_part.Nodes[974].SetValue(KratosMultiphysics.TEMPERATURE, 400.0)

        vector = KratosMultiphysics.Vector(3)
        vector[0] = 0.0
        vector[1] = 0.0
        vector[2] = 0.0
        # Initialize value
        for node in model_part.Nodes:
            node.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, vector)

        vector[0] = 1.0
        vector[1] = 1.0
        model_part.Nodes[1].SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, vector)
        model_part.Nodes[2].SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, vector)
        model_part.Nodes[973].SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, vector)
        model_part.Nodes[974].SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, vector)

        matrix = KratosMultiphysics.Matrix(2, 2)
        matrix[0, 0] = 0.0
        matrix[0, 1] = 0.0
        matrix[1, 0] = 0.0
        matrix[1, 1] = 0.0
        # Initialize value
        for node in model_part.Nodes:
            node.SetValue(KratosMultiphysics.INERTIA, matrix)

        matrix[0, 0] = 1.0
        matrix[1, 1] = 1.0
        model_part.Nodes[1].SetValue(KratosMultiphysics.INERTIA, matrix)
        model_part.Nodes[2].SetValue(KratosMultiphysics.INERTIA, matrix)
        model_part.Nodes[973].SetValue(KratosMultiphysics.INERTIA, matrix)
        model_part.Nodes[974].SetValue(KratosMultiphysics.INERTIA, matrix)

        self.__WriteOutput(model_part,"results_out")

        self.__Check("results_out_0.post.res","auxiliar_files_for_python_unittest/reference_files/results_out_ref.ref")

    def test_DoubleFreeError(self):
        current_model = KratosMultiphysics.Model()

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

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("all_active_out_0.post.msh")
        kratos_utils.DeleteFileIfExisting("all_active_out_0.post.res")
        kratos_utils.DeleteFileIfExisting("deactivated_out_0.post.msh")
        kratos_utils.DeleteFileIfExisting("deactivated_out_0.post.res")
        kratos_utils.DeleteFileIfExisting("results_out_0.post.msh")
        kratos_utils.DeleteFileIfExisting("results_out_0.post.res")
        kratos_utils.DeleteFileIfExisting("python_scripts.post.lst")
        kratos_utils.DeleteFileIfExisting("tests.post.lst")


if __name__ == '__main__':
    KratosUnittest.main()
