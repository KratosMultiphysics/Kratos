from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from compare_two_files_check_process import CompareTwoFilesCheckProcess

import os, math

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def CreateShellNodes(mp,element_name):
    mp.CreateNewNode(1, -0.5, - 0.5,  0.0)
    mp.CreateNewNode(2,  0.5,  -0.5,   0.0)
    mp.CreateNewNode(3,  0.5,  0.5,   0.0)
    mp.CreateNewNode(4, -0.5,  0.5,  0.0)
    mp.CreateNewNode(5,  0.0, 0.0, 0.0)

    if element_name.endswith("4N"): # create aditional nodes needed for quad-setup
        mp.CreateNewNode(6, -0.0, -0.5,   0.0)
        mp.CreateNewNode(7,  0.5,  0.0,  0.00)
        mp.CreateNewNode(8, -0.0,  0.5, 0.00)
        mp.CreateNewNode(9, -0.5, -0.0,   0.0)


def CreateShellElements(mp,element_name):
    if element_name.endswith("4N"): # Quadrilaterals
        mp.CreateNewElement(element_name, 1, [1,6,5,9], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [6,2,7,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [5,7,3,8], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [9,5,8,4], mp.GetProperties()[1])
    else: # Triangles
        mp.CreateNewElement(element_name, 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [4,1,5], mp.GetProperties()[1])


def WriteGiDOutput(model_part):
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(model_part,
        "local_axis_"+model_part.Name,
        KratosMultiphysics.Parameters("""
            {
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode"           : "GiD_PostAscii",
                        "WriteDeformedMeshFlag" : "WriteUndeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "nodal_results"       : [],
                    "gauss_point_results" : ["LOCAL_AXIS_1","LOCAL_AXIS_2","LOCAL_AXIS_3",
                                             "LOCAL_MATERIAL_AXIS_1", "LOCAL_MATERIAL_AXIS_2"]
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


class TestLocalAxisVisualization(KratosUnittest.TestCase):
    def tearDown(self):
        # delete all files leftover from the tests
        output_file_name = "local_axis_" + self.element_name + "_0.post.res"
        msh_file_name = "local_axis_" + self.element_name + "_0.post.msh"
        kratos_utils.DeleteFileIfExisting(msh_file_name)
        kratos_utils.DeleteFileIfExisting(output_file_name) # usually this is deleted by the check process but not if it fails

    def test_ThickQuadShellElement(self):
        self.element_name = "ShellThickElementCorotational3D4N"
        self.__ExecuteShellTest()

    def test_ThinQuadShellElement(self):
        self.element_name = "ShellThinElementCorotational3D4N"
        self.__ExecuteShellTest()

    def test_ThickTriShellElement(self):
        self.element_name = "ShellThickElementCorotational3D3N"
        self.__ExecuteShellTest()

    def test_ThinTriShellElement(self):
        self.element_name = "ShellThinElementCorotational3D3N"
        self.__ExecuteShellTest()

    def __ExecuteShellTest(self):
        model_part = KratosMultiphysics.ModelPart(self.element_name)

        CreateShellNodes(model_part, self.element_name)
        CreateShellElements(model_part, self.element_name)

        for i, elem in enumerate(model_part.Elements):
            angle = i*25*math.pi/180 # radians, every 25 degree
            elem.SetValue(StructuralMechanicsApplication.MATERIAL_ORIENTATION_ANGLE, angle)

        WriteGiDOutput(model_part)
        reference_file_name = "local_axis_" + self.element_name + "_0.post.res.ref"
        reference_file_name = os.path.join("local_axis_visualization_ref_result_files", reference_file_name)
        reference_file_name = GetFilePath(reference_file_name)
        output_file_name = "local_axis_" + self.element_name + "_0.post.res"
        self.__CheckResults(reference_file_name, output_file_name)

    def test_3DBeamElement(self):
        self.element_name = "CrBeamElement3D2N"

    def test_3DLinearBeamElement(self):
        self.element_name = "CrBeamLinearElement3D2N"

    def __CheckResults(self, ref_file_name, out_file_name):
        # check the results
        settings_check_process = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : \"""" + ref_file_name + """\",
            "output_file_name"      : \"""" + out_file_name + """\",
            "comparison_type"       : "post_res_file"
        }
        """)

        check_process = CompareTwoFilesCheckProcess(settings_check_process)

        check_process.ExecuteInitialize()
        check_process.ExecuteBeforeSolutionLoop()
        check_process.ExecuteInitializeSolutionStep()
        check_process.ExecuteFinalizeSolutionStep()
        check_process.ExecuteBeforeOutputStep()
        check_process.ExecuteAfterOutputStep()
        check_process.ExecuteFinalize()



if __name__ == '__main__':
    KratosUnittest.main()
