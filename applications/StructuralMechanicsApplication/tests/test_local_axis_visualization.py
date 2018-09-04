from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from postprocess_eigenvalues_process import PostProcessEigenvaluesProcess
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


def CreateGiDOutput(model_part):
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(model_part,
        "gid_output_"+model_part.Name,
        KratosMultiphysics.Parameters("""
            {
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostAscii",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "nodal_results"       : [],
                    "gauss_point_results" : ["LOCAL_AXIS_1","LOCAL_AXIS_2","LOCAL_AXIS_3", "LOCAL_MATERIAL_AXIS_1", "LOCAL_MATERIAL_AXIS_2"]
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
        kratos_utils.DeleteFileIfExisting("Structure_EigenResults_0.post.msh")
        kratos_utils.DeleteFileIfExisting("Structure_EigenResults_0.post.res") # usually this is deleted by the check process but not if it fails

    def test_ThickQuadShellLocalAxis(self):
        element_name = "ShellThickElementCorotational3D4N"
        test_model = KratosMultiphysics.Model()
        model_part = KratosMultiphysics.ModelPart(element_name)
        test_model.AddModelPart(model_part)

        CreateShellNodes(model_part, element_name)
        CreateShellElements(model_part, element_name)


        for i, elem in enumerate(model_part.Elements):
            elem.SetValue(StructuralMechanicsApplication.MATERIAL_ORIENTATION_ANGLE, (i+1)/18*math.pi) # every 10 degrees


        CreateGiDOutput(model_part)






    # def test_PostprocessEigenvaluesProcess(self):
    #     model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    #     model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    #     model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
    #     model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

    #     CreateNodes(model_part)

    #     # adding dofs is needed for the process internally
    #     KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
    #     KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
    #     KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)
    #     KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,model_part)
    #     KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,model_part)
    #     KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,model_part)

    #     test_model.AddModelPart(model_part)

    #     # set EigenValues and -Vectors
    #     num_eigenvalues = 4
    #     eigenval_vector = GetEigenValueVector(num_eigenvalues)
    #     model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR] = eigenval_vector

    #     for node in model_part.Nodes:
    #         node.SetValue(StructuralMechanicsApplication.EIGENVECTOR_MATRIX,
    #                       GetEigenVectorMatrix(num_eigenvalues, node.Id))

    #     # Use the process
    #     # here the minimum settings are specified to test the default values!
    #     settings_eigen_process = KratosMultiphysics.Parameters("""{"result_file_format_use_ascii" : true}""")

    #     post_eigen_process = PostProcessEigenvaluesProcess(test_model, settings_eigen_process)

    #     post_eigen_process.ExecuteInitialize()
    #     post_eigen_process.ExecuteBeforeSolutionLoop()
    #     post_eigen_process.ExecuteInitializeSolutionStep()
    #     post_eigen_process.ExecuteFinalizeSolutionStep()
    #     post_eigen_process.ExecuteBeforeOutputStep()
    #     post_eigen_process.ExecuteAfterOutputStep()
    #     post_eigen_process.ExecuteFinalize()

    #     # check the results
    #     settings_check_process = KratosMultiphysics.Parameters("""
    #     {
    #         "reference_file_name"   : "",
    #         "output_file_name"      : "",
    #         "remove_output_file"    : true,
    #         "comparison_type"       : "post_res_file"
    #     }
    #     """)

    #     settings_check_process["reference_file_name"].SetString(GetFilePath("test_postprocess_eigenvalues_process.ref"))
    #     settings_check_process["output_file_name"].SetString("Structure_EigenResults_0.post.res")

    #     check_process = CompareTwoFilesCheckProcess(test_model, settings_check_process)

    #     check_process.ExecuteInitialize()
    #     check_process.ExecuteBeforeSolutionLoop()
    #     check_process.ExecuteInitializeSolutionStep()
    #     check_process.ExecuteFinalizeSolutionStep()
    #     check_process.ExecuteBeforeOutputStep()
    #     check_process.ExecuteAfterOutputStep()
    #     check_process.ExecuteFinalize()


if __name__ == '__main__':
    KratosUnittest.main()
