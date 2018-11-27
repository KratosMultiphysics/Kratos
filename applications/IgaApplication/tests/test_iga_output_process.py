from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.IgaApplication as KratosIga
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from iga_output_process import IgaOutputProcess
from compare_two_files_check_process import CompareTwoFilesCheckProcess

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def CreateNodes(model_part):
    model_part.CreateNewNode(1, -0.5, - 0.45,  0.1)
    model_part.CreateNewNode(2,  0.7,  -0.5,   0.2)
    model_part.CreateNewNode(3,  0.55,  0.6,   0.15)
    model_part.CreateNewNode(4, -0.48,  0.65,  0.0)
    model_part.CreateNewNode(5,  0.02, -0.01, -0.15)
    model_part.CreateNewNode(6, -0.03, -0.5,   0.0)
    model_part.CreateNewNode(7,  0.51,  0.02,  0.03)
    model_part.CreateNewNode(8, -0.01,  0.52, -0.05)
    model_part.CreateNewNode(9, -0.49, -0.0,   0.0)

def SetSolutionSteps(model_part):
    model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0.1)
    model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [0,7.3,4.1])
    model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [3,4,0.1])
    model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [34,2,0.1])
    model_part.GetNode(5).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [0,2,0.435])
    model_part.GetNode(6).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [34,2,0.1])
    model_part.GetNode(7).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [7,2,34.1])
    model_part.GetNode(8).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [5,4,8.1])
    model_part.GetNode(9).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [4,24,92])

def CreateElements(model_part):
    element1 = model_part.CreateNewElement("ShellKLDiscreteElement",1,[1,2,3,4,5,6,7,8,9])
    element1.SetValue(SHAPE_FUNCTION_VALUES, [0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0])

class TestIgaOutputProcess(KratosUnittest.TestCase):
    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("test_nodal_results.post.res") # usually this is deleted by the check process but not if it fails


    def test_nodal_results(self):
        test_model = KratosMultiphysics.Model()
        model_part = test_model.CreateModelPart("Structure")
        comp_model_part = model_part.CreateSubModelPart("computing_domain")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        CreateNodes(comp_model_part)

        SetSolutionSteps(comp_model_part)

        # Use the process
        # here the minimum settings are specified to test the default values!
        settings_nodal_results_process = KratosMultiphysics.Parameters("""{
              "nodal_results": [ "DISPLACEMENT" ],
              "integration_point_results": [ ],
              "output_file_name": "nodal_results.post.res",
              "model_part_name": "computing_domain",
              "file_label": "step",
              "output_control_type": "time",
              "output_frequency": 0.1
          }""")

        post_eigen_process = IgaOutputProcess(test_model, settings_nodal_results_process)

        post_eigen_process.ExecuteInitialize()
        post_eigen_process.ExecuteBeforeSolutionLoop()
        post_eigen_process.ExecuteInitializeSolutionStep()
        post_eigen_process.ExecuteFinalizeSolutionStep()
        post_eigen_process.ExecuteBeforeOutputStep()
        post_eigen_process.PrintOutput()
        post_eigen_process.ExecuteAfterOutputStep()
        post_eigen_process.ExecuteFinalize()

        # check the results
        settings_check_process = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "",
            "output_file_name"      : "",
            "comparison_type"       : "post_res_file"
        }
        """)

        settings_check_process["reference_file_name"].SetString(GetFilePath("nodal_results.ref"))
        settings_check_process["output_file_name"].SetString("nodal_results.post.res")

        check_process = CompareTwoFilesCheckProcess(settings_check_process)

        check_process.ExecuteInitialize()
        check_process.ExecuteBeforeSolutionLoop()
        check_process.ExecuteInitializeSolutionStep()
        check_process.ExecuteFinalizeSolutionStep()
        check_process.ExecuteBeforeOutputStep()
        check_process.ExecuteAfterOutputStep()
        check_process.ExecuteFinalize()

    def test_elemental_results(self):
        print("test2")
        print(KratosMultiphysics.IgaApplication.COORDINATES)
        test_model = KratosMultiphysics.Model()
        model_part = test_model.CreateModelPart("Structure")
        comp_model_part = model_part.CreateSubModelPart("computing_domain")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        CreateNodes(comp_model_part)

        SetSolutionSteps(comp_model_part)

        #CreateElements(comp_model_part)

        # Use the process
        # here the minimum settings are specified to test the default values!
        settings_nodal_results_process = KratosMultiphysics.Parameters("""{
              "nodal_results": [ ],
              "integration_point_results": [ "KratosMultiphysics.IgaApplication.COORDINATES" ],
              "output_file_name": "elemental_results.post.res",
              "model_part_name": "computing_domain",
              "file_label": "step",
              "output_control_type": "time",
              "output_frequency": 0.1
          }""")

        post_eigen_process = IgaOutputProcess(test_model, settings_nodal_results_process)

    #    post_eigen_process.ExecuteInitialize()
    #    post_eigen_process.ExecuteBeforeSolutionLoop()
    #    post_eigen_process.ExecuteInitializeSolutionStep()
    #    post_eigen_process.ExecuteFinalizeSolutionStep()
    #    post_eigen_process.ExecuteBeforeOutputStep()
    #    post_eigen_process.PrintOutput()
    #    post_eigen_process.ExecuteAfterOutputStep()
    #    post_eigen_process.ExecuteFinalize()

        ## check the results
        #settings_check_process = KratosMultiphysics.Parameters("""
        #{
        #    "reference_file_name"   : "",
        #    "output_file_name"      : "",
        #    "comparison_type"       : "post_res_file"
        #}
        #""")

        #settings_check_process["reference_file_name"].SetString(GetFilePath("test_postprocess_eigenvalues_process.ref"))
        #settings_check_process["output_file_name"].SetString("Structure_EigenResults_0.post.res")

        #check_process = CompareTwoFilesCheckProcess(settings_check_process)

        #check_process.ExecuteInitialize()
        #check_process.ExecuteBeforeSolutionLoop()
        #check_process.ExecuteInitializeSolutionStep()
        #check_process.ExecuteFinalizeSolutionStep()
        #check_process.ExecuteBeforeOutputStep()
        #check_process.ExecuteAfterOutputStep()
        #check_process.ExecuteFinalize()


if __name__ == '__main__':
    KratosUnittest.main()
