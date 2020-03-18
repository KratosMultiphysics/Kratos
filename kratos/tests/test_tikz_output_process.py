from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.tikz_output_process as tikz_output_process

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestTikZOutputProcess(KratosUnittest.TestCase):
    @classmethod
    def test_ascii_tikz_output_2D(cls):
        ExecuteBasicTikZoutputProcessCheck()

    @classmethod
    def tearDown(cls):
        kratos_utils.DeleteDirectoryIfExisting("test_tikz_output")

def SetupModelPart2D(model_part):
    model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

    # Create nodes
    model_part.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
    model_part.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
    model_part.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
    model_part.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
    model_part.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
    model_part.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
    model_part.CreateNewNode(11, 1.50000, 0.50000, 0.00000)
    model_part.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
    model_part.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(14, 2.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
    model_part.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(17, 1.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(18, 1.00000, 0.00000, 0.00000)

    # Create elements
    model_part.CreateNewElement("Element2D4N", 1, [14, 11, 12, 15], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 2, [13, 10, 11, 14], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 3, [11, 17, 18, 12], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 4, [10, 16, 17, 11], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 5, [2, 4, 3, 1], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 6, [5, 8, 4, 2], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 7, [4, 7, 6, 3], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 8, [8, 9, 7, 4], model_part.GetProperties()[1])

def SetSolution(model_part):
    time = model_part.ProcessInfo[KratosMultiphysics.TIME] + 0.158
    step = model_part.ProcessInfo[KratosMultiphysics.STEP]

    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,node.X*time*step)

def SetupTikZOutputProcess(current_model, parameters):
    return tikz_output_process.Factory(parameters, current_model)

def Check(output_file,reference_file):
    import KratosMultiphysics.compare_two_files_check_process as compare_process

    ## Settings string in json format
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : "",
        "output_file_name"    : ""
    }""")
    params["reference_file_name"].SetString(GetFilePath(reference_file))
    params["output_file_name"].SetString(output_file)

    cmp_process = compare_process.CompareTwoFilesCheckProcess(params)

    cmp_process.ExecuteInitialize()
    cmp_process.ExecuteBeforeSolutionLoop()
    cmp_process.ExecuteInitializeSolutionStep()
    cmp_process.ExecuteFinalizeSolutionStep()
    cmp_process.ExecuteBeforeOutputStep()
    cmp_process.ExecuteAfterOutputStep()
    cmp_process.ExecuteFinalize()

def ExecuteBasicTikZoutputProcessCheck():
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    current_model = KratosMultiphysics.Model()
    model_part_name = "Main"
    model_part = current_model.CreateModelPart(model_part_name)
    SetupModelPart2D(model_part)

    tikz_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name" : "Main",
            "folder_name"     : "test_tikz_output"
        }
    }""")

    tikz_output_parameters["Parameters"]["model_part_name"].SetString(model_part_name)
    tikz_output_process = SetupTikZOutputProcess(current_model, tikz_output_parameters)

    # NOTE: Once postprocess works, activate this
    #time = 0.0
    #dt = 0.2
    step = 0
    #end_time = 0.1
    tikz_output_process.ExecuteInitialize()
    tikz_output_process.ExecuteBeforeSolutionLoop()

    Check(os.path.join("test_tikz_output","Main_STEP_" + str(step)+".tex"),\
                os.path.join("auxiliar_files_for_python_unittest", "tikz_output_process_ref_files", "Main_STEP_"+str(step)+".tex"))

    # NOTE: Once postprocess works, activate this
    #while (time <= end_time):
        #time = time + dt
        #step = step + 1
        #model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        #SetSolution(model_part)
        #tikz_output_process.ExecuteInitializeSolutionStep()
        #model_part.CloneTimeStep(time)
        #tikz_output_process.ExecuteFinalizeSolutionStep()
        #if tikz_output_process.IsOutputStep():
            #tikz_output_process.PrintOutput()

            #Check(os.path.join("test_tikz_output","Main_STEP_" + str(step)+".tex"),\
                #os.path.join("auxiliar_files_for_python_unittest", "tikz_output_process_ref_files", "Main_STEP_"+str(step)+".tex"))

if __name__ == '__main__':
    KratosUnittest.main()
