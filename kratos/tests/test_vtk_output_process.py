from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.vtk_output_process as vtk_output_process

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestVtkOutputProcess(KratosUnittest.TestCase):
    def test_ascii_vtk_output(self):
        ExecuteBasicVTKoutputProcessCheck("ascii")

    def test_binary_vtk_output(self):
        ExecuteBasicVTKoutputProcessCheck("binary")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_vtk_output")

def SetupModelPart(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
    #create nodes
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

    #create Element
    model_part.CreateNewElement("Element2D4N", 1,
                        [14, 11, 12, 15], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 2,
                        [13, 10, 11, 14], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 3,
                        [11, 17, 18, 12], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 4,
                        [10, 16, 17, 11], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 5,
                        [2, 4, 3, 1], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 6, [5, 8, 4, 2],
                        model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 7, [4, 7, 6, 3],
                        model_part.GetProperties()[1])
    model_part.CreateNewElement("Element2D4N", 8, [8, 9, 7, 4],
                        model_part.GetProperties()[1])

    model_part.CreateNewCondition("LineCondition2D2N", 1, [1,2], model_part.GetProperties()[1])
    model_part.CreateNewCondition("LineCondition2D2N", 2, [2,5], model_part.GetProperties()[1])
    model_part.CreateNewCondition("LineCondition2D2N", 3, [13,14], model_part.GetProperties()[1])
    model_part.CreateNewCondition("LineCondition2D2N", 4, [14,15], model_part.GetProperties()[1])

    #create a submodelpart for boundary conditions
    bcs = model_part.CreateSubModelPart("FixedEdgeNodes")
    bcs.AddNodes([1, 2, 5])
    bcs.AddConditions([1,2])

    bcmn = model_part.CreateSubModelPart("MovingNodes")
    bcmn.AddNodes([13, 14, 15])
    bcmn.AddConditions([3,4])

def SetSolution(model_part):
    time = model_part.ProcessInfo[KratosMultiphysics.TIME] + 0.158
    step = model_part.ProcessInfo[KratosMultiphysics.STEP]

    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,[node.X*time,node.Y,node.Z*step])
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,[2*node.X,2*node.Y,2*node.Z])
        node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,node.X*time*step)

    for i_elem, elem in enumerate(model_part.Elements):
        elem.SetValue(KratosMultiphysics.DETERMINANT, [i_elem*0.189,time,time*step])

    for i_cond, cond in enumerate(model_part.Conditions):
        cond.SetValue(KratosMultiphysics.DENSITY, i_cond*step)
        cond.SetValue(KratosMultiphysics.YOUNG_MODULUS, i_cond*time)

def SetupVtkOutputProcess(current_model, parameters):
    return vtk_output_process.Factory(parameters, current_model)

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

def ExecuteBasicVTKoutputProcessCheck(file_format):
    current_model = KratosMultiphysics.Model()
    model_part_name = "Main"
    model_part = current_model.CreateModelPart(model_part_name)
    SetupModelPart(model_part)

    vtk_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : \"""" + model_part_name + """\",
            "file_format"                        : \"""" + file_format + """\",
            "output_precision"                   : 8,
            "output_frequency"                   : 2,
            "folder_name"                        : "test_vtk_output",
            "nodal_solution_step_data_variables" : ["PRESSURE","DISPLACEMENT", "VELOCITY"],
            "element_data_value_variables"       : ["DETERMINANT"],
            "condition_data_value_variables"     : ["DENSITY", "YOUNG_MODULUS"]
        }
    }""")

    vtk_output_process = SetupVtkOutputProcess(current_model, vtk_output_parameters)

    time = 0.0
    dt = 0.2
    step = 0
    end_time = 1.0
    vtk_output_process.ExecuteInitialize()

    while (time <= end_time):
        time = time + dt
        step = step + 1
        model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        #print("STEP :: ", step, ", TIME :: ", time)
        SetSolution(model_part)
        vtk_output_process.ExecuteInitializeSolutionStep()
        model_part.CloneTimeStep(time)
        vtk_output_process.ExecuteFinalizeSolutionStep()
        if vtk_output_process.IsOutputStep():
            vtk_output_process.PrintOutput()

            Check(os.path.join("test_vtk_output","Main_0_" + str(step)+".vtk"),\
                os.path.join("vtk_output_process_ref_files", file_format, "Main_0_"+str(step)+".vtk"))

            Check(os.path.join("test_vtk_output","Main_FixedEdgeNodes_0_" + str(step)+".vtk"),\
                os.path.join("vtk_output_process_ref_files", file_format, "Main_FixedEdgeNodes_0_"+str(step)+".vtk"))

            Check(os.path.join("test_vtk_output","Main_MovingNodes_0_"+str(step)+".vtk"),\
                os.path.join("vtk_output_process_ref_files", file_format, "Main_MovingNodes_0_"+str(step)+".vtk"))


if __name__ == '__main__':
    KratosUnittest.main()