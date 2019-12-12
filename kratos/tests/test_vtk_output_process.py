from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.vtk_output_process as vtk_output_process
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestVtkOutputProcess(KratosUnittest.TestCase):
    def test_ascii_vtk_output_2D(self):
        ExecuteBasicVTKoutputProcessCheck("ascii", "2D")

    def test_binary_vtk_output_2D(self):
        ExecuteBasicVTKoutputProcessCheck("binary", "2D")

    def test_ascii_vtk_output_3D(self):
        ExecuteBasicVTKoutputProcessCheck("ascii", "3D")

    def test_binary_vtk_output_3D(self):
        ExecuteBasicVTKoutputProcessCheck("binary", "3D")

    def test_ascii_vtk_output_quad_3D(self):
        ExecuteBasicVTKoutputProcessCheck("ascii", "Quad3D")

    def test_binary_vtk_output_quad_3D(self):
        ExecuteBasicVTKoutputProcessCheck("binary", "Quad3D")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_vtk_output")

def SetupModelPart2D(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
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

    # Create conditions
    model_part.CreateNewCondition("LineCondition2D2N", 1, [1,2], model_part.GetProperties()[1])
    model_part.CreateNewCondition("LineCondition2D2N", 2, [2,5], model_part.GetProperties()[1])
    model_part.CreateNewCondition("LineCondition2D2N", 3, [13,14], model_part.GetProperties()[1])
    model_part.CreateNewCondition("LineCondition2D2N", 4, [14,15], model_part.GetProperties()[1])

    # Create a submodelpart for boundary conditions
    bcs = model_part.CreateSubModelPart("FixedEdgeNodes")
    bcs.AddNodes([1, 2, 5])
    bcs.AddConditions([1,2])

    # Setting flags
    for node in bcs.Nodes:
        node.Set(KratosMultiphysics.BOUNDARY)
    for cond in bcs.Conditions:
        cond.Set(KratosMultiphysics.BOUNDARY)

    bcmn = model_part.CreateSubModelPart("MovingNodes")
    bcmn.AddNodes([13, 14, 15])
    bcmn.AddConditions([3,4])

def SetupModelPart3D(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

    # Create nodes
    model_part.CreateNewNode(1, 0.0 , 1.0 , 1.0)
    model_part.CreateNewNode(2, 0.0 , 1.0 , 0.0)
    model_part.CreateNewNode(3, 0.0 , 0.0 , 1.0)
    model_part.CreateNewNode(4, 1.0 , 1.0 , 1.0)
    model_part.CreateNewNode(5, 0.0 , 0.0 , 0.0)
    model_part.CreateNewNode(6, 1.0 , 1.0 , 0.0)
    model_part.CreateNewNode(7, 1.0 , 0.0 , 1.0)
    model_part.CreateNewNode(8, 1.0 , 0.0 , 0.0)
    model_part.CreateNewNode(9, 2.0 , 1.0 , 1.0)
    model_part.CreateNewNode(10, 2.0 , 1.0 , 0.0)
    model_part.CreateNewNode(11, 2.0 , 0.0 , 1.0)
    model_part.CreateNewNode(12, 2.0 , 0.0 , 0.0)
    model_part.CreateNewNode(13, 0.0 , 0.0 , 2.0)
    model_part.CreateNewNode(14, 1.0 , 0.0 , 2.0)
    model_part.CreateNewNode(15, 1.0 , 1.0 , 2.0)

    # Create elements
    model_part.CreateNewElement("Element3D4N", 1, [12, 10, 8, 9], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 2, [4, 6, 9, 7], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 3, [11, 7, 9, 8], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 4, [5, 3, 8, 6], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 5, [4, 6, 7, 3], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 6, [2, 3, 5, 6], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 7, [10, 9, 6, 8], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 8, [7, 8, 3, 6], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 9, [7, 8, 6, 9], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 10, [4, 1, 6, 3], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 11, [9, 12, 11, 8], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D4N", 12, [3, 2, 1, 6], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D6N", 13, [3, 7, 4, 13, 14, 15], model_part.GetProperties()[1])

    # Create conditions
    model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
    model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [3,4,5], model_part.GetProperties()[1])
    model_part.CreateNewCondition("SurfaceCondition3D3N", 3, [10,11,12], model_part.GetProperties()[1])
    model_part.CreateNewCondition("SurfaceCondition3D3N", 4, [8,9,10], model_part.GetProperties()[1])

    # Create a submodelpart for boundary conditions
    bcs = model_part.CreateSubModelPart("FixedEdgeNodes")
    bcs.AddNodes([1, 2, 5])
    bcs.AddConditions([1,2])

    # Setting flags
    for node in bcs.Nodes:
        node.Set(KratosMultiphysics.BOUNDARY)
    for cond in bcs.Conditions:
        cond.Set(KratosMultiphysics.BOUNDARY)

    bcmn = model_part.CreateSubModelPart("MovingNodes")
    bcmn.AddNodes([13, 14, 15])
    bcmn.AddConditions([3,4])

def SetupModelPartQuadratic3D(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

    # Create nodes
    model_part.CreateNewNode(1, 1.0000000000, 0.0000000000, 1.0000000000)
    model_part.CreateNewNode(2, 0.5000000000, 0.0000000000, 1.0000000000)
    model_part.CreateNewNode(3, 1.0000000000, 0.5000000000, 1.0000000000)
    model_part.CreateNewNode(4, 1.0000000000, 0.0000000000, 0.5000000000)
    model_part.CreateNewNode(5, 0.5000000000, 0.0000000000, 0.5000000000)
    model_part.CreateNewNode(6, 1.0000000000, 0.5000000000, 0.5000000000)
    model_part.CreateNewNode(7, 0.5000000000, 0.5000000000, 1.0000000000)
    model_part.CreateNewNode(8, 0.5000000000, 0.5000000000, 0.5000000000)
    model_part.CreateNewNode(9, 0.2019246055, 0.3959160307, 0.6930668948)
    model_part.CreateNewNode(10, 0.7019246055, 0.8959160307, 0.6930668948)
    model_part.CreateNewNode(11, 1.0000000000, 0.0000000000, 0.0000000000)
    model_part.CreateNewNode(12, 0.0000000000, 0.0000000000, 1.0000000000)
    model_part.CreateNewNode(13, 1.0000000000, 1.0000000000, 1.0000000000)
    model_part.CreateNewNode(14, 0.5000000000, 0.0000000000, 0.0000000000)
    model_part.CreateNewNode(15, 1.0000000000, 0.5000000000, 0.0000000000)
    model_part.CreateNewNode(16, 0.5000000000, 1.0000000000, 1.0000000000)
    model_part.CreateNewNode(17, 0.0000000000, 0.5000000000, 1.0000000000)
    model_part.CreateNewNode(18, 0.0000000000, 0.0000000000, 0.5000000000)
    model_part.CreateNewNode(19, 1.0000000000, 1.0000000000, 0.5000000000)
    model_part.CreateNewNode(20, 0.4038492111, 0.7918320615, 0.3861337896)
    model_part.CreateNewNode(21, 0.2019246055, 0.3959160307, 0.1930668948)
    model_part.CreateNewNode(22, 0.5000000000, 0.5000000000, 0.0000000000)
    model_part.CreateNewNode(23, 0.5000000000, 1.0000000000, 0.5000000000)
    model_part.CreateNewNode(24, 0.0000000000, 0.5000000000, 0.5000000000)
    model_part.CreateNewNode(25, 0.2019246055, 0.8959160307, 0.6930668948)
    model_part.CreateNewNode(26, 0.7019246055, 0.8959160307, 0.1930668948)
    model_part.CreateNewNode(27, 0.0000000000, 0.0000000000, 0.0000000000)
    model_part.CreateNewNode(28, 1.0000000000, 1.0000000000, 0.0000000000)
    model_part.CreateNewNode(29, 0.0000000000, 1.0000000000, 1.0000000000)
    model_part.CreateNewNode(30, 0.2019246055, 0.8959160307, 0.1930668948)
    model_part.CreateNewNode(31, 0.5000000000, 1.0000000000, 0.0000000000)
    model_part.CreateNewNode(32, 0.0000000000, 0.5000000000, 0.0000000000)
    model_part.CreateNewNode(33, 0.0000000000, 1.0000000000, 0.5000000000)
    model_part.CreateNewNode(34, 0.0000000000, 1.0000000000, 0.0000000000)

    # Create elements
    model_part.CreateNewElement("Element3D10N", 1,  [27, 11, 28, 13, 14, 15, 22,  8,  6, 19], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 2,  [34, 12, 27, 20, 24, 18, 32, 30,  9, 21], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 3,  [27, 34, 20, 28, 32, 30, 21, 22, 31, 26], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 4,  [34, 20, 28, 29, 30, 26, 31, 33, 25, 23], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 5,  [34, 20, 29, 12, 30, 25, 33, 24,  9, 17], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 6,  [20, 27, 28, 13, 21, 22, 26, 10,  8, 19], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 7,  [28, 20, 13, 29, 26, 10, 19, 23, 25, 16], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 8,  [20, 13, 29, 12, 10, 16, 25,  9,  7, 17], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 9,  [27,  1, 11, 13,  5,  4, 14,  8,  3,  6], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 10, [12, 13,  1, 27,  7,  3,  2, 18,  8,  5], model_part.GetProperties()[1])
    model_part.CreateNewElement("Element3D10N", 11, [12, 27, 20, 13, 18, 21,  9,  7,  8, 10], model_part.GetProperties()[1])

    # Create conditions
    model_part.CreateNewCondition("SurfaceCondition3D6N", 1, [12, 1, 13, 2, 3, 7], model_part.GetProperties()[1])
    model_part.CreateNewCondition("SurfaceCondition3D6N", 2, [13, 29, 12, 16, 17, 7], model_part.GetProperties()[1])
    model_part.CreateNewCondition("SurfaceCondition3D6N", 3, [27, 11, 28, 14, 15, 21], model_part.GetProperties()[1])
    model_part.CreateNewCondition("SurfaceCondition3D6N", 4, [28, 34, 27, 31, 32, 21], model_part.GetProperties()[1])

    # Create a submodelpart for boundary conditions
    bcs = model_part.CreateSubModelPart("FixedEdgeNodes")
    bcs.AddNodes([1,2,3,7,12,13,16,17,29])
    bcs.AddConditions([1,2])

    bcmn = model_part.CreateSubModelPart("MovingNodes")
    bcmn.AddNodes([11,14,15,21,27,28,31,32,34])
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

def Check(output_file,reference_file, file_format):
    ## Settings string in json format
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : "",
        "output_file_name"    : ""
    }""")
    params["reference_file_name"].SetString(GetFilePath(reference_file))
    params["output_file_name"].SetString(output_file)
    if file_format == "ascii":
        params.AddEmptyValue("comparison_type").SetString("vtk")

    CompareTwoFilesCheckProcess(params).Execute()

def ExecuteBasicVTKoutputProcessCheck(file_format = "ascii", setup = "2D"):
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    current_model = KratosMultiphysics.Model()
    model_part_name = "Main"
    model_part = current_model.CreateModelPart(model_part_name)
    if setup == "2D":
        SetupModelPart2D(model_part)
    elif setup == "3D":
        SetupModelPart3D(model_part)
    else:
        SetupModelPartQuadratic3D(model_part)

    vtk_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : "Main",
            "file_format"                        : "ascii",
            "output_precision"                   : 8,
            "output_frequency"                   : 2,
            "output_sub_model_parts"             : true,
            "folder_name"                        : "test_vtk_output",
            "nodal_solution_step_data_variables" : ["PRESSURE","DISPLACEMENT", "VELOCITY"],
            "nodal_flags"                        : ["BOUNDARY"],
            "element_data_value_variables"       : ["DETERMINANT"],
            "condition_data_value_variables"     : ["DENSITY", "YOUNG_MODULUS"],
            "condition_flags"                    : ["BOUNDARY"]
        }
    }""")

    vtk_output_parameters["Parameters"]["model_part_name"].SetString(model_part_name)
    vtk_output_parameters["Parameters"]["file_format"].SetString(file_format)
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
        SetSolution(model_part)
        vtk_output_process.ExecuteInitializeSolutionStep()
        model_part.CloneTimeStep(time)
        vtk_output_process.ExecuteFinalizeSolutionStep()
        if vtk_output_process.IsOutputStep():
            vtk_output_process.PrintOutput()

            Check(os.path.join("test_vtk_output","Main_0_" + str(step)+".vtk"),\
                os.path.join("auxiliar_files_for_python_unittest", "vtk_output_process_ref_files", file_format + setup, "Main_0_"+str(step)+".vtk"), file_format)

            Check(os.path.join("test_vtk_output","Main_FixedEdgeNodes_0_" + str(step)+".vtk"),\
                os.path.join("auxiliar_files_for_python_unittest", "vtk_output_process_ref_files", file_format + setup, "Main_FixedEdgeNodes_0_"+str(step)+".vtk"), file_format)

            Check(os.path.join("test_vtk_output","Main_MovingNodes_0_"+str(step)+".vtk"),\
                os.path.join("auxiliar_files_for_python_unittest", "vtk_output_process_ref_files", file_format + setup, "Main_MovingNodes_0_"+str(step)+".vtk"), file_format)


if __name__ == '__main__':
    KratosUnittest.main()
