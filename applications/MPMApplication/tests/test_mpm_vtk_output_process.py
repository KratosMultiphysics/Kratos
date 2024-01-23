import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.MPMApplication.mpm_vtk_output_process as mpm_vtk_output_process
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import os
import pathlib

class TestMPMVtkOutputProcess(KratosUnittest.TestCase):
    def test_ascii_mpm_conditions_vtk_output_2D(self):
        ExecuteBasicMPMVTKOutputProcessCheck("ascii", "condition")

    def test_ascii_mpm_elements_vtk_output_2D(self):
        ExecuteBasicMPMVTKOutputProcessCheck("ascii", "element")

    def test_binary_mpm_conditions_vtk_output_2D(self):
        if os.name == "nt":
            self.skipTest("Binary output currently not working on Windows")
        ExecuteBasicMPMVTKOutputProcessCheck("binary", "condition")

    def test_binary_mpm_elements_vtk_output_2D(self):
        if os.name == "nt":
            self.skipTest("Binary output currently not working on Windows")
        ExecuteBasicMPMVTKOutputProcessCheck("binary", "element")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_mpm_vtk_output")

def GetFilePath(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName

def SetupModel2D(grid_model_part, initial_mesh_model_part, mpm_model_part):
    # Define Initial Mesh Model Part (used for defining material point elements)
    initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
    initial_mesh_sub_model_part = initial_mesh_model_part.CreateSubModelPart("SubInitialMesh")
    initial_mesh_sub_model_part.GetProperties()[1].SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)

    initial_mesh_sub_model_part.CreateNewNode(1, -0.25, -0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(2, -0.25,  0.00, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(3, -0.25, +0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(4,  0.00, -0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(5,  0.00,  0.00, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(6,  0.00, +0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(7, +0.25, -0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(8, +0.25,  0.00, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(9, +0.25, +0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 1, [1, 2, 5, 4], initial_mesh_sub_model_part.GetProperties()[1])
    initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 2, [2, 3, 6, 5], initial_mesh_sub_model_part.GetProperties()[1])
    initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 3, [4, 5, 8, 7], initial_mesh_sub_model_part.GetProperties()[1])
    initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 4, [5, 6, 9, 8], initial_mesh_sub_model_part.GetProperties()[1])

    # Define Background Grid model part
    grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
    grid_sub_model_part = grid_model_part.CreateSubModelPart("SubBackgroundGrid")
    # Nodes for background grid
    grid_sub_model_part.CreateNewNode(1,  -1.0, -0.5, 0.0)
    grid_sub_model_part.CreateNewNode(2,  -1.0,  0.0, 0.0)
    grid_sub_model_part.CreateNewNode(3,  -1.0, +0.5, 0.0)
    grid_sub_model_part.CreateNewNode(4,  -0.5, -0.5, 0.0)
    grid_sub_model_part.CreateNewNode(5,  -0.5,  0.0, 0.0)
    grid_sub_model_part.CreateNewNode(6,  -0.5, +0.5, 0.0)
    grid_sub_model_part.CreateNewNode(7,   0.0, -0.5, 0.0)
    grid_sub_model_part.CreateNewNode(8,   0.0,  0.0, 0.0)
    grid_sub_model_part.CreateNewNode(9,   0.0, +0.5, 0.0)
    grid_sub_model_part.CreateNewNode(10, +0.5, -0.5, 0.0)
    grid_sub_model_part.CreateNewNode(11, +0.5,  0.0, 0.0)
    grid_sub_model_part.CreateNewNode(12, +0.5, +0.5, 0.0)
    grid_sub_model_part.CreateNewNode(13, +1.0, -0.5, 0.0)
    grid_sub_model_part.CreateNewNode(14, +1.0,  0.0, 0.0)
    grid_sub_model_part.CreateNewNode(15, +1.0, +0.5, 0.0)
    # Elements for background grid
    grid_sub_model_part.CreateNewElement("Element2D4N", 1, [ 1,  2,  5,  4], grid_sub_model_part.GetProperties()[1])
    grid_sub_model_part.CreateNewElement("Element2D4N", 2, [ 2,  3,  6,  5], grid_sub_model_part.GetProperties()[1])
    grid_sub_model_part.CreateNewElement("Element2D4N", 3, [ 4,  5,  8,  7], grid_sub_model_part.GetProperties()[1])
    grid_sub_model_part.CreateNewElement("Element2D4N", 4, [ 5,  6,  9,  8], grid_sub_model_part.GetProperties()[1])
    grid_sub_model_part.CreateNewElement("Element2D4N", 5, [ 7,  8, 11, 10], grid_sub_model_part.GetProperties()[1])
    grid_sub_model_part.CreateNewElement("Element2D4N", 6, [ 8,  9, 12, 11], grid_sub_model_part.GetProperties()[1])
    grid_sub_model_part.CreateNewElement("Element2D4N", 7, [10, 11, 14, 13], grid_sub_model_part.GetProperties()[1])
    grid_sub_model_part.CreateNewElement("Element2D4N", 8, [11, 12, 15, 14], grid_sub_model_part.GetProperties()[1])

    # Interface
    grid_interface = grid_model_part.CreateSubModelPart("InterfaceConditions")
    # Nodes for interface condition
    grid_interface.CreateNewNode(16, -0.2, 0.0, 0.0)
    grid_interface.CreateNewNode(17, -0.1, 0.0, 0.0)
    grid_interface.CreateNewNode(18,  0.0, 0.0, 0.0)
    grid_interface.CreateNewNode(19, +0.1, 0.0, 0.0)
    grid_interface.CreateNewNode(20, +0.2, 0.0, 0.0)
    # Conditions
    grid_interface.CreateNewCondition("LineCondition2D2N", 1, [16, 17], grid_interface.GetProperties()[1])
    grid_interface.CreateNewCondition("LineCondition2D2N", 2, [17, 18], grid_interface.GetProperties()[1])
    grid_interface.CreateNewCondition("LineCondition2D2N", 3, [18, 19], grid_interface.GetProperties()[1])
    grid_interface.CreateNewCondition("LineCondition2D2N", 4, [19, 20], grid_interface.GetProperties()[1])
    KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, grid_interface.Conditions)
    for condition in grid_interface.Conditions:
        condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, 1)
        condition.SetValue(KratosMPM.MPC_BOUNDARY_CONDITION_TYPE, 1)

    # Define Material Point Model Part (this will contain mp elements)
    mpm_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

    # Activate elements of initial mesh model part
    mpm_model_part.SetNodes(grid_model_part.GetNodes())
    KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_model_part.Elements)
    # Generate Material Point Elements
    KratosMPM.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, mpm_model_part, False)
    KratosMPM.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, mpm_model_part)

def SetSolution(model_part):
    time = model_part.ProcessInfo[KratosMultiphysics.TIME] + 0.150
    step = model_part.ProcessInfo[KratosMultiphysics.STEP]

    for elem in model_part.Elements:
        coord = elem.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, model_part.ProcessInfo)[0]
        elem.SetValuesOnIntegrationPoints(KratosMPM.MP_COORD, [[coord[0]*time, coord[1]+step, coord[2]]], model_part.ProcessInfo)
        density = elem.CalculateOnIntegrationPoints(KratosMPM.MP_DENSITY, model_part.ProcessInfo)[0]
        elem.SetValuesOnIntegrationPoints(KratosMPM.MP_DENSITY, [density+0.2], model_part.ProcessInfo)
        displacement = elem.CalculateOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, model_part.ProcessInfo)[0]
        elem.SetValuesOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, [[displacement[0]+0.1*time, displacement[1]+1, displacement[2]+step/10]], model_part.ProcessInfo)

    for condition in model_part.Conditions:
        coord = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD, model_part.ProcessInfo)[0]
        condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_COORD, [[coord[0]*time, coord[1], coord[2]+step]], model_part.ProcessInfo)
        displacement = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, model_part.ProcessInfo)[0]
        condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, [[displacement[0]*time, displacement[1]+1, displacement[2]+step/10]], model_part.ProcessInfo)
        area = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_AREA, model_part.ProcessInfo)[0]
        condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_AREA, [area*2], model_part.ProcessInfo)

def Check(output_path, reference_files_path, file_name, file_format, entity_type):
    output_file = output_path/file_name
    reference_file = reference_files_path/f"{file_format}_{entity_type}"/file_name
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : "",
        "output_file_name"    : ""
    }""")
    params["reference_file_name"].SetString(str(GetFilePath(reference_file)))
    params["output_file_name"].SetString(str(output_file))
    if file_format == "ascii":
        params.AddEmptyValue("comparison_type").SetString("vtk")
    CompareTwoFilesCheckProcess(params).Execute()

def SetupMPMVtkOutputProcess(parameters, model):
    return mpm_vtk_output_process.Factory(parameters, model)

def ExecuteBasicMPMVTKOutputProcessCheck(file_format, entity_type):
    model = KratosMultiphysics.Model()
    initial_mesh = model.CreateModelPart("InitialMesh")
    background_grid = model.CreateModelPart("Background_Grid")
    mpm_model_part = model.CreateModelPart("MPMModelPart")
    SetupModel2D(background_grid, initial_mesh, mpm_model_part)

    mpm_vtk_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : "MPMModelPart",
            "file_format"                        : "ascii",
            "entity_type"                        : "",
            "output_precision"                   : 7,
            "output_interval"                    : 2,
            "output_control_type"                : "step",
            "output_sub_model_parts"             : true,
            "output_path"                        : "test_mpm_vtk_output",
            "save_output_files_in_folder"        : true,
            "gauss_point_variables_in_elements"  : [],
            "element_flags"                      : ["BOUNDARY"],
            "condition_flags"                    : ["BOUNDARY"]
        }
    }""")

    mpm_vtk_output_parameters["Parameters"]["file_format"].SetString(file_format)
    mpm_vtk_output_parameters["Parameters"]["entity_type"].SetString(entity_type)
    if entity_type == "element":
        mpm_vtk_output_parameters["Parameters"]["gauss_point_variables_in_elements"].SetStringArray(["MP_DISPLACEMENT","MP_DENSITY"])
    elif entity_type == "condition":
        mpm_vtk_output_parameters["Parameters"]["gauss_point_variables_in_elements"].SetStringArray(["MPC_DISPLACEMENT","MPC_AREA"])
    mpm_vtk_output_process = SetupMPMVtkOutputProcess(mpm_vtk_output_parameters, model)

    output_path = pathlib.Path(mpm_vtk_output_parameters["Parameters"]["output_path"].GetString())
    reference_files_path = pathlib.Path("mpm_vtk_output_process_files")

    time = 0.0
    dt = 0.2
    step = 0
    end_time = 1.0
    mpm_vtk_output_process.ExecuteInitialize()
    mpm_vtk_output_process.ExecuteBeforeSolutionLoop()
    Check(output_path, reference_files_path, f"MPMModelPart_0_{step}.vtk", file_format, entity_type)
    Check(output_path, reference_files_path, f"MPMModelPart_InterfaceConditions_0_{step}.vtk", file_format, entity_type)
    Check(output_path, reference_files_path, f"MPMModelPart_SubInitialMesh_0_{step}.vtk", file_format, entity_type)
    Check(output_path, reference_files_path, f"Background_Grid_0_{step}.vtk", file_format, entity_type)
    Check(output_path, reference_files_path, f"Background_Grid_InterfaceConditions_0_{step}.vtk", file_format, entity_type)
    Check(output_path, reference_files_path, f"Background_Grid_SubBackgroundGrid_0_{step}.vtk", file_format, entity_type)

    while (time < end_time):
        time += dt
        step += 1
        mpm_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        SetSolution(mpm_model_part)
        mpm_model_part.CloneTimeStep(time)
        mpm_vtk_output_process.ExecuteInitializeSolutionStep()
        mpm_vtk_output_process.ExecuteFinalizeSolutionStep()
        if mpm_vtk_output_process.IsOutputStep():
            mpm_vtk_output_process.ExecuteBeforeOutputStep()
            mpm_vtk_output_process.PrintOutput()
            mpm_vtk_output_process.ExecuteAfterOutputStep()
            # Compare output file with reference file
            Check(output_path, reference_files_path, f"MPMModelPart_0_{step}.vtk", file_format, entity_type)
            Check(output_path, reference_files_path, f"MPMModelPart_InterfaceConditions_0_{step}.vtk", file_format, entity_type)
            Check(output_path, reference_files_path, f"MPMModelPart_SubInitialMesh_0_{step}.vtk", file_format, entity_type)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
