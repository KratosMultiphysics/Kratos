import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.MPMApplication.mpm_multiple_points_output_process as mpm_multiple_points_output_process
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import os
import pathlib

class TestMPMPointOutputProcess(KratosUnittest.TestCase):
    def test_mpm_output_process_condition_2D(self):
        ExecuteBasicMPMPointOutputProcess("condition")

    def test_mpm_output_process_element_2D(self):
        ExecuteBasicMPMPointOutputProcess("element")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_material_point_output")

def GetFilePath(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName

def SetupModel2D(grid_model_part, initial_mesh_model_part, mpm_model_part):
    # Define Initial Mesh Model Part (used for defining material point elements)
    initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
    initial_mesh_sub_model_part = initial_mesh_model_part.CreateSubModelPart("SubInitialMesh")
    # Number of material point elements for each element of the initial mesh
    initial_mesh_sub_model_part.GetProperties()[1].SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)
    # Nodes for the initial mesh
    initial_mesh_sub_model_part.CreateNewNode(1, -0.25, -0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(2, -0.25,  0.00, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(3, -0.25, +0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(4,  0.00, -0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(5,  0.00,  0.00, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(6,  0.00, +0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(7, +0.25, -0.25, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(8, +0.25,  0.00, 0.0)
    initial_mesh_sub_model_part.CreateNewNode(9, +0.25, +0.25, 0.0)
    # Elements for initial mesh
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
    # Generate Material Point Elements and Conditions
    KratosMPM.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, mpm_model_part, False)
    KratosMPM.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, mpm_model_part)

def SetSolution(model_part):
    time = model_part.ProcessInfo[KratosMultiphysics.TIME] + 0.150
    step = model_part.ProcessInfo[KratosMultiphysics.STEP] + 0.2

    for element in model_part.Elements:
        id_elem = element.Id

        coord = element.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, model_part.ProcessInfo)[0]
        updated_coord = [coord[0]*(time+id_elem), coord[1]+(step*id_elem), coord[2]]
        element.SetValuesOnIntegrationPoints(KratosMPM.MP_COORD, [updated_coord], model_part.ProcessInfo)

        density = element.CalculateOnIntegrationPoints(KratosMPM.MP_DENSITY, model_part.ProcessInfo)[0]
        updated_density = density + 0.2/id_elem
        element.SetValuesOnIntegrationPoints(KratosMPM.MP_DENSITY, [updated_density], model_part.ProcessInfo)

        displacement = element.CalculateOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, model_part.ProcessInfo)[0]
        updated_displacement = [displacement[0]/id_elem+0.1*time, displacement[1]+id_elem/(2+time), displacement[2]]
        element.SetValuesOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, [updated_displacement], model_part.ProcessInfo)

    for condition in model_part.Conditions:
        id_cond = condition.Id

        coord = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD, model_part.ProcessInfo)[0]
        updated_coord = [coord[0]*time, coord[1]+id_elem, coord[2]]
        condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_COORD, [updated_coord], model_part.ProcessInfo)

        displacement = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, model_part.ProcessInfo)[0]
        updated_displacement = [displacement[0]*time, displacement[1]+1, displacement[2]+step/10]
        condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, [updated_displacement], model_part.ProcessInfo)

        area = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_AREA, model_part.ProcessInfo)[0]
        updated_area = area + id_cond
        condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_AREA, [updated_area], model_part.ProcessInfo)

def Check(file_name, output_path, reference_files_path):
    output_file = pathlib.Path(output_path)/pathlib.Path(file_name)
    reference_file = pathlib.Path(reference_files_path)/pathlib.Path(file_name)
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : "",
        "output_file_name"    : "",
        "comparison_type"     : "deterministic"
    }""")
    params["reference_file_name"].SetString(str(GetFilePath(reference_file)))
    params["output_file_name"].SetString(str(output_file))
    CompareTwoFilesCheckProcess(params).Execute()

def ExecuteBasicMPMPointOutputProcess(entity_type):
    model = KratosMultiphysics.Model()
    initial_mesh = model.CreateModelPart("InitialMesh")
    background_grid = model.CreateModelPart("Background_Grid")
    mpm_model_part = model.CreateModelPart("MPMModelPart")
    SetupModel2D(background_grid, initial_mesh, mpm_model_part)

    parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"      : "MPMModelPart",
            "entity_type"          : "",
            "interval"             : [0.0, 1e30],
            "positions"            : [],
            "output_variables"     : [],
            "search_tolerance"     : 1e-3,
            "print_format"         : ".5e",
            "output_file_settings" : {}
        }
    }""")

    file_name = f"mp_process_output_{entity_type}"
    output_path = "test_material_point_output"

    parameters["Parameters"]["output_file_settings"].AddString("file_name",file_name)
    parameters["Parameters"]["output_file_settings"].AddString("output_path",output_path)
    parameters["Parameters"]["output_file_settings"].AddString("file_extension","dat")

    parameters["Parameters"]["entity_type"].SetString(entity_type)

    if entity_type == "element":
        positions = KratosMultiphysics.Matrix(2,3)
        # Coordinates first point
        positions[0,0] = -0.125
        positions[0,1] = -0.125
        positions[0,2] = +0.0
        # Coordinates second point
        positions[1,0] = -0.12499
        positions[1,1] = +0.12501
        positions[1,2] = +0.0
        parameters["Parameters"]["positions"].SetMatrix(positions)
        parameters["Parameters"]["output_variables"].SetStringArray(["MP_COORD","MP_DISPLACEMENT","MP_DENSITY"])

    elif entity_type == "condition":
        positions = KratosMultiphysics.Matrix(2,3)
        # Coordinates first point
        positions[0,0] = -0.0501
        positions[0,1] = +0.0
        positions[0,2] = +0.0
        # Coordinates second point
        positions[1,0] = +0.15
        positions[1,1] = -0.00001
        positions[1,2] = +0.0
        parameters["Parameters"]["positions"].SetMatrix(positions)
        parameters["Parameters"]["output_variables"].SetStringArray(["MPC_COORD","MPC_DISPLACEMENT","MPC_AREA"])

    test_process = mpm_multiple_points_output_process.Factory(parameters, model)

    time = 0.0
    dt = 0.2
    step = 0
    end_time = 1.0
    test_process.ExecuteInitialize()
    test_process.ExecuteBeforeSolutionLoop()

    while (time < end_time):
        time += dt
        step += 1
        mpm_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        SetSolution(mpm_model_part)
        mpm_model_part.CloneTimeStep(time)
        test_process.ExecuteInitializeSolutionStep()
        test_process.ExecuteFinalizeSolutionStep()
        if test_process.IsOutputStep():
            test_process.ExecuteBeforeOutputStep()
            test_process.PrintOutput()
            test_process.ExecuteAfterOutputStep()

    reference_files_path = "material_point_output_process_files"

    # Compare output file with reference file
    Check(file_name + "_1.dat", output_path, reference_files_path)
    Check(file_name + "_2.dat", output_path, reference_files_path)

    test_process.ExecuteFinalize()

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.DETAIL)
    KratosUnittest.main()
