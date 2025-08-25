import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.ensight_output_process as ensight_output_process
import test_vtk_output_process
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestEnsightOutputProcess(KratosUnittest.TestCase):
    def test_ascii_ensight_output_2D(self):
        ExecuteBasicVTKoutputProcessCheck("ascii", "2D")

    # def test_binary_ensight_output_2D(self):
    #     ExecuteBasicVTKoutputProcessCheck("binary", "2D")

    # def test_ascii_ensight_output_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("ascii", "3D")

    # def test_binary_ensight_output_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("binary", "3D")

    # def test_ascii_ensight_output_quad_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("ascii", "Quad3D")

    # def test_binary_ensight_output_quad_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("binary", "Quad3D")

    # def test_ascii_ensight_output_quad_prism_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("ascii", "QuadraticPrism3D")

    # def test_binary_ensight_output_quad_prism_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("binary", "QuadraticPrism3D")

    # def test_ascii_ensight_output_quad_hexahedra_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("ascii", "QuadraticHexahedra3D")

    # def test_binary_ensight_output_quad_hexahedra_3D(self):
    #     ExecuteBasicVTKoutputProcessCheck("binary", "QuadraticHexahedra3D")

    # def test_ascii_ensight_output_quad_hexahedra_3D_27N(self):
    #     ExecuteBasicVTKoutputProcessCheck("ascii", "QuadraticHexahedra3D27N")

    # def test_binary_ensight_output_quad_hexahedra_3D_27N(self):
    #     ExecuteBasicVTKoutputProcessCheck("binary", "QuadraticHexahedra3D27N")

    def tearDown(self):
        # kratos_utils.DeleteDirectoryIfExisting("test_ensight_output")
        pass

def SetupModelPart2D(model_part):
    test_vtk_output_process.SetupModelPart2D(model_part)

def SetupModelPart3D(model_part):
    test_vtk_output_process.SetupModelPart3D(model_part)

def SetupModelPartQuadratic3D(model_part):
    test_vtk_output_process.SetupModelPartQuadratic3D(model_part)

def SetupModelPartQuadraticPrism3D(model_part):
    test_vtk_output_process.SetupModelPartQuadraticPrism3D(model_part)

def SetupModelPartQuadraticHexahedra3D(model_part):
    test_vtk_output_process.SetupModelPartQuadraticHexahedra3D(model_part)

def SetupModelPartHexahedra3D27N(model_part):
    test_vtk_output_process.SetupModelPartHexahedra3D27N(model_part)

def SetSolution(model_part):
    test_vtk_output_process.SetSolution(model_part)

def SetupEnsightOutputProcess(current_model, parameters):
    return ensight_output_process.Factory(parameters, current_model)

def Check(output_file,reference_file, file_format):
    ## Settings string in json format
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : "",
        "output_file_name"    : ""
    }""")
    params["reference_file_name"].SetString(GetFilePath(reference_file))
    params["output_file_name"].SetString(output_file)
    if file_format == "ascii":
        params.AddEmptyValue("comparison_type").SetString("ensight")

    CompareTwoFilesCheckProcess(params).Execute()

def ExecuteBasicVTKoutputProcessCheck(file_format = "ascii", setup = "2D"):
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    current_model = KratosMultiphysics.Model()
    model_part_name = "Main"
    model_part = current_model.CreateModelPart(model_part_name)
    bc_defined = True
    if setup == "2D":
        SetupModelPart2D(model_part)
    elif setup == "3D":
        SetupModelPart3D(model_part)
    elif setup == "Quad3D":
        SetupModelPartQuadratic3D(model_part)
    elif setup == "QuadraticPrism3D":
        SetupModelPartQuadraticPrism3D(model_part)
        bc_defined = False
    elif setup == "QuadraticHexahedra3D":
        SetupModelPartQuadraticHexahedra3D(model_part)
        bc_defined = False
    elif setup == "QuadraticHexahedra3D27N":
        SetupModelPartHexahedra3D27N(model_part)
        bc_defined = False
    else:
        raise Exception("Unknown setup: " + setup)

    ensight_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : "Main",
            "ensight_file_format"                : "6", // Options: "6", "gold"
            "file_format"                        : "ascii",
            "output_precision"                   : 8,
            "output_interval"                    : 2,
            "output_sub_model_parts"             : true,
            "evolving_geometry"                  : true,
            "output_path"                        : "test_ensight_output",
            "nodal_solution_step_data_variables" : ["PRESSURE","DISPLACEMENT", "VELOCITY"],
            "nodal_flags"                        : ["BOUNDARY"]//, // NOTE: Implemented but giving issues in Paraview
            //"element_data_value_variables"       : ["DETERMINANT"],
            //"condition_data_value_variables"     : ["DENSITY", "YOUNG_MODULUS"],
            //"condition_flags"                    : ["BOUNDARY"]
        }
    }""")

    ensight_output_parameters["Parameters"]["model_part_name"].SetString(model_part_name)
    ensight_output_parameters["Parameters"]["file_format"].SetString(file_format)
    ensight_output_process = SetupEnsightOutputProcess(current_model, ensight_output_parameters)

    time = 0.0
    dt = 0.2
    step = 0
    end_time = 1.0
    ensight_output_process.ExecuteInitialize()

    while (time <= end_time):
        time = time + dt
        step = step + 1
        model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        SetSolution(model_part)
        ensight_output_process.ExecuteInitializeSolutionStep()
        model_part.CloneTimeStep(time)
        ensight_output_process.ExecuteFinalizeSolutionStep()
        if ensight_output_process.IsOutputStep():
            ensight_output_process.PrintOutput()

            # Check(os.path.join("test_ensight_output","Main_0_" + str(step)+".case"),\
            #     os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main_0_"+str(step)+".case"), file_format)

            # if bc_defined:
            #     Check(os.path.join("test_ensight_output","Main_FixedEdgeNodes_0_" + str(step)+".case"),\
            #         os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main_FixedEdgeNodes_0_"+str(step)+".case"), file_format)

            #     Check(os.path.join("test_ensight_output","Main_MovingNodes_0_"+str(step)+".case"),\
            #         os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main_MovingNodes_0_"+str(step)+".case"), file_format)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
