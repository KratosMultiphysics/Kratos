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
    """
    Unit tests for the Ensight output process in Kratos.
    This test class contains methods to verify the correct functionality of the Ensight output process
    for various element types and dimensions in ASCII format. The tests cover 2D, 3D, quadratic 3D,
    quadratic prism 3D, and quadratic hexahedra 3D cases. Binary output tests are present but commented out,
    indicating that binary support is implemented but not yet tested.
    Each test method calls `ExecuteBasicEnsightOutputProcessCheck` with the appropriate format and geometry type.
    The `tearDown` method ensures that the output directory ("test_ensight_output") is deleted after each test run
    to maintain a clean test environment.
    """
    def test_ascii_ensight_output_2D(self):
        ExecuteBasicEnsightOutputProcessCheck("ascii", "2D")

    def test_ascii_ensight_output_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("ascii", "3D")

    def test_ascii_ensight_output_quad_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("ascii", "Quad3D")

    def test_ascii_ensight_output_quad_prism_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("ascii", "QuadraticPrism3D")

    def test_ascii_ensight_output_quad_hexahedra_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("ascii", "QuadraticHexahedra3D")

    # TODO: Binary is implemented, but requires testing
    # def test_binary_ensight_output_2D(self):
    #     ExecuteBasicEnsightOutputProcessCheck("binary", "2D")

    # def test_binary_ensight_output_3D(self):
    #     ExecuteBasicEnsightOutputProcessCheck("binary", "3D")

    # def test_binary_ensight_output_quad_3D(self):
    #     ExecuteBasicEnsightOutputProcessCheck("binary", "Quad3D")

    # def test_binary_ensight_output_quad_prism_3D(self):
    #     ExecuteBasicEnsightOutputProcessCheck("binary", "QuadraticPrism3D")

    # def test_binary_ensight_output_quad_hexahedra_3D(self):
    #     ExecuteBasicEnsightOutputProcessCheck("binary", "QuadraticHexahedra3D")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_ensight_output")

def SetupModelPart2D(model_part):
    """
    Sets up a ModelPart with a 2D mesh.

    Args:
        model_part (kratos.ModelPart): The model part to be initialized with the 2D mesh.
    """
    test_vtk_output_process.SetupModelPart2D(model_part)

def SetupModelPart3D(model_part):
    """
    Sets up a ModelPart with a 3D mesh.

    Args:
        model_part (kratos.ModelPart): The model part to be initialized with the 3D mesh.
    """
    test_vtk_output_process.SetupModelPart3D(model_part)

def SetupModelPartQuadratic3D(model_part):
    """
    Sets up a ModelPart with a 3D quadratic mesh.

    Args:
        model_part (kratos.ModelPart): The model part to be initialized with the quadratic mesh.
    """
    test_vtk_output_process.SetupModelPartQuadratic3D(model_part)

def SetupModelPartQuadraticPrism3D(model_part):
    """
    Sets up a ModelPart with a 3D quadratic prism mesh.

    Args:
        model_part (kratos.ModelPart): The model part to be initialized with the quadratic prism mesh.
    """
    test_vtk_output_process.SetupModelPartQuadraticPrism3D(model_part)

def SetupModelPartQuadraticHexahedra3D(model_part):
    """
    Sets up a ModelPart with a 3D quadratic hexahedral mesh.

    Args:
        model_part (kratos.ModelPart): The model part to be initialized with the quadratic hexahedral mesh.
    """
    test_vtk_output_process.SetupModelPartQuadraticHexahedra3D(model_part)

def SetSolution(model_part):
    """
    Sets a sample solution on the provided model part.

    Args:
        model_part (kratos.ModelPart): The model part on which to set the solution.
    """
    test_vtk_output_process.SetSolution(model_part)

def SetupEnsightOutputProcess(current_model, parameters):
    """
    Sets up the Ensight output process with the given model and parameters.

    Args:
        current_model (kratos.Model): The current Kratos model.
        parameters (kratos.Parameters): The parameters for configuring the Ensight output process.
    """
    return ensight_output_process.Factory(parameters, current_model)

def Check(output_file,reference_file, file_format, extension):
    """
    Compares an output file with a reference file using the specified file format and extension.
    Args:
        output_file (str): The path to the output file to be checked.
        reference_file (str): The path to the reference file to compare against.
        file_format (str): The format of the files (e.g., "ascii").
        extension (str): The file extension or comparison type to use for the check.
    """
    ## Settings string in json format
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : "",
        "output_file_name"    : ""
    }""")
    params["reference_file_name"].SetString(GetFilePath(reference_file))
    params["output_file_name"].SetString(output_file)
    if file_format == "ascii":
        params.AddEmptyValue("comparison_type").SetString(extension)

    CompareTwoFilesCheckProcess(params).Execute()

def ExecuteBasicEnsightOutputProcessCheck(file_format = "ascii", setup = "2D"):
    """
    Executes a basic check of the Ensight output process for various model setups and file formats.
    This function sets up a Kratos model part according to the specified geometry (2D, 3D, or various quadratic elements),
    configures the Ensight output process with a set of parameters, and runs a simple time-stepping loop.
    At each output step, it generates Ensight files and compares them against reference files to validate correctness.
    Parameters
    """
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    current_model = KratosMultiphysics.Model()
    model_part_name = "Main"
    model_part = current_model.CreateModelPart(model_part_name)
    if setup == "2D":
        SetupModelPart2D(model_part)
    elif setup == "3D":
        SetupModelPart3D(model_part)
    elif setup == "Quad3D":
        SetupModelPartQuadratic3D(model_part)
    elif setup == "QuadraticPrism3D":
        SetupModelPartQuadraticPrism3D(model_part)
    elif setup == "QuadraticHexahedra3D":
        SetupModelPartQuadraticHexahedra3D(model_part)
    else:
        raise Exception("Unknown setup: " + setup)

    ensight_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : "Main",
            "ensight_file_format"                : "6", // Options: "6", "gold"
            "file_format"                        : "ascii",
            "output_precision"                   : 6,
            "output_interval"                    : 2,
            "output_sub_model_parts"             : true,
            "evolving_geometry"                  : true,
            "output_path"                        : "test_ensight_output",
            "nodal_solution_step_data_variables" : ["PRESSURE","DISPLACEMENT", "VELOCITY"],
            "nodal_flags"                        : ["BOUNDARY"],
            "element_data_value_variables"       : ["DETERMINANT"],
            "condition_data_value_variables"     : ["DENSITY", "YOUNG_MODULUS"],
            "condition_flags"                    : ["BOUNDARY"]
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

    while time <= end_time:
        time = time + dt
        model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        SetSolution(model_part)
        ensight_output_process.ExecuteInitializeSolutionStep()
        model_part.CloneTimeStep(time)
        ensight_output_process.ExecuteFinalizeSolutionStep()
        if ensight_output_process.IsOutputStep():
            ensight_output_process.PrintOutput()

            label_step = "000" + str(step)

            Check(os.path.join("test_ensight_output","Main." + label_step + ".geo"),\
                os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main." + label_step + ".geo"), file_format, "geo")

            # Check node variables
            for variable_name in ensight_output_parameters["Parameters"]["nodal_solution_step_data_variables"].GetStringArray():
                is_scalar = KratosMultiphysics.KratosGlobals.Kernel.HasDoubleVariable(variable_name)
                name_extension = "scl" if is_scalar else "vec"
                Check(os.path.join("test_ensight_output","Main." + variable_name + "." + label_step + ".node." + name_extension),\
                    os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main." + variable_name + "." + label_step + ".node." + name_extension), file_format, "ensight_solution")

            # Check nodal flags
            for flag_name in ensight_output_parameters["Parameters"]["nodal_flags"].GetStringArray():
                Check(os.path.join("test_ensight_output","Main." + flag_name + "." + label_step + ".node.scl"),\
                    os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main." + flag_name + "." + label_step + ".node.scl"), file_format, "ensight_solution")

            # Check elemental variables
            for variable_name in ensight_output_parameters["Parameters"]["element_data_value_variables"].GetStringArray():
                is_scalar = KratosMultiphysics.KratosGlobals.Kernel.HasDoubleVariable(variable_name)
                name_extension = "scl" if is_scalar else "vec"
                Check(os.path.join("test_ensight_output","Main." + variable_name + "." + label_step + ".elem." + name_extension),\
                    os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main." + variable_name + "." + label_step + ".elem." + name_extension), file_format, "ensight_solution")

            # Check conditional variables
            for variable_name in ensight_output_parameters["Parameters"]["condition_data_value_variables"].GetStringArray():
                is_scalar = KratosMultiphysics.KratosGlobals.Kernel.HasDoubleVariable(variable_name)
                name_extension = "scl" if is_scalar else "vec"
                Check(os.path.join("test_ensight_output","Main." + variable_name + "." + label_step + ".cond." + name_extension),\
                    os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main." + variable_name + "." + label_step + ".cond." + name_extension), file_format, "ensight_solution")

            # Check conditional flags
            for flag_name in ensight_output_parameters["Parameters"]["condition_flags"].GetStringArray():
                Check(os.path.join("test_ensight_output","Main." + flag_name + "." + label_step + ".cond.scl"),\
                    os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main." + flag_name + "." + label_step + ".cond.scl"), file_format, "ensight_solution")

            step = step + 1

    Check(os.path.join("test_ensight_output","Main.case"),\
          os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main.case"), file_format, "case")

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
