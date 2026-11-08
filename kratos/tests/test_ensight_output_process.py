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
    for various element types and dimensions in both ASCII and binary formats. The tests cover 2D, 3D,
    quadratic 3D, quadratic prism 3D, and quadratic hexahedra 3D cases.
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

    def test_binary_ensight_output_2D(self):
        ExecuteBasicEnsightOutputProcessCheck("binary", "2D")

    def test_binary_ensight_output_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("binary", "3D")

    def test_binary_ensight_output_quad_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("binary", "Quad3D")

    def test_binary_ensight_output_quad_prism_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("binary", "QuadraticPrism3D")

    def test_binary_ensight_output_quad_hexahedra_3D(self):
        ExecuteBasicEnsightOutputProcessCheck("binary", "QuadraticHexahedra3D")

    # --- Symmetric tensor (CAUCHY_STRESS_TENSOR, 3x3 symmetric Matrix) ---

    def test_ascii_ensight6_symmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "6", "2D", "symmetric")

    def test_binary_ensight6_symmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "6", "2D", "symmetric")

    def test_ascii_ensightgold_symmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "gold", "2D", "symmetric")

    def test_binary_ensightgold_symmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "gold", "2D", "symmetric")

    def test_ascii_ensight6_symmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "6", "3D", "symmetric")

    def test_binary_ensight6_symmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "6", "3D", "symmetric")

    def test_ascii_ensightgold_symmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "gold", "3D", "symmetric")

    def test_binary_ensightgold_symmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "gold", "3D", "symmetric")

    # --- Asymmetric tensor (DEFORMATION_GRADIENT, 3x3 asymmetric Matrix) ---

    def test_ascii_ensight6_asymmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "6", "2D", "asymmetric")

    def test_binary_ensight6_asymmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "6", "2D", "asymmetric")

    def test_ascii_ensightgold_asymmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "gold", "2D", "asymmetric")

    def test_binary_ensightgold_asymmetric_tensor_2D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "gold", "2D", "asymmetric")

    def test_ascii_ensight6_asymmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "6", "3D", "asymmetric")

    def test_binary_ensight6_asymmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "6", "3D", "asymmetric")

    def test_ascii_ensightgold_asymmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("ascii", "gold", "3D", "asymmetric")

    def test_binary_ensightgold_asymmetric_tensor_3D(self):
        ExecuteEnsightTensorOutputProcessCheck("binary", "gold", "3D", "asymmetric")

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
        file_format (str): The format of the files (e.g., "ascii" or "binary").
        extension (str): The file extension or comparison type to use for the check.
    """
    ## Settings string in json format
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : "",
        "output_file_name"    : ""
    }""")
    params["reference_file_name"].SetString(GetFilePath(reference_file))
    params["output_file_name"].SetString(output_file)
    # The .case file is always ASCII text regardless of file_format.
    # For ascii format, use the named comparison type for all files.
    # For binary format, binary data files use exact byte comparison (default).
    if file_format == "ascii" or extension == "case":
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
            "ensight_file_format"                : "6", // Options: "5", "6", "gold"
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
          os.path.join("auxiliar_files_for_python_unittest", "ensight_output_process_ref_files", file_format + setup, "Main.case"), "ascii", "case")

def SetSolutionWithMatrixVariables(model_part, tensor_type):
    """
    Sets matrix-valued variables on elements and conditions for tensor output testing.

    For 'symmetric' tensor_type, sets CAUCHY_STRESS_TENSOR (symmetric 3x3 Matrix) on elements.
    For 'asymmetric' tensor_type, sets DEFORMATION_GRADIENT (asymmetric 3x3 Matrix) on elements.

    Args:
        model_part (kratos.ModelPart): The model part on which to set the solution.
        tensor_type (str): Either "symmetric" or "asymmetric".
    """
    time = model_part.ProcessInfo[KratosMultiphysics.TIME] + 0.158
    step = model_part.ProcessInfo[KratosMultiphysics.STEP]

    # Also set standard solution step variables so all nodes have values
    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, [node.X * time, node.Y, node.Z * step])
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, [2 * node.X, 2 * node.Y, 2 * node.Z])
        node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, node.X * time * step)

    for i_elem, elem in enumerate(model_part.Elements):
        m = KratosMultiphysics.Matrix(3, 3)
        cx = float(i_elem + 1) * 0.1
        if tensor_type == "symmetric":
            # CAUCHY_STRESS_TENSOR: symmetric 3x3
            m[0, 0] = cx;       m[0, 1] = cx * 0.5; m[0, 2] = cx * 0.2
            m[1, 0] = cx * 0.5; m[1, 1] = cx * 2.0; m[1, 2] = cx * 0.3
            m[2, 0] = cx * 0.2; m[2, 1] = cx * 0.3; m[2, 2] = cx * 3.0
            elem.SetValue(KratosMultiphysics.CAUCHY_STRESS_TENSOR, m)
        else:
            # DEFORMATION_GRADIENT: asymmetric 3x3
            m[0, 0] = 1.0 + cx; m[0, 1] = cx * 0.1; m[0, 2] = cx * 0.2
            m[1, 0] = cx * 0.3; m[1, 1] = 1.0 + cx; m[1, 2] = cx * 0.4
            m[2, 0] = cx * 0.5; m[2, 1] = cx * 0.6; m[2, 2] = 1.0 + cx
            elem.SetValue(KratosMultiphysics.DEFORMATION_GRADIENT, m)

    for i_cond, cond in enumerate(model_part.Conditions):
        cond.SetValue(KratosMultiphysics.DENSITY, i_cond * step)
        cond.SetValue(KratosMultiphysics.YOUNG_MODULUS, i_cond * time)

def CheckFileNonEmpty(filepath):
    """Asserts that the file at filepath exists and has non-zero size."""
    if not os.path.isfile(filepath):
        raise AssertionError(f"Expected output file not found: {filepath}")
    if os.path.getsize(filepath) == 0:
        raise AssertionError(f"Output file is empty: {filepath}")

def ExecuteEnsightTensorOutputProcessCheck(file_format="ascii", ensight_format="gold", setup="2D", tensor_type="symmetric"):
    """
    Tests EnSight tensor (symmetric and asymmetric) variable output for nodes and elements.

    Runs the EnSight output process with a Matrix variable (CAUCHY_STRESS_TENSOR for symmetric,
    DEFORMATION_GRADIENT for asymmetric) on elements and verifies that the output files are
    created and non-empty. This covers the 32-bit float binary casting path as well as the
    ASCII path for both EnSight 6 and EnSight Gold formats.

    Args:
        file_format (str): "ascii" or "binary".
        ensight_format (str): "5", "6", or "gold".
        setup (str): Geometry setup — "2D" or "3D".
        tensor_type (str): "symmetric" or "asymmetric".
    """
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    current_model = KratosMultiphysics.Model()
    model_part_name = "Main"
    model_part = current_model.CreateModelPart(model_part_name)

    if setup == "2D":
        SetupModelPart2D(model_part)
    elif setup == "3D":
        SetupModelPart3D(model_part)
    else:
        raise ValueError(f"Unknown setup: {setup}")

    tensor_variable = "CAUCHY_STRESS_TENSOR" if tensor_type == "symmetric" else "DEFORMATION_GRADIENT"

    ensight_output_parameters = KratosMultiphysics.Parameters(f"""{{
        "Parameters" : {{
            "model_part_name"                    : "{model_part_name}",
            "ensight_file_format"                : "{ensight_format}",
            "file_format"                        : "{file_format}",
            "output_precision"                   : 6,
            "output_interval"                    : 1,
            "output_sub_model_parts"             : false,
            "evolving_geometry"                  : true,
            "output_path"                        : "test_ensight_output",
            "nodal_solution_step_data_variables" : ["PRESSURE", "DISPLACEMENT"],
            "element_data_value_variables"       : ["{tensor_variable}"]
        }}
    }}""")

    ensight_proc = SetupEnsightOutputProcess(current_model, ensight_output_parameters)

    time = 0.0
    dt = 1.0
    end_time = 2.0
    ensight_proc.ExecuteInitialize()

    step = 0
    while time < end_time:
        time += dt
        model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        SetSolutionWithMatrixVariables(model_part, tensor_type)
        ensight_proc.ExecuteInitializeSolutionStep()
        model_part.CloneTimeStep(time)
        ensight_proc.ExecuteFinalizeSolutionStep()
        if ensight_proc.IsOutputStep():
            ensight_proc.PrintOutput()

            label_step = "000" + str(step)

            # Verify geometry file was produced
            CheckFileNonEmpty(os.path.join("test_ensight_output", f"Main.{label_step}.geo"))

            # Verify tensor variable file was produced (.ten extension)
            CheckFileNonEmpty(os.path.join(
                "test_ensight_output",
                f"Main.{tensor_variable}.{label_step}.elem.ten"
            ))

            step += 1

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
