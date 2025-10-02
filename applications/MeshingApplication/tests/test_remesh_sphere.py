# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

# Import stuff
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def VTKDebug(model_part):
    """
    Exports the given Kratos model part to VTK format for debugging purposes.

    This function initializes a VTK output process for the provided model part,
    configures the output settings (such as file mode, mesh writing options, and nodal results),
    and executes the necessary steps to write the output files. The output includes
    the nodal "DISTANCE" variable and is written in binary format as a single file.
    """
    from KratosMultiphysics.vtk_output_process import VtkOutputProcess
    vtk_output = VtkOutputProcess(model_part.GetModel(),
                                KratosMultiphysics.Parameters("""{
                                        "model_part_name"                    : "MainModelPart",
                                        "nodal_solution_step_data_variables" : ["DISTANCE", "DISTANCE_GRADIENT"]
                                    }
                                    """)
                                )

    vtk_output.ExecuteInitialize()
    vtk_output.ExecuteBeforeSolutionLoop()
    vtk_output.ExecuteInitializeSolutionStep()
    vtk_output.PrintOutput()
    vtk_output.ExecuteFinalizeSolutionStep()
    vtk_output.ExecuteFinalize()

def GiDDebug(model_part):
    """
    Exports the given Kratos model part to GiD format for debugging purposes.

    This function initializes a GiD output process for the provided model part,
    configures the output settings (such as file mode, mesh writing options, and nodal results),
    and executes the necessary steps to write the output files. The output includes
    the nodal "DISTANCE" variable and is written in binary format as a single file.
    """
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(main_model_part,
                                "gid_output",
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "nodal_results"       : ["DISTANCE", "DISTANCE_GRADIENT"]
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

def GenerateSolution(current_model, output_file_name):
    """
    Generates a solution database by exporting specified variables from the given model to a JSON file.

    This function initializes and executes a JSON output process for the provided Kratos model,
    saving the specified variables (e.g., "DISTANCE") to the given output file.
    """
    # The following is used to create the solution database
    from KratosMultiphysics import json_output_process

    out_parameters = KratosMultiphysics.Parameters("""
    {
        "output_variables"     : ["DISTANCE"],
        "output_file_name"     : "",
        "model_part_name"      : "MainModelPart",
        "time_frequency"       : 0.0
    }
    """)
    out_parameters["output_file_name"].SetString(output_file_name)

    out = json_output_process.JsonOutputProcess(current_model, out_parameters)
    out.ExecuteInitialize()
    out.ExecuteBeforeSolutionLoop()
    out.ExecuteFinalizeSolutionStep()

def CheckSolution(reference_file_name, output_file_name, dimension = 3):
    """
    Checks the solution by comparing the reference and output files.
    """
    check_parameters = KratosMultiphysics.Parameters("""
    {
        "reference_file_name"   : "",
        "output_file_name"      : "",
        "dimension"             : 3,
        "comparison_type"       : "sol_file"
    }
    """)
    check_parameters["reference_file_name"].SetString(reference_file_name)
    check_parameters["output_file_name"].SetString(output_file_name)
    check_parameters["dimension"].SetInt(dimension)
    check_files = CompareTwoFilesCheckProcess(check_parameters)

    check_files.ExecuteInitialize()
    check_files.ExecuteBeforeSolutionLoop()
    check_files.ExecuteInitializeSolutionStep()
    check_files.ExecuteFinalizeSolutionStep()
    check_files.ExecuteFinalize()

def CheckDistanceExtrapolation(current_model, input_filename, mmg_version):
    check_parameters = KratosMultiphysics.Parameters("""
    {
        "check_variables"      : ["DISTANCE"],
        "input_file_name"      : "",
        "model_part_name"      : "MainModelPart",
        "time_frequency"       : 0.0
    }
    """)

    input_file_path = GetFilePath(input_filename)
    check_parameters["input_file_name"].SetString(input_file_path)
    if mmg_version == VersionToNumber("5.8.0"):
        check_parameters["input_file_name"].SetString(input_file_path + "_5_8.json")
    elif mmg_version == VersionToNumber("5.7.0"):
        check_parameters["input_file_name"].SetString(input_file_path + "_5_7.json")
    elif mmg_version >= VersionToNumber("5.5.0") and mmg_version < VersionToNumber("5.7.0"):
        check_parameters["input_file_name"].SetString(input_file_path + "_5_5.json")
    else:
        check_parameters["input_file_name"].SetString(input_file_path + ".json")
    check = FromJsonCheckResultProcess(current_model, check_parameters)
    check.ExecuteInitialize()
    check.ExecuteBeforeSolutionLoop()
    check.ExecuteFinalizeSolutionStep()

def VersionToNumber(version):
    """
    Converts a version string to a numerical representation for easy comparison.

    This function takes a version string formatted as "X.Y.Z" and converts it into
    a single integer by multiplying the major version (X) by 10000, the minor version (Y)
    by 100, and adding the patch version (Z). This allows for straightforward numerical
    comparisons between different version strings.

    Args:
        version (str): The version string to convert, formatted as "X.Y.Z".

    Returns:
        int: The numerical representation of the version.
    """
    parts = version.split('.')
    major = int(parts[0]) if len(parts) > 0 else 0
    minor = int(parts[1]) if len(parts) > 1 else 0
    patch = int(parts[2]) if len(parts) > 2 else 0
    return major * 10000 + minor * 100 + patch

class TestRemeshMMG3D(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def test_remesh_sphere(self):
        """
        Tests the remeshing of a spherical geometry using the MMG library within Kratos Multiphysics.

        This test performs the following steps:
            1. Imports a coarse sphere model part from file.
            2. Computes the nodal H (characteristic length) and the gradient of the DISTANCE variable.
            3. Initializes the metric tensor for all nodes to zero.
            4. Defines and applies a metric tensor based on the level set solution using ComputeLevelSetSolMetricProcess.
            5. Configures and executes the MMG remeshing process.
            6. Optionally exports the result to GiD (commented out).
            7. Checks the remeshing results by comparing the solution files.
            8. Performs additional checks using FromJsonCheckResultProcess to validate the DISTANCE variable against reference data.

        Raises:
            AssertionError: If the remeshing or solution checks fail.
        """
        # We fill the model
        self._FillModel("mmg_eulerian_test/coarse_sphere_test")

        # We calculate the gradient of the distance variable
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.main_model_part.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        for i in range(6):
            ZeroVector[i] = 0.0

        for node in self.main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

        # We define a metric using the ComputeLevelSetSolMetricProcess
        MetricParameters = KratosMultiphysics.Parameters("""
        {
            "minimal_size"                      : 1.0e-1,
            "enforce_current"                   : false,
            "anisotropy_remeshing"              : false,
            "anisotropy_parameters"             :{
                "hmin_over_hmax_anisotropic_ratio"  : 0.15,
                "boundary_layer_max_distance"       : 1.0e-4,
                "interpolation"                     : "Linear"
            }
        }
        """)
        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(self.main_model_part, KratosMultiphysics.DISTANCE_GRADIENT, MetricParameters)
        metric_process.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_eulerian_test/coarse_sphere_test",
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(GetFilePath(mmg_parameters["filename"].GetString()))
        mmg_process = MeshingApplication.MmgProcess3D(self.main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        # # Finally we export to GiD (for debugging purposes)
        # GiDDebug(self.main_model_part)

        # # Export to VTK (for debugging purposes)
        # VTKDebug(self.main_model_part)

        # We check the results
        CheckSolution(GetFilePath("mmg_eulerian_test/coarse_sphere_test_result.sol"),
                      GetFilePath("mmg_eulerian_test/coarse_sphere_test_step=0.sol"))
        CheckDistanceExtrapolation(self.current_model, "mmg_eulerian_test/distante_extrapolation", VersionToNumber(mmg_process.GetMmgVersion()))

        ## The following is used to create the solution database
        # GenerateSolution(self.current_model, "mmg_eulerian_test/distante_extrapolation.json")

    def test_remesh_sphere_skin(self):
        """
        Tests the remeshing of a sphere's skin using the MMG surface remeshing process.

        This test performs the following steps:
            1. Loads a coarse sphere skin model part from file.
            2. Initializes the DISTANCE variable for each node based on the absolute X coordinate.
            3. Calculates the nodal characteristic length (NODAL_H) for all nodes.
            4. Sets a predefined metric tensor for all nodes to guide the remeshing.
            5. Configures and executes the MMG surface remeshing process.
            6. Compares the remeshed solution with a reference solution file.
            7. Performs additional checks on the DISTANCE variable using a JSON-based result checker.
            8. Optionally generates a solution database for further validation.
        """
        # We fill the model
        self._FillModel("mmg_eulerian_test/coarse_sphere_skin_test")

        # We set manually a distance function
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, abs(node.X))

        # We calculate the gradient of the distance variable
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, self.main_model_part.Nodes)
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()

        # We set to zero the metric
        metric_vector = KratosMultiphysics.Vector(6)
        for i in range(3):
            metric_vector[i] = 1.0
        for i in range(3, 6):
            metric_vector[i] = 0.0

        for node in self.main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_3D, metric_vector)

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_eulerian_test/coarse_sphere_skin_test",
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(GetFilePath(mmg_parameters["filename"].GetString()))
        mmg_process = MeshingApplication.MmgProcess3DSurfaces(self.main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        # # Finally we export to GiD
        # GiDDebug(self.main_model_part)

        # # Export to VTK (for debugging purposes)
        # VTKDebug(self.main_model_part)

        # We check the results
        CheckSolution(GetFilePath("mmg_eulerian_test/coarse_sphere_skin_test_result.sol"),
                      GetFilePath("mmg_eulerian_test/coarse_sphere_skin_test_step=0.sol"))
        CheckDistanceExtrapolation(self.current_model, "mmg_eulerian_test/distante_extrapolation_skin", VersionToNumber(mmg_process.GetMmgVersion()))

        ## The following is used to create the solution database
        # GenerateSolution(self.current_model, "mmg_eulerian_test/distante_extrapolation_skin.json")

    def test_remesh_sphere_skin_prisms(self):
        """
        Tests the remeshing of a sphere's skin with prism elements using the MMG surface remeshing process.

        This test performs the following steps:
            1. Loads a coarse sphere model part with skin prisms from file.
            2. Initializes the DISTANCE variable for each node as the absolute value of its X coordinate.
            3. Sets the METRIC_SCALAR variable to zero for all nodes.
            4. Calculates the non-historical nodal H (characteristic length) for all nodes.
            5. Sets a metric tensor for each node, with the first three components set to 1.0 and the last three to 0.0.
            6. Configures MMG remeshing parameters, enabling prism element collapse and external file saving.
            7. Executes the MMG surface remeshing process.
            8. Optionally exports the result to GiD for debugging (commented out).
            9. Compares the remeshed solution to a reference solution to validate correctness.
        """
        # We fill the model
        self._FillModel("mmg_eulerian_test/coarse_sphere_skin_prisms_test")

        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, abs(node.X))
            node.SetValue(MeshingApplication.METRIC_SCALAR, 0.0)

        # We calculate the gradient of the distance variable
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, self.main_model_part.Nodes)
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()

        # We set to zero the metric
        metric_vector = KratosMultiphysics.Vector(6)
        for i in range(3):
            metric_vector[i] = 1.0
        for i in range(3, 6):
            metric_vector[i] = 0.0

        for node in self.main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_3D, metric_vector)

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_eulerian_test/coarse_sphere_skin_prisms_test",
            "collapse_prisms_elements"         : true,
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(GetFilePath(mmg_parameters["filename"].GetString()))
        mmg_process = MeshingApplication.MmgProcess3DSurfaces(self.main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        ## Finally we export to GiD
        # GiDDebug(self.main_model_part)

        # We check the results
        CheckSolution(GetFilePath("mmg_eulerian_test/coarse_sphere_skin_prisms_test_result.sol"),
                      GetFilePath("mmg_eulerian_test/coarse_sphere_skin_prisms_test_step=0.sol"))

    def test_isosurface_remesh_sphere(self):
        """
        Tests the isosurface remeshing of a sphere using the MMG library in Kratos Multiphysics.

        This test performs the following steps:
            1. Loads a model part from a specified file.
            2. Manually sets a signed distance function on the nodes to represent a sphere.
            3. Calculates the nodal characteristic length (NodalH) for the mesh.
            4. Configures MMG remeshing parameters for isosurface discretization.
            5. Executes the MMG remeshing process.
            6. Optionally exports the result to GiD for debugging (commented out).
            7. Checks the remeshed solution against a reference solution.

        The test validates that the isosurface remeshing process produces the expected mesh for a sphere defined by the distance function.
        """
        # We fill the model
        self._FillModel("mmg_eulerian_test/test_sphere_isosurface")

        # Set manually a distance function
        circle_radious = 0.25
        center_coordinates = [0.5, 0.5, 0.5]

        for node in self.main_model_part.Nodes:
            distance = ((node.X-center_coordinates[0])**2+(node.Y-center_coordinates[1])**2+(node.Z-center_coordinates[2])**2)**0.5 - circle_radious
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

        # We calculate the gradient of the distance variable
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"              : "Isosurface",
            "filename"                         : "mmg_eulerian_test/test_sphere_isosurface",
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(GetFilePath(mmg_parameters["filename"].GetString()))
        mmg_process = MeshingApplication.MmgProcess3D(self.main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        # # Finally we export to GiD
        # GiDDebug(self.main_model_part)

        # # Export to VTK (for debugging purposes)
        # VTKDebug(self.main_model_part)

        # We check the results
        if VersionToNumber(mmg_process.GetMmgVersion()) >= VersionToNumber("5.7.0"):
            CheckSolution(GetFilePath("mmg_eulerian_test/test_sphere_isosurface_result_5_7.sol"),
                          GetFilePath("mmg_eulerian_test/test_sphere_isosurface_step=0.sol"))
        else:
            CheckSolution(GetFilePath("mmg_eulerian_test/test_sphere_isosurface_result.sol"),
                          GetFilePath("mmg_eulerian_test/test_sphere_isosurface_step=0.sol"))

    def _FillModel(self, mdpa_name):
        """
        Initializes the main model part for the remeshing test.

        This method creates a new KratosMultiphysics model, sets up the main model part,
        and configures essential process information such as domain size, time, and delta time.
        It also adds the required nodal solution step variables for distance and distance gradient.
        """
        # We create the model part
        self.current_model = KratosMultiphysics.Model()
        self.main_model_part = self.current_model.CreateModelPart("MainModelPart")
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)

        # We add the variables needed
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        # We import the model self.main_model_part
        KratosMultiphysics.ModelPartIO(GetFilePath(mdpa_name)).ReadModelPart(self.main_model_part)

if __name__ == '__main__':
    # Configure logging level and start the test runner
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
