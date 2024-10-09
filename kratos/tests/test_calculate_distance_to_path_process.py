import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

# Analytic distance function
def DistanceFunction(node, radius):
    x = node.X
    y = node.Y
    node_radius = (x**2 + y**2)**0.5 - 1.0
    z = node.Z
    sub_radius = (node_radius**2 + z**2)**0.5
    return sub_radius - radius

def VTKDebug(model):
    import KratosMultiphysics.vtk_output_process as vtk_output_process
    vtk_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : "Torus",
            "file_format"                        : "ascii",
            "output_precision"                   : 8,
            "output_interval"                    : 2,
            "output_sub_model_parts"             : false,
            "output_path"                        : "vtk_output_torus",
            "nodal_solution_step_data_variables" : ["DISTANCE"]
        }
    }""")

    vtk_output_process_torus = vtk_output_process.Factory(vtk_output_parameters, model)
    vtk_output_process_torus.ExecuteInitializeSolutionStep()
    vtk_output_process_torus.ExecuteFinalizeSolutionStep()
    vtk_output_process_torus.PrintOutput()

    vtk_output_parameters = KratosMultiphysics.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : "Circle",
            "file_format"                        : "ascii",
            "output_precision"                   : 8,
            "output_interval"                    : 2,
            "output_sub_model_parts"             : false,
            "output_path"                        : "vtk_output_circle",
            "nodal_solution_step_data_variables" : []
        }
    }""")

    vtk_output_process_circle = vtk_output_process.Factory(vtk_output_parameters, model)
    vtk_output_process_circle.ExecuteInitializeSolutionStep()
    vtk_output_process_circle.ExecuteFinalizeSolutionStep()
    vtk_output_process_circle.PrintOutput()

def GidDebug(model):
    from KratosMultiphysics.gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(
        model.GetModelPart("Torus"),
        "Torus",
        KratosMultiphysics.Parameters("""{
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostBinary",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "file_label": "time",
                "output_control_type": "step",
                "output_interval": 1.0,
                "body_output": true,
                "nodal_results": ["DISTANCE"]
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

    gid_output = GiDOutputProcess(
        model.GetModelPart("Circle"),
        "Circle",
        KratosMultiphysics.Parameters("""{
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostBinary",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "file_label": "time",
                "output_control_type": "step",
                "output_interval": 1.0,
                "body_output": true,
                "nodal_results": []
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

class TestCalculateDistanceToPathProcess(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def __base_test_calculate_distance_to_path_process(self, radius = 0.0, brute_force_calculation = False, check_tolerance = 5e-3, debug = False):
        # Define model
        self.current_model = KratosMultiphysics.Model()

        # Import torus
        model_part_torus = self.current_model.CreateModelPart("Torus")
        model_part_torus.ProcessInfo[KratosMultiphysics.STEP] = 0
        model_part_torus.CloneTimeStep(0.0)
        model_part_torus.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        input_mdpa = GetFilePath("test_files/mdpa_files/torus3d")
        model_part_io_torus = KratosMultiphysics.ModelPartIO(input_mdpa)
        model_part_io_torus.ReadModelPart(model_part_torus)

        # Import circle
        model_part_circle = self.current_model.CreateModelPart("Circle")
        model_part_circle.ProcessInfo[KratosMultiphysics.STEP] = 0
        model_part_circle.CloneTimeStep(0.0)
        input_mdpa = GetFilePath("test_files/mdpa_files/circle1d")
        model_part_io_circle = KratosMultiphysics.ModelPartIO(input_mdpa)
        model_part_io_circle.ReadModelPart(model_part_circle)

        # Set the distance function
        # ## Analytic distance function
        # for node in model_part_torus.Nodes:
        #     node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, DistanceFunction(node, radius))

        ## Compute distance
        parameters = KratosMultiphysics.Parameters("""{
            "distance_model_part_name" : "Torus",
            "path_model_part_name"     : "Circle",
            "distance_variable_name"   : "DISTANCE",
            "brute_force_calculation"  : false,
            "radius_path"              : 0.0,
            "distance_tolerance"       : 1.0e-9
        }""")

        parameters["radius_path"].SetDouble(radius)
        parameters["brute_force_calculation"].SetBool(brute_force_calculation)
        calculate_distance = KratosMultiphysics.CalculateDistanceToPathProcess(self.current_model, parameters)
        calculate_distance.Execute()

        # Check results
        for node in model_part_torus.Nodes:
            if debug:
                print(node.Id, "\t", node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), "\t", DistanceFunction(node, radius))
            else:
                self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), DistanceFunction(node, radius), delta=check_tolerance)

        # Debug
        if debug:
            VTKDebug(self.current_model)
            GidDebug(self.current_model)

    def test_calculate_distance_to_path_process_brute_force_zero_radius(self):
        self.__base_test_calculate_distance_to_path_process(0.0, True, 5e-3, False)

    def test_calculate_distance_to_path_process_brute_force_radius(self):
        self.__base_test_calculate_distance_to_path_process(0.1, True, 5e-3, False)

    def test_calculate_distance_to_path_process_zero_radius(self):
        self.__base_test_calculate_distance_to_path_process(0.0, False, 5e-3, False)

    def test_calculate_distance_to_path_process_radius(self):
        self.__base_test_calculate_distance_to_path_process(0.1, False, 5e-3, False)


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
