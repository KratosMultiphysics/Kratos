import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics import process_factory

import math
import os

try:
    import numpy
    missing_numpy = False
except ImportError:
    missing_numpy = True

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestLineGraphOutputProcess(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(missing_numpy,"Missing python libraries (numpy)")
    def test_line_graph_output_process(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_files/mdpa_files/test_processes"))
        model_part_io.ReadModelPart(model_part)

        reference_file_name_1 = GetFilePath("test_files/point_output_process_ref_files/line_graph_output_ref_1.dat")
        reference_file_name_2 = GetFilePath("test_files/point_output_process_ref_files/line_graph_output_ref_2.dat")

        settings = KratosMultiphysics.Parameters("""{
            "process_list" : [{
                "python_module"  : "line_graph_output_process",
                "kratos_module"  : "KratosMultiphysics.ShallowWaterApplication.postprocess",
                "process_name"   : "LineGraphOutputProcess",
                "Parameters"            : {
                    "model_part_name"      : "Main",
                    "start_point"          : [0.55, 0.2, 0.0],
                    "end_point"            : [0.65, 0.4, 0.0],
                    "output_variables"     : ["DISPLACEMENT", "VISCOSITY"],
                    "entity_type"          : "element",
                    "search_configuration" : "current",
                    "output_file_settings" : {
                        "file_name"   : "Main",
                        "output_path" : ""
                    },
                    "output_control_settings" : {
                        "output_control_type"     : "time",
                        "time_frequency"          : 3.0
                    }
                }
            },{
                "python_module"  : "compare_two_files_check_process",
                "kratos_module"  : "KratosMultiphysics",
                "process_name"   : "CompareTwoFilesCheckProcess",
                "Parameters"            : {
                    "reference_file_name"   : "",
                    "output_file_name"      : "Main_0.300.dat",
                    "comparison_type"       : "dat_file"
                }
            },{
                "python_module"  : "compare_two_files_check_process",
                "kratos_module"  : "KratosMultiphysics",
                "process_name"   : "CompareTwoFilesCheckProcess",
                "Parameters"            : {
                    "reference_file_name"   : "",
                    "output_file_name"      : "Main_3.300.dat",
                    "comparison_type"       : "dat_file"
                }
            }]
        }""")

        settings["process_list"][1]["Parameters"]["reference_file_name"].SetString(reference_file_name_1)
        settings["process_list"][2]["Parameters"]["reference_file_name"].SetString(reference_file_name_2)

        end_time = 5.0
        delta_time = 0.3

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0

        SolutionLoopPointOutputProcesses(model_part, settings, end_time, delta_time)


def SetNodalValuesForPointOutputProcesses(model_part):
    time = model_part.ProcessInfo[KratosMultiphysics.TIME]
    vec = KratosMultiphysics.Vector(3)
    for node in model_part.Nodes:
        vec[0] = round(math.sqrt(node.X**2+node.Y**2)*time ,6)
        vec[1] = round(node.X**2+node.Y**2 + time ,6)
        vec[2] = round(node.X+node.Y + time ,6)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, vec)
        node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION, vec*time)
        node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY, time**2 + 1.038)

def SolutionLoopPointOutputProcesses(model_part, settings, end_time, delta_time):
    current_model = model_part.GetModel()
    list_of_processes = process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses(
        settings["process_list"] )

    for process in list_of_processes:
        process.ExecuteInitialize()

    for process in list_of_processes:
        process.Check()

    for process in list_of_processes:
        process.ExecuteBeforeSolutionLoop()

    while model_part.ProcessInfo[KratosMultiphysics.TIME] < end_time:
        model_part.ProcessInfo[KratosMultiphysics.TIME] += delta_time

        SetNodalValuesForPointOutputProcesses(model_part)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        for process in list_of_processes:
            process.ExecuteBeforeOutputStep()

        for process in list_of_processes:
            if issubclass(type(process), KratosMultiphysics.OutputProcess):
                if process.IsOutputStep():
                    process.PrintOutput()

        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()
