import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM
import KratosMultiphysics.print_info_in_file_process as PrintProcess
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestPrintInfoInFile(KratosUnittest.TestCase):

    current_model = KM.Model()
    model_part = current_model.CreateModelPart("test")
    properties = KM.Properties(1)
    model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

    model_part.ProcessInfo = KM.ProcessInfo()
    model_part.ProcessInfo[KM.STEP] = 1

    node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
    node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
    nodes = [node1, node2]

    geom = KM.Tetrahedra3D4(node1,node2,node3,node4)

    submodel_1 = model_part.CreateSubModelPart("to_plot")
    submodel_1.AddNodes([1,2])

    submodel_2 = model_part.CreateSubModelPart("to_plot_2")
    submodel_2.AddNodes([3,4])

    node1.SetSolutionStepValue(KM.DISPLACEMENT, [1.0,0,0])
    node2.SetSolutionStepValue(KM.DISPLACEMENT, [2.0,0,0])
    node3.SetSolutionStepValue(KM.DISPLACEMENT, [3.0,0,0])
    node4.SetSolutionStepValue(KM.DISPLACEMENT, [4.0,0,0])

    settings = KM.Parameters("""{
        "model_part_name"                    : "test.to_plot",
        "variable_name"                      : "DISPLACEMENT",
        "results_type"                       : "nodal_historical",
        "file_name"                          : "info_file.txt",
        "output_control_type"                : "step",
        "erase_previous_info"                : true,
        "output_interval"                    : 1,
        "sum_results_from_multiple_entities" : true,
        "write_buffer_size"                  : 1,
        "output_path"                        : ""}""")
    process = PrintProcess.PrintInfoInFileProcess(current_model, settings)

























if __name__ == '__main__':
    KratosUnittest.main()
