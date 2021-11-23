import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM
import KratosMultiphysics.print_info_in_file_process as PrintProcess
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestPrintInfoInFile(KratosUnittest.TestCase):

    current_path = os.path.abspath(os.getcwd())

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

    node1.SetSolutionStepValue(KM.DISPLACEMENT, [1.1,0.0,0.0])
    node2.SetSolutionStepValue(KM.DISPLACEMENT, [2.6,0.0,0.0])
    node3.SetSolutionStepValue(KM.DISPLACEMENT, [3.0,0.1,0.0])
    node4.SetSolutionStepValue(KM.DISPLACEMENT, [4.3,0.0,0.3])

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
    process.PrintOutput()

    ref_file_name = os.path.abspath(settings["file_name"].GetString())

    if not os.path.isfile(ref_file_name):
                err_msg  = 'The specified reference file name "'
                err_msg += ref_file_name
                err_msg += '" is not valid!'
                raise Exception(err_msg)

    with open(ref_file_name, "r") as plot_file:
        contents = plot_file.readlines()
        for line in contents:
            if line[0] != "#":
                numbers = line.split("\t")
                if len(numbers) > 2:
                    if float(numbers[0]) != 0.0 or float(numbers[1]) != 3.7 or float(numbers[2]) != 0.0 or float(numbers[3]) != 0.0:
                        raise Exception("The print does not give the expected result...")

    # now we test the second submodel
    settings["model_part_name"].SetString("test.to_plot_2")
    process = PrintProcess.PrintInfoInFileProcess(current_model, settings)
    process.PrintOutput()

    with open(ref_file_name, "r") as plot_file:
        contents = plot_file.readlines()
        for line in contents:
            if line[0] != "#":
                numbers = line.split("\t")
                if len(numbers) > 2:
                    if float(numbers[0]) != 0.0 or float(numbers[1]) != 7.3 or float(numbers[2]) != 0.1 or float(numbers[3]) != 0.3:
                        raise Exception("The print does not give the expected result...")
    process.ExecuteFinalize()
    os.remove(ref_file_name)

if __name__ == '__main__':
    KratosUnittest.main()
