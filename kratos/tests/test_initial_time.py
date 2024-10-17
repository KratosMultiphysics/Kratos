import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtils
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

if KratosUtils.CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


@KratosUnittest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
class TestInitialTime(KratosUnittest.TestCase):

    def testInitialTime(self):
        model = KratosMultiphysics.Model()
        parameters_name = GetFilePath("test_files/parameters_files/test_initial_time.json")
        with open(parameters_name, 'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self._SetModelPartName(parameters, "First")
        self._SetInitialTime(parameters, 0.0)

        simulation1 = ConvectionDiffusionAnalysis(model, parameters)
        simulation1.Run()

        self._SetModelPartName(parameters, "Second")
        self._SetInitialTime(parameters, 1000.0)

        simulation2 = ConvectionDiffusionAnalysis(model, parameters)
        simulation2.Run()

        model_part1 = model.GetModelPart("First")
        model_part2 = model.GetModelPart("Second")
        for node1, node2 in zip (model_part1.Nodes, model_part2.Nodes):
            temperature1 = node1.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            temperature2 = node2.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            self.assertAlmostEqual(temperature1, temperature2)


    @staticmethod
    def _SetModelPartName(parameters, model_part_name):
        parameters["solver_settings"]["model_part_name"].SetString(model_part_name)

        materials_pattern = "test_files/materials_files/initial_time_{}.json"
        materials_filename = GetFilePath(materials_pattern.format(model_part_name.lower()))
        parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString(materials_filename)

        for processes_list in parameters["processes"].values():
            for process in processes_list.values():
                process["Parameters"]["model_part_name"].SetString(model_part_name)

        if parameters.Has("output_processes"):
            for processes_list in parameters["output_processes"].values():
                for process in processes_list.values():
                    process["Parameters"]["model_part_name"].SetString(model_part_name)
                    process["Parameters"]["output_name"].SetString(model_part_name)

    @staticmethod
    def _SetInitialTime(parameters, initial_time):
        parameters["problem_data"]["start_time"].SetDouble(initial_time)
        parameters["problem_data"]["end_time"].SetDouble(initial_time + 1)

        for process in parameters["processes"]["initial_conditions_process_list"].values():
            process["Parameters"]["interval"].SetVector([0, initial_time])


if __name__ == '__main__':
    KratosUnittest.main()
