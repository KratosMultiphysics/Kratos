from __future__ import print_function, absolute_import, division

import os

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("FluidDynamicsApplication")
if dependencies_are_available:
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestSerializer(KratosUnittest.TestCase):
    def _prepare_fluid_test(self):
        # Define a model and load the parameters
        self.pre_serialized_model = KratosMultiphysics.Model()
        with open(GetFilePath("auxiliar_files_for_python_unittest/parameters_files/test_serializer.json"),'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        file_name = parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(GetFilePath(file_name))
        # First the model is initialized
        self.pre_serialized_simulation = FluidDynamicsAnalysis(self.pre_serialized_model, parameters)
        self.pre_serialized_simulation.Initialize()

        # Before serializing the model, main model part is set to RESTARTED
        self.main_model_part_name = parameters["solver_settings"]["model_part_name"].GetString()
        self.pre_serialized_model.GetModelPart(self.main_model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,True)
        self.serialized_model = KratosMultiphysics.StreamSerializer()
        self.serialized_model.Save("ModelSerialization",self.pre_serialized_model)

        with open(GetFilePath("auxiliar_files_for_python_unittest/parameters_files/test_serializer.json"),'r') as parameter_file:
            self.project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        # Parameters are read again and input type set to use_input_model_part since the serialized model already has the mdpa loaded
        self.project_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString("use_input_model_part")

        # Deserialize and store the new model
        self.current_model = KratosMultiphysics.Model()
        self.serialized_model.Load("ModelSerialization",self.current_model)

    def _check_results(self):
        pre_serialized_model_part = self.pre_serialized_model.GetModelPart(self.main_model_part_name)
        pre_serialized_pressure_results = [node.GetSolutionStepValue(KratosMultiphysics.PRESSURE) for node in pre_serialized_model_part.Nodes]
        pre_serialized_velocity_results = [node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) for node in pre_serialized_model_part.Nodes]

        serialized_model_part = self.current_model.GetModelPart(self.main_model_part_name)
        serialized_pressure_results = [node.GetSolutionStepValue(KratosMultiphysics.PRESSURE) for node in serialized_model_part.Nodes]
        serialized_velocity_results = [node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) for node in serialized_model_part.Nodes]

        # Comparing results before and after serializing
        for pre_serialized_result, serialized_result in zip(pre_serialized_pressure_results,serialized_pressure_results):
            self.assertAlmostEqual(pre_serialized_result, serialized_result)
        for pre_serialized_result, serialized_result in zip(pre_serialized_velocity_results,serialized_velocity_results):
            for value_pre_seralized, value_serialized in zip(pre_serialized_result, serialized_result):
                self.assertAlmostEqual(value_pre_seralized, value_serialized)

    @KratosUnittest.skipUnless(dependencies_are_available,"FluidDynamicsApplication is not available")
    def test_serializer_fluid_analysis(self):
        self._prepare_fluid_test()
        # Solving simulation before serializing to later check the results
        self.pre_serialized_simulation.RunSolutionLoop()
        self.pre_serialized_simulation.Finalize()
        # Solving simulation after serializing
        self.serialized_simulation = FluidDynamicsAnalysis(self.current_model, self.project_parameters)
        self.serialized_simulation.Run()
        self._check_results()

    def test_serializer_loading(self):
        # Creating a model with nodes and variables
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateNewNode(1,0.0,0.0,0.0)
        other = current_model.CreateModelPart("Other")
        other.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        other.CreateNewNode(1,0.0,0.0,0.0)

        # Serializing model
        serialized_model = KratosMultiphysics.StreamSerializer()
        serialized_model.Save("ModelSerialization", current_model)

        # Loading model several times
        first_model = KratosMultiphysics.Model()
        serialized_model.Load("ModelSerialization", first_model)
        second_model = KratosMultiphysics.Model()
        serialized_model.LoadFromBeginning("ModelSerialization", second_model)
        third_model = KratosMultiphysics.Model()
        serialized_model.LoadFromBeginning("ModelSerialization", third_model)

        self.assertTrue(third_model["Main"].HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertTrue(third_model["Other"].HasNodalSolutionStepVariable(KratosMultiphysics.PRESSURE))
        self.assertTrue(third_model.HasModelPart("Main.Inlets"))
        self.assertTrue(third_model.HasModelPart("Main.Temp"))
        self.assertTrue(1 in third_model["Main"].Nodes)
        self.assertTrue(1 in third_model["Other"].Nodes)

if __name__ == '__main__':
    KratosUnittest.main()
