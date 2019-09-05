from __future__ import print_function, absolute_import, division

import os
import sys

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("FluidDynamicsApplication")
if dependencies_are_available:
    import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis



def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestSerializer(KratosUnittest.TestCase):
    def _prepare_test(self):
        # Define a Model
        self.pre_serialized_model = KratosMultiphysics.Model()
        with open("auxiliar_files_for_python_unnitest/parameters_files/test_serializer.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        ## First the model is initialized
        self.pre_serialized_simulation = FluidDynamicsAnalysis(self.pre_serialized_model, parameters)
        self.pre_serialized_simulation.Initialize()

        ## Before serializing the model, main model part is set to RESTARTED
        self.main_model_part_name = parameters["solver_settings"]["model_part_name"].GetString()
        self.pre_serialized_model.GetModelPart(self.main_model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,True)
        serialized_model = KratosMultiphysics.StreamSerializer()
        serialized_model.Save("ModelSerialization",self.pre_serialized_model)

        with open("auxiliar_files_for_python_unnitest/parameters_files/test_serializer.json",'r') as parameter_file:
            self.project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        ## Parameters are read again and input type set to use_input_model_part since the serialized model already has the mdpa loaded
        self.project_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString("use_input_model_part")

        ## Deserialize and store the new model
        self.current_model = KratosMultiphysics.Model()
        serialized_model.Load("ModelSerialization",self.current_model)
        del(serialized_model)

    def _check_results(self):
        pre_serialized_model_part = self.pre_serialized_model.GetModelPart(self.main_model_part_name).GetSubModelPart("NoSlip2D_structure")
        pre_serialized_results = []
        for node in pre_serialized_model_part.Nodes:
            pressure = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
            pre_serialized_results.append(pressure)
        serialized_model_part = self.current_model.GetModelPart(self.main_model_part_name).GetSubModelPart("NoSlip2D_structure")
        serialized_results = []
        for node in serialized_model_part.Nodes:
            pressure = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
            serialized_results.append(pressure)
        for i_node in range(len(serialized_results)):
            self.assertAlmostEqual(pre_serialized_results[i_node], serialized_results[i_node])

    @KratosUnittest.skipUnless(dependencies_are_available,"FluidDynamicsApplication is not available")
    def test_serializer_fluid_analysis(self):
        self._prepare_test()
        # Solving simulation before serializing to later check the results
        self.pre_serialized_simulation.RunSolutionLoop()
        self.pre_serialized_simulation.Finalize()
        # Solving simulation after serializing
        self.serialized_simulation = FluidDynamicsAnalysis(self.current_model, self.project_parameters)
        self.serialized_simulation.Run()
        self._check_results()

if __name__ == '__main__':
    KratosUnittest.main()
