# We import the libraries
import KratosMultiphysics
from KratosMultiphysics import TIME, DELTA_TIME
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.MeshingApplication import gradual_variable_interpolation_process
import numpy as np

import os

class TestGradualVariableInterpolationProcess(KratosUnittest.TestCase):

    def _SimulateStep(cls, modelpart, process):
        # Updating the time step by incrementing the DELTA_TIME
        new_t = modelpart.ProcessInfo[TIME] + modelpart.ProcessInfo[DELTA_TIME]
        # Cloning the time step to keep track of the time evolution
        modelpart.CloneTimeStep(new_t)
        # Running the interpolation process for each solution step
        process.ExecuteInitializeSolutionStep()

    def test_gradual_variable_interpolation_process(self):
        # Set logging level to warning to avoid cluttering the output
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # Create the model part
        current_model = KratosMultiphysics.Model()
        destination_model_part = current_model.CreateModelPart("destination_model_part")
        destination_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        # Import the destination model part from a file
        file_path = os.path.dirname(os.path.realpath(__file__))
        model_part_io = KratosMultiphysics.ModelPartIO(f"{file_path}/gradual_variable_interpolation_test_files/destination_model_part")
        model_part_io.ReadModelPart(destination_model_part)

        # Initialize the process info
        destination_model_part.ProcessInfo.SetValue(TIME, 0.0)
        destination_model_part.ProcessInfo.SetValue(DELTA_TIME, 1.0)

        # Define the settings for the interpolation process
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "destination_model_part_name": "destination_model_part",
                "origin_model_part_file_name": "gradual_variable_interpolation_test_files/origin_model_part",
                "interpolation_variables_list": ["TEMPERATURE"],
                "constrain_variables": true,
                "steps_for_rampup": 2,
                "alpha_rampup_increment": 0.0
            }
        }""")

        # Create the interpolation process
        process = gradual_variable_interpolation_process.Factory(settings, current_model)
        process.ExecuteInitialize()

        # Run the interpolation process for each time step and collect the interpolated values
        obtained_output = []
        for i in range(2):
            self._SimulateStep(destination_model_part, process)
            for node in destination_model_part.Nodes:
                obtained_output.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        
        # Load the expected output and compare it with the obtained output
        expected_output = np.load("gradual_variable_interpolation_test_files/ExpectedOutput.npy")
        self.assertVectorAlmostEqual(expected_output, obtained_output)
        
if __name__ == '__main__':
    KratosUnittest.main()
