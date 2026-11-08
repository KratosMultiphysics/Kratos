import os
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
if kratos_utilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    from  KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis



@KratosUnittest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
class TestNumpyOutputProcess(KratosUnittest.TestCase):

    def testNumpyOutputProcess(self):
        self.work_folder = "thermal_dynamic_test_files/FOM/"
        parameters_filename = "../ProjectParameters.json"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                self.full_parameters = KratosMultiphysics.Parameters(parameter_file.read())

            process_settings = KratosMultiphysics.Parameters("""{
                "Parameters": {
                    "model_part_name": "ThermalModelPart",
                    "output_path": "../../numpy_output_process_test_files/numpy_output",
                    "nodal_results": ["TEMPERATURE"],
                    "output_control_type": "step",
                    "output_interval": 40.0
                },
                "kratos_module": "KratosMultiphysics.RomApplication",
                "process_name": "NumpyOutputProcess",
                "python_module": "numpy_output_process"
            }""")

            self.full_parameters["output_processes"]["numpy_output"].Append(process_settings)
            model = KratosMultiphysics.Model()
            self.simulation = ConvectionDiffusionAnalysis(model, self.full_parameters)

            # test step-based outputs
            self.simulation.Run()
            self.__CheckResults_step()

            model = KratosMultiphysics.Model()
            self.full_parameters["output_processes"]["numpy_output"][0]["Parameters"]["output_control_type"].SetString("time")
            self.full_parameters["output_processes"]["numpy_output"][0]["Parameters"]["output_interval"].SetDouble(1000.0)
            self.simulation = ConvectionDiffusionAnalysis(model, self.full_parameters)

            # test time-based outputs
            self.simulation.Run()
            self.__CheckResults_time()


    def __CheckResults_step(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):

            output_interval = self.full_parameters["output_processes"]["numpy_output"][0]["Parameters"]["output_interval"].GetDouble()
            end_time = self.full_parameters["problem_data"]["end_time"].GetDouble()
            time_stepping = self.full_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()

            number_of_steps = int(end_time/(time_stepping*output_interval))
            for i in range(1,number_of_steps+1):
                this_step = int(i*output_interval)
                # Load outputed solution
                self.output_folder = self.full_parameters["output_processes"]["numpy_output"][0]["Parameters"]["output_path"].GetString()
                outputed_solution = os.path.join(self.output_folder, f"solution_{this_step}.npy")
                rom_coeffs_obtained_output = np.load(outputed_solution)

                # Load reference solution
                reference_output_name = f"ExpectedOutput_numpy_output_step/solution_{this_step}.npy"
                reference_output = np.load(reference_output_name)

                self.assertMatrixAlmostEqual(KratosMultiphysics.Matrix(rom_coeffs_obtained_output), KratosMultiphysics.Matrix(reference_output))


    def __CheckResults_time(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Check results
            output_interval = self.full_parameters["output_processes"]["numpy_output"][0]["Parameters"]["output_interval"].GetDouble()
            end_time = self.full_parameters["problem_data"]["end_time"].GetDouble()
            number_of_outputs = int(end_time / output_interval)

            for i in range(1, number_of_outputs + 1):
                this_time = float(i * output_interval)

                # Format the time string to have decimal places as per the numpy output process
                time_label = f"{this_time:.7f}"

                # Load outputed solution
                self.output_folder = self.full_parameters["output_processes"]["numpy_output"][0]["Parameters"]["output_path"].GetString()
                outputed_solution = os.path.join(self.output_folder, f"solution_{time_label}.npy")
                rom_coeffs_obtained_output = np.load(outputed_solution)

                # Load reference solution
                reference_output_name = f"ExpectedOutput_numpy_output_time/solution_{time_label}.npy"
                reference_output = np.load(reference_output_name)

                self.assertMatrixAlmostEqual(KratosMultiphysics.Matrix(rom_coeffs_obtained_output), KratosMultiphysics.Matrix(reference_output))


    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            output_directory_path = os.path.join(os.getcwd(), self.output_folder)
            kratos_utilities.DeleteDirectoryIfExisting(output_directory_path)


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
