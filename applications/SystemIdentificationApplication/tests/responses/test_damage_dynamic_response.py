import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.StructuralMechanicsApplication
import math
import shutil, os
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
import KratosMultiphysics.kratos_utilities as kratos_utils

class TestDamageDetectionDynamicResponse(kratos_unittest.TestCase):
    def test_DamageDynamicResponse_OneSolidElement(self):
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_8/damaged_system/hdf5_output")
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_8/hdf5_output")
        self.addCleanup(kratos_utils.DeleteFileIfExisting, "auxiliary_files_8/one_element.time")

        with kratos_unittest.WorkFolderScope(".", __file__):
            
            # First run the damaged system to generate measured data
            with open("auxiliary_files_8/damaged_system/project_parameters.json", "r") as file_input1:
                damaged_system_parameters = Kratos.Parameters(file_input1.read())
            
            damaged_system_model = Kratos.Model()
            damaged_system_analysis = StructuralMechanicsAnalysis(damaged_system_model, damaged_system_parameters)
            damaged_system_analysis.Run()

            num_of_timesteps_ds = math.ceil((damaged_system_parameters["problem_data"]["end_time"].GetDouble() - damaged_system_parameters["problem_data"]["start_time"].GetDouble()) / damaged_system_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble())
            
            # Run the system identification
            with open("auxiliary_files_8/optimization_parameters_p_norm.json", "r") as file_input2:
                parameters = Kratos.Parameters(file_input2.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            num_of_timesteps_si = response.num_of_timesteps
            # Making sure damaged_system and system_identifcation have the same time steps and output files
            self.assertAlmostEqual(num_of_timesteps_ds, num_of_timesteps_si, 8)

            ref_value = response.CalculateValue()

            self.assertAlmostEqual(ref_value, 3.1199266806031143e-06, 9)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                Kratos.VariableUtils().ClearHistoricalData(model_part)
                E_orig = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig * (1+delta)
                sensitivity = ((response.CalculateValue() - ref_value) / delta)
                self.assertAlmostEqual(gradients[index]*E_orig, sensitivity/num_of_timesteps_si, 10)
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig

    def test_DamageDynamicResponse_CuboidSolid(self):
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_7/damaged_system/hdf5_output")
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_7/hdf5_output")
        self.addCleanup(kratos_utils.DeleteFileIfExisting, "auxiliary_files_7/Structure.time")

        with kratos_unittest.WorkFolderScope(".", __file__):
            
            # First run the damaged system to generate measured data
            with open("auxiliary_files_7/damaged_system/project_parameters.json", "r") as file_input1:
                damaged_system_parameters = Kratos.Parameters(file_input1.read())
            
            damaged_system_model = Kratos.Model()
            damaged_system_analysis = StructuralMechanicsAnalysis(damaged_system_model, damaged_system_parameters)
            damaged_system_analysis.Run()

            num_of_timesteps_ds = math.ceil((damaged_system_parameters["problem_data"]["end_time"].GetDouble() - damaged_system_parameters["problem_data"]["start_time"].GetDouble()) / damaged_system_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble())
            
            # Run the system identification
            with open("auxiliary_files_7/optimization_parameters_p_norm.json", "r") as file_input2:
                parameters = Kratos.Parameters(file_input2.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            num_of_timesteps_si = response.num_of_timesteps
            # Making sure damaged_system and system_identifcation have the same time steps and output files
            self.assertAlmostEqual(num_of_timesteps_ds, num_of_timesteps_si, 8)

            ref_value = response.CalculateValue()

            self.assertAlmostEqual(ref_value, 3.0560791665536015e-05, 9)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                Kratos.VariableUtils().ClearHistoricalData(model_part)
                E_orig = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig * (1+delta)
                sensitivity = ((response.CalculateValue() - ref_value) / delta)
                self.assertAlmostEqual(gradients[index]*E_orig, sensitivity/num_of_timesteps_si, 10)
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig

    def test_DamageDynamicResponse_OneShellElement(self):
        #self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_9/damaged_system/hdf5_output")
        #self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_9/hdf5_output")
        self.addCleanup(kratos_utils.DeleteFileIfExisting, "auxiliary_files_9/one_element.time")

        with kratos_unittest.WorkFolderScope(".", __file__):
            
            # First run the damaged system to generate measured data
            with open("auxiliary_files_9/damaged_system/project_parameters.json", "r") as file_input1:
                damaged_system_parameters = Kratos.Parameters(file_input1.read())
            
            damaged_system_model = Kratos.Model()
            damaged_system_analysis = StructuralMechanicsAnalysis(damaged_system_model, damaged_system_parameters)
            damaged_system_analysis.Run()

            num_of_timesteps_ds = math.ceil((damaged_system_parameters["problem_data"]["end_time"].GetDouble() - damaged_system_parameters["problem_data"]["start_time"].GetDouble()) / damaged_system_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble())
            
            # Run the system identification
            with open("auxiliary_files_9/optimization_parameters_p_norm.json", "r") as file_input2:
                parameters = Kratos.Parameters(file_input2.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            num_of_timesteps_si = response.num_of_timesteps
            # Making sure damaged_system and system_identifcation have the same time steps and output files
            self.assertAlmostEqual(num_of_timesteps_ds, num_of_timesteps_si, 8)

            ref_value = response.CalculateValue()

            value_outpt_file = "cost_function.dat"
            cost_function_ref = []

            with open(value_outpt_file, "r") as fi:
                line = fi.readline()
                while line:
                    line = fi.readline()
                    if line == "":
                        break
                    sline = line.strip()
                    cost_function_ref.append(float(sline))

            print(cost_function_ref)
            dir_name = os.path.dirname(value_outpt_file)
            shutil.copy(value_outpt_file, os.path.join(dir_name, f"auxiliary_files_9/cost_function_reference.dat"))

            #self.assertAlmostEqual(ref_value, 1.4744389948359893e-08, 10)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()

            delta = 1e-6
            for index, element in enumerate(model_part.Elements):
                Kratos.VariableUtils().ClearHistoricalData(model_part)
                E_orig = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig * (1+delta)
                val = response.CalculateValue()
                cost_function_curr = []
                with open(value_outpt_file, "r") as fi:
                    line = fi.readline()
                    while line:
                        line = fi.readline()
                        if line == "":
                            break
                        sline = line.strip()
                        cost_function_curr.append(float(sline))

                print("current cost function for pertrub element ", element.Id, " is ", cost_function_curr)
                shutil.copy(value_outpt_file, os.path.join(dir_name, f"auxiliary_files_9/cost_function_fd_elem_{element.Id}.dat"))
                #print("current val is ", val)

                corrected_cost_function_curr = list(cost_function_curr)
                for i in range(1, len(corrected_cost_function_curr)):
                    corrected_cost_function_curr[i] = cost_function_curr[i] - cost_function_curr[i-1]

                corrected_cost_function_ref = list(cost_function_ref)
                for i in range(1, len(corrected_cost_function_ref)):
                    corrected_cost_function_ref[i] = cost_function_ref[i] - cost_function_ref[i-1]

                print("cost_function_curr ", cost_function_curr)
                print("corrected_cost_function_curr ", corrected_cost_function_curr)
                print("cost_function_ref", cost_function_ref)
                print("corrected_cost_function_ref", corrected_cost_function_ref)

                sensitivity = ((np.array(corrected_cost_function_curr) - np.array(corrected_cost_function_ref)) / (delta))
                #sensitivity = ((response.CalculateValue() - ref_value) / delta)
                print("Element ", index+1 ," , Adjoint Sensitivity : " , gradients[index]*E_orig, " vs, FD : ", np.sum(sensitivity)/num_of_timesteps_si)
                self.assertAlmostEqual(gradients[index]*E_orig, np.sum(sensitivity)/num_of_timesteps_si, 9)
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig

    def test_DamageDynamicResponse_RectangleShell(self):
        #self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_10/damaged_system/hdf5_output")
        #self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "auxiliary_files_10/hdf5_output")
        self.addCleanup(kratos_utils.DeleteFileIfExisting, "auxiliary_files_10/one_element.time")

        with kratos_unittest.WorkFolderScope(".", __file__):
            
            # First run the damaged system to generate measured data
            with open("auxiliary_files_10/damaged_system/project_parameters.json", "r") as file_input1:
                damaged_system_parameters = Kratos.Parameters(file_input1.read())
            
            damaged_system_model = Kratos.Model()
            damaged_system_analysis = StructuralMechanicsAnalysis(damaged_system_model, damaged_system_parameters)
            damaged_system_analysis.Run()

            num_of_timesteps_ds = math.ceil((damaged_system_parameters["problem_data"]["end_time"].GetDouble() - damaged_system_parameters["problem_data"]["start_time"].GetDouble()) / damaged_system_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble())
            
            # Run the system identification
            with open("auxiliary_files_10/optimization_parameters_p_norm.json", "r") as file_input2:
                parameters = Kratos.Parameters(file_input2.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            num_of_timesteps_si = response.num_of_timesteps
            print("num_time_steps ", num_of_timesteps_ds)
            # Making sure damaged_system and system_identifcation have the same time steps and output files
            self.assertAlmostEqual(num_of_timesteps_ds, num_of_timesteps_si, 8)

            ref_value = response.CalculateValue()

            value_outpt_file = "cost_function.dat"
            cost_function_ref = []

            with open(value_outpt_file, "r") as fi:
                line = fi.readline()
                while line:
                    line = fi.readline()
                    if line == "":
                        break
                    sline = line.strip()
                    cost_function_ref.append(float(sline))

            print(cost_function_ref)
            dir_name = os.path.dirname(value_outpt_file)
            shutil.copy(value_outpt_file, os.path.join(dir_name, f"auxiliary_files_10/cost_function_reference.dat"))

            #self.assertAlmostEqual(ref_value, 2.034700695402209e-07, 10)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                Kratos.VariableUtils().ClearHistoricalData(model_part)
                E_orig = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig * (1+delta)

                val = response.CalculateValue()
                cost_function_curr = []
                with open(value_outpt_file, "r") as fi:
                    line = fi.readline()
                    while line:
                        line = fi.readline()
                        if line == "":
                            break
                        sline = line.strip()
                        cost_function_curr.append(float(sline))

                print("current cost function for pertrub element ", element.Id, " is ", cost_function_curr)
                shutil.copy(value_outpt_file, os.path.join(dir_name, f"auxiliary_files_10/cost_function_fd_elem_{element.Id}.dat"))

                corrected_cost_function_curr = list(cost_function_curr)
                for i in range(1, len(corrected_cost_function_curr)):
                    corrected_cost_function_curr[i] = cost_function_curr[i] - cost_function_curr[i-1]

                corrected_cost_function_ref = list(cost_function_ref)
                for i in range(1, len(corrected_cost_function_ref)):
                    corrected_cost_function_ref[i] = cost_function_ref[i] - cost_function_ref[i-1]

                print("cost_function_curr ", cost_function_curr)
                print("corrected_cost_function_curr ", corrected_cost_function_curr)
                print("cost_function_ref", cost_function_ref)
                print("corrected_cost_function_ref", corrected_cost_function_ref)

                sensitivity = ((np.array(corrected_cost_function_curr) - np.array(corrected_cost_function_ref)) / (delta))

                print("Element ", index+1 ," , Adjoint Sensitivity : " , gradients[index]*E_orig, " vs, FD : ", np.sum(sensitivity)/num_of_timesteps_si)
                self.assertAlmostEqual(gradients[index]*E_orig, np.sum(sensitivity)/num_of_timesteps_si, 9)


                element.Properties[Kratos.YOUNG_MODULUS] = E_orig

if __name__ == "__main__":
    kratos_unittest.main()