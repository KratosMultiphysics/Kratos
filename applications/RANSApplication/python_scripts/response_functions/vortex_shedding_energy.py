"""
This module contains an interface to the available response functions
"""
import math
import time as timer
from pathlib import Path
from csv import reader

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import SolveAdjointProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import RecursiveCopy
from KratosMultiphysics.RANSApplication import RansAuxiliaryUtilities

class VortexSheddingEnergy(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"       : "vortex_shedding",
            "problem_setup_folder": "PLEASE_SPECIFY_PROBLEM_SETUP_FOLDER",
            "problem_setup_files" : {
                "primal_project_parameters_file" : "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE",
                "adjoint_project_parameters_file": "PLEASE_SPECIFY_ADJOINT_PROJECT_PARAMETERS_FILE"
            },
            "vortex_shedding_frequency_bins_list": [1, 2, 3, 4, 5]
        }
        """)

        self.response_settings.ValidateAndAssignDefaults(default_parameters)

        self.problem_setup_folder = Path(self.response_settings["problem_setup_folder"].GetString()).absolute()
        if (not self.problem_setup_folder.is_dir()):
            raise Exception("The provided \"problem_setup_folder\" is not a directory [ \"problem_setup_folder\" = {:s} ].".format(str(self.problem_setup_folder)))

        self.problem_setup_file_settings = self.response_settings["problem_setup_files"]

        # checks for adjoint response settings
        adjoint_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["adjoint_project_parameters_file"].GetString()
        with open(adjoint_project_parameters_file, "r") as file_input:
            input_adjoint_project_parameters_lines = file_input.read()

        self.adjoint_project_parameters_list = []
        self.frequency_bin_indices = []

        for frequency_bin_index in self.response_settings["vortex_shedding_frequency_bins_list"].GetVector():
            current_frequency_bin_index = int(frequency_bin_index)
            self.frequency_bin_indices.append(current_frequency_bin_index)
            current_parameters = Kratos.Parameters(input_adjoint_project_parameters_lines.replace("<VORTEX_SHEDDING_FREQUENCY_BIN_INDEX>", str(current_frequency_bin_index)))

            adjoint_response_function_settings = current_parameters["solver_settings"]["response_function_settings"]["custom_settings"]
            self.model_part_name = adjoint_response_function_settings["model_part_name"].GetString()
            self.point_coordinates = adjoint_response_function_settings["point_coordinates"].GetVector()
            self.velocity_direction = adjoint_response_function_settings["velocity_direction"].GetVector()
            self.total_number_of_time_steps = adjoint_response_function_settings["total_number_of_time_steps"].GetInt()
            self.number_of_windowing_steps = adjoint_response_function_settings["window_size"].GetInt()

            current_parameters["frequency_bin_index"].SetInt(frequency_bin_index)

            # add the real component adjoint parameters
            current_parameters["is_real_component"].SetBool(True)
            self.adjoint_project_parameters_list.append(current_parameters)

            # add the imaginary component adjoint parameters
            current_parameters = current_parameters.Clone()
            current_parameters["is_real_component"].SetBool(False)
            self.adjoint_project_parameters_list.append(current_parameters)

        self.frequency_real_components = []
        self.frequency_imag_components = []
        self.frequency_amplitudes = []

    def UpdateDesign(self, updated_model_part, variable):
        self.updated_model_part = updated_model_part

    def InitializeSolutionStep(self):
        self.frequency_real_components = []
        self.frequency_imag_components = []
        self.frequency_amplitudes = []

        startTime = timer.time()

        # copy data to the working directory
        RecursiveCopy(str(self.problem_setup_folder), ".")

        # write the new shape mdpa
        Kratos.ModelPartIO(self.mdpa_name, Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.updated_model_part)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed to copy data = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()

        # open primal settings
        primal_settings_file_name = self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        primal_parameters = SolvePrimalProblem(primal_settings_file_name)

        # get primal output process and settings
        primal_frequency_bin_indices_list, primal_output_process_parameters = self._GetResponseFunctionOutputProcess(primal_parameters)
        if (primal_output_process_parameters is None):
            raise Exception("No FrequencyBinOutputProcess is found in primal project parameters.")

        # read data from the csv file
        output_file_name = primal_output_process_parameters["Parameters"]["output_file_settings"]["file_name"].GetString()
        real_data, imag_data = self._GetFrequencyBinComponentValues(output_file_name, primal_frequency_bin_indices_list)

        # apply Hann window
        self.real_windowed_data = Kratos.Vector()
        self.imag_windowed_data = Kratos.Vector()
        for real_data_vector, imag_data_vector in zip(real_data, imag_data):
            RansAuxiliaryUtilities.ApplyHannWindow(self.real_windowed_data, real_data_vector, 0, self.number_of_windowing_steps - 1)
            RansAuxiliaryUtilities.ApplyHannWindow(self.imag_windowed_data, imag_data_vector, 0, self.number_of_windowing_steps - 1)

        for real_data_vector, imag_data_vector in zip(self.real_windowed_data, self.imag_windowed_data):
            self.frequency_real_components.append(RansAuxiliaryUtilities.VectorSummation(real_data_vector))
            self.frequency_imag_components.append(RansAuxiliaryUtilities.VectorSummation(imag_data_vector))
            self.frequency_amplitudes.append(math.sqrt(self.frequency_real_components[-1] ** 2 + self.frequency_imag_components[-1] ** 2))

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint lift problem
        start_time = timer.time()

        adjoint_base_file_name = self.problem_setup_file_settings["adjoint_project_parameters_file"].GetString()
        adjoint_base_file_name = adjoint_base_file_name[:-5]
        self.gradients = {}
        for index, adjoint_parameters in enumerate(self.adjoint_project_parameters_list):
            start_adjoint_simulation = timer.time()

            # write the modified adjoint parameters
            adjoint_parameters_file_name = "{:s}_FBin_{:d}.json".format(adjoint_base_file_name, index + 1)
            adjoint_log_file_name = "{:s}_FBin_{:d}.log".format(adjoint_base_file_name, index + 1)
            with open(adjoint_parameters_file_name, "w") as file_output:
                file_output.write(adjoint_parameters.PrettyPrintJsonString())

            model = Kratos.Model()
            _ = SolveAdjointProblem(model, adjoint_parameters_file_name, adjoint_log_file_name)
            adjoint_model_part = model[self.self.model_part_name]

            if index % 2 == 0:
                # compute contributions from real part of DFT
                coeff = self.frequency_real_components[index] * 4 / ((self.total_number_of_time_steps ** 2) * (self.frequency_amplitudes[index]))
                component_type = "real"
            else:
                # compute contributions from imaginary part of DFT
                coeff = self.frequency_imag_components[index] * 4 / ((self.total_number_of_time_steps ** 2) * (self.frequency_amplitudes[index]))
                component_type = "imaginary"

            for node in adjoint_model_part.Nodes:
                self.gradients[node.Id] += node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY) * coeff

            Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint for frequency bin index {:d} {:s} component = {:f} s".format(index + 1, component_type, round(timer.time() - start_adjoint_simulation,2)))

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the total adjoint analysis = ",round(timer.time() - start_time,2),"s")

    def GetValue(self):
        if len(self.frequency_amplitudes) == 0:
            raise Exception("Please execute CalculateValue method first.")

        frequency_amplitude_sum = 0.0
        for v in self.frequency_amplitudes:
            frequency_amplitude_sum += v

        return frequency_amplitude_sum

    def GetNodalGradient(self, variable):
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        return self.gradients

    def IsEvaluatedInFolder(self):
        return True

    def _GetFrequencyBinComponentValues(self, output_csv_file_name, indices_list):
        with open(output_csv_file_name, "r") as file_input:
            data = reader(filter(lambda row: row[0]!='#', file_input), delimiter=",")

        # here we need to get the last time step data to apply the windowing
        real_data_values = []
        imag_data_values = []
        for _ in indices_list:
            real_data_values.append(Kratos.Vector(self.number_of_windowing_steps))
            imag_data_values.append(Kratos.Vector(self.number_of_windowing_steps))

        for row_index, row in enumerate(data):
            current_offsetted_index = row_index - self.total_number_of_time_steps + self.number_of_windowing_steps
            if (current_offsetted_index >= 0):
                if (current_offsetted_index >= self.number_of_windowing_steps):
                    raise Exception("Number of time steps read in the file is larger than the provided value in the settings.")

                for value_index, primal_index in enumerate(indices_list):
                    real_data_values[value_index][current_offsetted_index] = float(row[1 + primal_index * 2])
                    imag_data_values[value_index][current_offsetted_index] = float(row[1 + primal_index * 2 + 1])

        return real_data_values, imag_data_values

    def _GetResponseFunctionOutputProcess(self, kratos_parameters):
        output_process_list = kratos_parameters["output_processes"]
        for output_category in output_process_list:
            for output_process_settings in output_category:
                if (
                    output_process_settings.Has("python_module") and output_process_settings["python_module"].GetString() == "frequency_bin_output_process" and
                    output_process_settings.Has("kratos_module") and output_process_settings["kratos_module"].GetString() == "KratosMultiphysics.RANSApplication" and
                    output_process_settings.Has["Parameters"].Has("model_part_name") and output_process_settings.Has["Parameters"]["model_part_name"].GetString() == self.model_part_name and
                    output_process_settings.Has["Parameters"].Has("point_coordinates") and output_process_settings.Has["Parameters"]["point_coordinates"].GetVector() == self.point_coordinates and
                    output_process_settings.Has["Parameters"].Has("velocity_direction") and output_process_settings.Has["Parameters"]["velocity_direction"].GetVector() == self.velocity_direction and
                    output_process_settings.Has["Parameters"].Has("total_number_of_time_steps") and output_process_settings.Has["Parameters"]["total_number_of_time_steps"].GetInt() == self.total_number_of_time_steps
                ):
                    f_bin_indices = output_process_settings.Has["Parameters"]["frequency_bin_indices"].GetVector()

                    found_all_bin_indices = True
                    primal_bin_index_list = []
                    for v_adjoint in self.frequency_bin_indices:
                        found_adjoint_bin_index = False
                        for index, v_primal in enumerate(f_bin_indices):
                            if int(v_primal) == v_adjoint:
                                primal_bin_index_list.append(index)
                                found_adjoint_bin_index = True
                                break

                        if (not found_adjoint_bin_index):
                            found_all_bin_indices = False
                            break

                    if (found_all_bin_indices):
                        return primal_bin_index_list, output_process_settings
                    else:
                        return None

        return None

    @staticmethod
    def _GetLabel():
        return "AdjointVortexSheddingEnergy"

