"""
This module contains an interface to the available response functions
"""
import math
import time as timer
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import SolveAdjointProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import GetDragValues
from KratosMultiphysics.RANSApplication.response_functions.utilities import RecursiveCopy

class DragFrequency(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"       : "drag_frequency",
            "frequency_bins_list" : [1, 2, 3, 4, 5],
            "problem_setup_folder": "PLEASE_SPECIFY_PROBLEM_SETUP_FOLDER",
            "problem_setup_files" : {
                "primal_project_parameters_file" : "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE",
                "adjoint_project_parameters_file": "PLEASE_SPECIFY_ADJOINT_PROJECT_PARAMETERS_FILE"
            }
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

        for frequency_bin_index in self.response_settings["frequency_bins_list"].GetVector():
            current_frequency_bin_index = int(frequency_bin_index)
            self.frequency_bin_indices.append(current_frequency_bin_index)
            current_parameters = Kratos.Parameters(input_adjoint_project_parameters_lines.replace("<FREQUENCY_BIN_INDEX>", str(current_frequency_bin_index)))

            adjoint_response_function_settings = current_parameters["solver_settings"]["response_function_settings"]["custom_settings"]
            self.drag_model_part_name = adjoint_response_function_settings["structure_model_part_name"].GetString()
            self.drag_direction = adjoint_response_function_settings["drag_direction"].GetVector()
            self.total_number_of_time_steps = adjoint_response_function_settings["total_number_of_time_steps"].GetInt()
            self.number_of_windowing_steps = adjoint_response_function_settings["window_size"].GetInt()

            current_parameters["solver_settings"]["response_function_settings"]["custom_settings"]["frequency_bin_index"].SetInt(current_frequency_bin_index)

            # add the real component adjoint parameters
            current_parameters["solver_settings"]["response_function_settings"]["custom_settings"]["is_real_component"].SetBool(True)
            self.adjoint_project_parameters_list.append(current_parameters)

            # add the imaginary component adjoint parameters
            current_parameters = current_parameters.Clone()
            current_parameters["solver_settings"]["response_function_settings"]["custom_settings"]["is_real_component"].SetBool(False)
            self.adjoint_project_parameters_list.append(current_parameters)

        self.frequency_real_components = []
        self.frequency_imag_components = []
        self.frequency_amplitudes = []

        # check for model part name settings
        primal_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        with open(primal_project_parameters_file, "r") as file_input:
            primal_settings = Kratos.Parameters(file_input.read())

        self.mdpa_name = primal_settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()

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

        # read time step reaction values
        _, reactions = GetDragValues(primal_parameters, self.drag_model_part_name)
        number_of_steps = len(reactions)

        if (number_of_steps != self.total_number_of_time_steps):
            Kratos.Logger.PrintWarning(self.__class__.__name__, "Total number of steps given in json settings mismatches with number of steps found from primal evaluation. [ Total number of steps in json = {:d}, total number of steps found from primal = {:d} ].".format(
                number_of_steps, self.total_number_of_time_steps))

        number_of_frequencies = len(self.frequency_bin_indices)

        self.frequency_real_components = [0.0] * number_of_frequencies
        self.frequency_imag_components = [0.0] * number_of_frequencies

        # compute the windowed frequency distribution real and imaginary components
        for index, reaction in enumerate(reactions[number_of_steps - self.number_of_windowing_steps:]):
            current_step = index + (number_of_steps - self.number_of_windowing_steps)
            windowing_value = 0.5 * (1.0 - math.cos(2.0 * math.pi * index / self.number_of_windowing_steps))
            drag_value = reaction[0] * self.drag_direction[0] + reaction[1] * self.drag_direction[1] + reaction[2] * self.drag_direction[2]

            for index, k in enumerate(self.frequency_bin_indices):
                self.frequency_real_components[index] += windowing_value * math.cos(2.0 * math.pi * current_step * k / number_of_steps) * drag_value
                self.frequency_imag_components[index] += windowing_value * math.sin(2.0 * math.pi * current_step * k / number_of_steps) * drag_value

        # compute windowed frequency distribution amplitudes
        for i in range(number_of_frequencies):
            self.frequency_amplitudes[i] = math.sqrt(self.frequency_real_components[i] ** 2 + self.frequency_imag_components[i] ** 2) * 2 / number_of_steps

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint frequency bin problems
        start_time = timer.time()

        adjoint_base_file_name = self.problem_setup_file_settings["adjoint_project_parameters_file"].GetString()
        adjoint_base_file_name = adjoint_base_file_name[:-5]
        self.gradients = {}
        for index, adjoint_parameters in enumerate(self.adjoint_project_parameters_list):
            start_adjoint_simulation = timer.time()

            current_frequency_bin_index = index // 2
            current_frequency_bin = self.frequency_bin_indices[current_frequency_bin_index]
            if (index % 2 == 0):
                component_type = "real"
                # compute contributions from real part of DFT
                coeff = self.frequency_real_components[current_frequency_bin_index] * 4.0 / ((self.total_number_of_time_steps ** 2) * (self.frequency_amplitudes[current_frequency_bin_index]))
            else:
                component_type = "imaginary"
                # compute contributions from imaginary part of DFT
                coeff = self.frequency_imag_components[current_frequency_bin_index] * 4.0 / ((self.total_number_of_time_steps ** 2) * (self.frequency_amplitudes[current_frequency_bin_index]))

            # write the modified adjoint parameters
            adjoint_parameters_file_name = "{:s}_FBin_{:d}_{:s}.json".format(adjoint_base_file_name, current_frequency_bin, component_type)
            adjoint_log_file_name = "{:s}_FBin_{:d}_{:s}.log".format(adjoint_base_file_name, current_frequency_bin,component_type)
            with open(adjoint_parameters_file_name, "w") as file_output:
                file_output.write(adjoint_parameters.PrettyPrintJsonString())

            model = Kratos.Model()
            _ = SolveAdjointProblem(model, adjoint_parameters_file_name, adjoint_log_file_name)
            adjoint_model_part = model[self.self.model_part_name]

            for node in adjoint_model_part.Nodes:
                self.gradients[node.Id] += node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY) * coeff

            Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint for frequency bin {:d} {:s} component = {:f} s".format(current_frequency_bin, component_type, round(timer.time() - start_adjoint_simulation,2)))

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

    @staticmethod
    def _GetLabel():
        return "AdjointDragFrequency"

