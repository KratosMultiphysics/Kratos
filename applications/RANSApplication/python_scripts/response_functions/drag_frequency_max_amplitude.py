"""
This module contains an interface to the available response functions
"""
import math, shutil
import time as timer
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import SolveAdjointProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import GetDragValues
from KratosMultiphysics.RANSApplication.response_functions.utilities import RecursiveCopy

from KratosMultiphysics.FluidDynamicsApplication.adjoint_stabilization_utilities import CalculateDragFrequencyDistribution

class DragFrequencyMaxAmplitude(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"              : "drag_frequency_max_amplitude",
            "frequency_range"            : [1e-12, 10.0],
            "problem_setup_folder"       : "PLEASE_SPECIFY_PROBLEM_SETUP_FOLDER",
            "problem_setup_files"        : {
                "primal_project_parameters_file" : "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE",
                "adjoint_project_parameters_file": "PLEASE_SPECIFY_ADJOINT_PROJECT_PARAMETERS_FILE"
            },
            "clean_primal_solution": false,
            "primal_solution_folder_name": "PLEASE_SPECIFY_PRIMAL_SOLUTION_FOLDER_NAME"
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

        input_adjoint_project_parameters_lines = input_adjoint_project_parameters_lines.replace("<FREQUENCY_BIN_INDEX>", str(-1))
        dummy_adjoint_parameters = Kratos.Parameters(input_adjoint_project_parameters_lines)

        self.frequency_real_components = []
        self.frequency_imag_components = []
        self.frequency_amplitudes = []
        self.max_frequency_bin_index = -1

        # check for model part name settings
        primal_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        with open(primal_project_parameters_file, "r") as file_input:
            primal_settings = Kratos.Parameters(file_input.read())

        self.mdpa_name = primal_settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()
        self.main_model_part_name = primal_settings["solver_settings"]["model_part_name"].GetString()
        # assumes always simulations start with start time 0.0
        self.total_length = primal_settings["problem_data"]["end_time"].GetDouble()

        self.windowing_length = dummy_adjoint_parameters["solver_settings"]["response_function_settings"]["custom_settings"]["window_time_length"].GetDouble()
        if (self.total_length < self.windowing_length):
            raise RuntimeError("Total duration of the simulation should be greater than or equal to windowing length. [ Total simulation duration = {:f}s, windowing length = {:f}s ]".format(self.total_length, self.windowing_length))

        self.drag_model_part_name = self.main_model_part_name + "." + dummy_adjoint_parameters["solver_settings"]["response_function_settings"]["custom_settings"]["structure_model_part_name"].GetString()
        self.drag_direction = dummy_adjoint_parameters["solver_settings"]["response_function_settings"]["custom_settings"]["drag_direction"].GetVector()

        frequency_range = self.response_settings["frequency_range"].GetVector()
        if (frequency_range.Size() != 2):
            raise RuntimeError("\"frequency_range\" should have two components corresponding to minimum and maximum of the frequencies to consider in Hz.")
        self.min_frequency = min(frequency_range[0], frequency_range[1])
        self.max_frequency = max(frequency_range[0], frequency_range[1])

    def UpdateDesign(self, updated_model_part, variable):
        self.updated_model_part = updated_model_part

    def InitializeSolutionStep(self):
        self.frequency_real_components = []
        self.frequency_imag_components = []
        self.frequency_amplitudes = []
        self.max_frequency_bin_index = -1

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
        primal_parameters = SolvePrimalProblem(primal_settings_file_name, "primal_evaluation.log")

        # read time step reaction values
        time_steps, reactions = GetDragValues(primal_parameters, self.drag_model_part_name)
        frequency_list, self.frequency_real_components, self.frequency_imag_components, self.frequency_amplitudes = CalculateDragFrequencyDistribution(time_steps, reactions, self.drag_direction, self.windowing_length)
        delta_time = time_steps[1] - time_steps[0]

        number_of_steps = len(reactions)
        number_of_frequencies = number_of_steps // 2
        windowing_steps = int(self.windowing_length / delta_time)
        frequency_resolution = 1.0 / (delta_time * number_of_steps)

        # now get the maximum amplitude frequency from the chosen range
        self.max_frequency_bin_index = -1
        for i in range(number_of_frequencies):
            if (frequency_list[i] >= self.min_frequency and frequency_list[i] <= self.max_frequency):
                if (self.max_frequency_bin_index == -1):
                    self.max_frequency_bin_index = i

                if (self.frequency_amplitudes[i] > self.frequency_amplitudes[self.max_frequency_bin_index]):
                    self.max_frequency_bin_index = i

        if (self.max_frequency_bin_index == -1):
            raise RuntimeError("No frequencies were found in the given range. Please try reducing time step to have more resolution in frequency domain.")

        max_frequency_bin_value = self.max_frequency_bin_index * frequency_resolution

        header = ""
        header += "Primal evaluation summary:\n"
        header += "   delta_time                        : {:f} s\n".format(delta_time)
        header += "   Total time length                 : {:f} s\n".format(time_steps[-1])
        header += "   Total time steps                  : {:d}\n".format(int(time_steps[-1] / delta_time))
        header += "   Reaction steps                    : {:d}\n".format(number_of_steps)
        header += "   Windowing length                  : {:f} s\n".format(self.windowing_length)
        header += "   Windowing steps                   : {:d}\n".format(windowing_steps)
        header += "   Frequency resolution              : {:f} Hz\n".format(frequency_resolution)
        header += "   Max frequency                     : {:f} Hz\n".format(frequency_resolution * number_of_frequencies)
        header += "   Max drag amplitude                : {:f} N\n".format(self.frequency_amplitudes[self.max_frequency_bin_index])
        header += "   Max drag amplitude frequency      : {:f} Hz\n".format(max_frequency_bin_value)
        header += "   Max drag amplitude frequency index: {:d}\n".format(self.max_frequency_bin_index)

        # now write frequency data to an output file
        with open("drag_frequency_{:s}.csv".format(self.drag_model_part_name), "w") as file_output:
            file_output.write("# Drag frequency output for {:s}".format(self.drag_model_part_name))
            file_output.write("# " + header.replace("\n", "\n# "))
            file_output.write("\n")
            file_output.write("# Index, Frequency [Hz], DFT Real Component [N], DFT Imaginary Component [N], DFT Amplitude\n")
            for i in range(number_of_frequencies):
                file_output.write("{:d}, {:f}, {:f}, {:f}\n".format(i, frequency_list[i], self.frequency_real_components[i], self.frequency_imag_components[i], self.frequency_amplitudes[i]))

        # now print information about found max frequency
        Kratos.Logger.PrintInfo(self._GetLabel(), "\n\n" + header + "\n\n")

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint frequency bin problems
        start_time = timer.time()

        # reset gradients
        self.gradients = {}

        # run adjoints for solution components of DFT
        model = self._RunAdjointProblem("real", self.problem_setup_file_settings, False)

        adjoint_model_part = model[self.main_model_part_name]
        for node in adjoint_model_part.Nodes:
            self.gradients[node.Id] = Kratos.Array3(0.0)

        self._CalculateSensitivities(model, self.frequency_real_components)
        model = self._RunAdjointProblem("imag", self.problem_setup_file_settings, False)
        self._CalculateSensitivities(model, self.frequency_imag_components)

        if self.response_settings["clean_primal_solution"].GetBool():
            shutil.rmtree(self.response_settings["primal_solution_folder_name"].GetString())

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the total adjoint analysis = ", round(timer.time() - start_time,2),"s")

    def GetValue(self):
        if len(self.frequency_amplitudes) == 0:
            raise Exception("Please execute CalculateValue method first.")

        return self.frequency_amplitudes[self.max_frequency_bin_index]

    def GetNodalGradient(self, variable):
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        return self.gradients

    def IsEvaluatedInFolder(self):
        return True

    def _RunAdjointProblem(self, component_type, settings, skip_if_already_run):
        start_adjoint_simulation = timer.time()

        adjoint_project_parameters_file = settings["adjoint_project_parameters_file"].GetString()
        adjoint_base_name = adjoint_project_parameters_file[:-5]
        with open(adjoint_project_parameters_file, "r") as file_input:
            input_adjoint_project_parameters_lines = file_input.read()

        # set the frequency bin index in text input
        # this is done so that if user uses <FREQUENCY_BIN_INDEX> or <FREQUENCY_COMPONENT_TYPE>
        # flag in some other settings, they will also be replaced properly
        input_adjoint_project_parameters_lines = input_adjoint_project_parameters_lines.replace("<FREQUENCY_BIN_INDEX>", str(self.max_frequency_bin_index))
        input_adjoint_project_parameters_lines = input_adjoint_project_parameters_lines.replace("<FREQUENCY_COMPONENT_TYPE>", component_type)
        adjoint_parameters = Kratos.Parameters(input_adjoint_project_parameters_lines)

        # now set the parameters just to be sure
        adjoint_response_settings = adjoint_parameters["solver_settings"]["response_function_settings"]["custom_settings"]
        adjoint_response_settings["is_real_component"].SetBool(component_type == "real")
        adjoint_response_settings["frequency_bin_index"].SetInt(self.max_frequency_bin_index)

        adjoint_filled_project_parameters_file_name = adjoint_base_name + "_final_{:s}.json".format(component_type)
        adjoint_log_file_name = adjoint_base_name + "_{:s}.log".format(component_type)

        DragFrequencyMaxAmplitude._WriteKratosParameters(adjoint_filled_project_parameters_file_name, adjoint_parameters)

        model = Kratos.Model()
        _ = SolveAdjointProblem(model, adjoint_filled_project_parameters_file_name, adjoint_log_file_name, skip_if_already_run)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint for frequency bin {:d} {:s} component = {:f} s using {:s}".format(self.max_frequency_bin_index, component_type, round(timer.time() - start_adjoint_simulation,2), adjoint_project_parameters_file))

        return model

    def _CalculateSensitivities(self, model, components_list):
        adjoint_model_part = model[self.main_model_part_name]
        delta_time = -1.0 * adjoint_model_part.ProcessInfo[Kratos.DELTA_TIME]
        windowing_steps = int(self.windowing_length / delta_time)

        coeff = components_list[self.max_frequency_bin_index] * 8.0 / (windowing_steps ** 2)

        for node in adjoint_model_part.Nodes:
            self.gradients[node.Id] += node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY) * coeff

    @staticmethod
    def _GetLabel():
        return "DragFrequencyMaxAmplitude"

    @staticmethod
    def _WriteKratosParameters(file_name, parameters):
        with open(file_name, "w") as file_output:
            file_output.write(parameters.PrettyPrintJsonString())

