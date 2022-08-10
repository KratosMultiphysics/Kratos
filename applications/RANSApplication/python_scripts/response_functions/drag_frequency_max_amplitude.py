"""
This module contains an interface to the available response functions
"""
import math, shutil
import time as timer
import json
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import GetDragValues
from KratosMultiphysics.RANSApplication.response_functions.utilities import RecursiveCopy
from KratosMultiphysics.RANSApplication.response_functions.utilities import CalculateShapeSensitivity
from KratosMultiphysics.RANSApplication.response_functions.utilities import UpdateFilesWithPlaceHolders
from KratosMultiphysics.RANSApplication.response_functions.utilities import UpdateStringWithPlaceHolders

def CalculateMaxFrequencyData(
        parameters: Kratos.Parameters,
        min_frequency: float,
        max_frequency: float,
        windowing_length: float,
        drag_direction,
        drag_model_part_name: str):
    # read time step reaction values
    time_steps, reactions = GetDragValues(parameters, drag_model_part_name)
    delta_time = time_steps[1] - time_steps[0]

    # compute the drag
    number_of_steps = len(time_steps)
    drag_values = [0.0] * number_of_steps
    for i, reaction in enumerate(reactions):
        drag_values[i] = reaction[0] * drag_direction[0] + reaction[1] * drag_direction[1] + reaction[2] * drag_direction[2]

    fluid_fft_utilities = KratosCFD.FluidFFTUtilities(time_steps[-1], windowing_length, delta_time)

    frequency_list, frequency_real_components, frequency_imag_components, frequency_amplitudes_squares = fluid_fft_utilities.CalculateFFTFrequencyDistribution(drag_values)
    frequency_amplitudes = [math.sqrt(_v) for _v in frequency_amplitudes_squares]

    # now get the maximum amplitude frequency from the chosen range
    max_frequency_bin_index = -1
    for index, (frequency, frequency_amplitude)  in enumerate(zip(frequency_list, frequency_amplitudes)):
        if (frequency >= min_frequency and frequency <= max_frequency):
            if (max_frequency_bin_index == -1):
                max_frequency_bin_index = index

            if (frequency_amplitude > frequency_amplitudes[max_frequency_bin_index]):
                max_frequency_bin_index = index

    if (max_frequency_bin_index == -1):
        raise RuntimeError("No frequencies were found in the given range. Please try reducing time step to have more resolution in frequency domain.")

    results_summary_dict = {
        "Frequency information": {
            "Frequency resolution [Hz]": fluid_fft_utilities.GetFrequencyResolution(),
            "Max recorded frequency [Hz]": fluid_fft_utilities.GetMaximumFrequency(),
            "Considered frequency range": {
                "Min frequency [Hz]" : min_frequency,
                "Max frequency [Hz]" : max_frequency
            },
        },
        "Total length information": {
            "Delta time [s]": delta_time,
            "Total length [s]": time_steps[-1],
            "Total steps [#]": fluid_fft_utilities.GetTotalNumberOfSteps()
        },
        "Windowing length information": {
            "Windowing length [s]": windowing_length,
            "Windowing steps [#]": fluid_fft_utilities.GetNumberOfWindowingSteps()
        },
        "Max amplitude data": {
            "Frequency [Hz]": frequency_list[max_frequency_bin_index],
            "Amplitude [N]": frequency_amplitudes[max_frequency_bin_index],
            "Real component [N]": frequency_real_components[max_frequency_bin_index],
            "Imag component [N]": frequency_imag_components[max_frequency_bin_index],
            "Bin index [#]": max_frequency_bin_index
        }
    }

    return fluid_fft_utilities, frequency_list, frequency_real_components, frequency_imag_components, frequency_amplitudes, results_summary_dict

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
                "adjoint_project_parameters_file": "PLEASE_SPECIFY_ADJOINT_PROJECT_PARAMETERS_FILE",
                "place_holder_settings": {
                    "frequency_bin_index_place_holder": "<FREQUENCY_BIN_INDEX>",
                    "frequency_component_place_holder": "<FREQUENCY_COMPONENT_TYPE>",
                    "list_of_files_to_update"         : [
                        {
                            "original_file_name": "",
                            "updated_file_name" : ""
                        }
                    ]
                }
            },
            "clean_primal_solution": false,
            "primal_solution_folder_name": "PLEASE_SPECIFY_PRIMAL_SOLUTION_FOLDER_NAME"
        }
        """)

        self.response_settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.max_frequency_bin_index = -1

        self.problem_setup_folder = Path(self.response_settings["problem_setup_folder"].GetString()).absolute()
        if (not self.problem_setup_folder.is_dir()):
            raise Exception("The provided \"problem_setup_folder\" is not a directory [ \"problem_setup_folder\" = {:s} ].".format(str(self.problem_setup_folder)))

        self.problem_setup_file_settings = self.response_settings["problem_setup_files"]

        # checks for adjoint response settings
        adjoint_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["adjoint_project_parameters_file"].GetString()
        with open(adjoint_project_parameters_file, "r") as file_input:
            input_adjoint_project_parameters_lines = file_input.read()

        input_adjoint_project_parameters_lines = input_adjoint_project_parameters_lines.replace(self.problem_setup_file_settings["place_holder_settings"]["frequency_bin_index_place_holder"].GetString(), str(-1))
        dummy_adjoint_parameters = Kratos.Parameters(input_adjoint_project_parameters_lines)

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
        startTime = timer.time()

        # copy data to the working directory
        RecursiveCopy(str(self.problem_setup_folder), ".")

        write_mesh = True
        if Path(self.mdpa_name + ".mdpa").is_file():
            model = Kratos.Model()
            dummy_model_part = model.CreateModelPart("dummy_model_part")
            Kratos.ModelPartIO(self.mdpa_name, Kratos.IO.READ | Kratos.IO.MESH_ONLY | Kratos.IO.IGNORE_VARIABLES_ERROR).ReadModelPart(dummy_model_part)

            write_mesh = False
            write_mesh = write_mesh or dummy_model_part.NumberOfNodes() != self.updated_model_part.NumberOfNodes()
            write_mesh = write_mesh or dummy_model_part.NumberOfConditions() != self.updated_model_part.NumberOfConditions()
            write_mesh = write_mesh or dummy_model_part.NumberOfElements() != self.updated_model_part.NumberOfElements()

            Kratos.Logger.PrintInfo(self._GetLabel(), "Found existing {:s}.".format(self.mdpa_name + ".mdpa"))
            if write_mesh:
                Kratos.Logger.PrintInfo(self._GetLabel(), "Existing {:s} does not match with the optimization model part. This will be overwritten by optimization model part.".format(self.mdpa_name + ".mdpa"))

        if write_mesh:
            # write the new shape mdpa
            Kratos.ModelPartIO(self.mdpa_name, Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.updated_model_part)
            Kratos.Logger.PrintInfo(self._GetLabel(), "Writing optimization model part to {:s}.".format(self.mdpa_name + ".mdpa"))
        else:
            for node in self.updated_model_part.Nodes:
                node.X = dummy_model_part.GetNode(node.Id).X
                node.Y = dummy_model_part.GetNode(node.Id).Y
                node.Z = dummy_model_part.GetNode(node.Id).Z
            Kratos.Logger.PrintInfo(self._GetLabel(), "Updating optimization model part from {:s}.".format(self.mdpa_name + ".mdpa"))

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed to copy data = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()

        # open primal settings
        primal_settings_file_name = self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        primal_parameters = SolvePrimalProblem(primal_settings_file_name, "primal_evaluation.log")

        self.fluid_fft_utilities, frequency_list, frequency_real_components, frequency_imag_components, frequency_amplitudes, summary = CalculateMaxFrequencyData(primal_parameters, self.min_frequency, self.max_frequency, self.windowing_length, self.drag_direction, self.drag_model_part_name)

        max_amplitude_data = summary["Max amplitude data"]
        self.max_frequency_real_component = max_amplitude_data["Real component [N]"]
        self.max_frequency_imag_component = max_amplitude_data["Imag component [N]"]
        self.max_frequency_amplitude = max_amplitude_data["Amplitude [N]"]
        self.max_frequency_bin_index = max_amplitude_data["Bin index [#]"]

        header = ""
        header += "Primal evaluation summary:\n"
        header += json.dumps(summary, indent=4)

        # now write frequency data to an output file
        with open("drag_frequency_{:s}.csv".format(self.drag_model_part_name), "w") as file_output:
            file_output.write("# Drag frequency output for {:s}".format(self.drag_model_part_name))
            file_output.write("# " + header.replace("\n", "\n# "))
            file_output.write("\n")
            file_output.write("# Index, Frequency [Hz], DFT Real Component [N], DFT Imaginary Component [N], DFT Amplitude [N]\n")
            for i, (frequency, frequency_real_component, frequency_imag_component, frequency_amplitude) in enumerate(zip(frequency_list, frequency_real_components, frequency_imag_components, frequency_amplitudes)):
                file_output.write("{:d}, {:f}, {:f}, {:f}\n".format(i, frequency, frequency_real_component, frequency_imag_component, frequency_amplitude))

        # now print information about found max frequency
        Kratos.Logger.PrintInfo(self._GetLabel(), "\n\n" + header + "\n\n")

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        if (self.max_frequency_bin_index == -1):
            # this is required in the case of iteration restart
            self.CalculateValue()

        # solve adjoint frequency bin problems
        start_time = timer.time()

        real_sensitivities = self._CalculateSensitivities("real", self.problem_setup_file_settings)
        imag_sensitivities = self._CalculateSensitivities("imag", self.problem_setup_file_settings)

        if sorted(list(real_sensitivities.keys())) != sorted(list(imag_sensitivities.keys())):
            raise Exception("Mismatching sensitivities found for real and imaginary component sensitivity calculation.")

        # reset gradients
        self.gradients = {}
        for node_id, real_sensitivity in real_sensitivities.items():
            imag_sensitivity = imag_sensitivities[node_id]

            frequency_amplitude_square_sensitivity = Kratos.Array3(0.0)
            for i in range(3):
                frequency_amplitude_square_sensitivity[i] = self.fluid_fft_utilities.CalculateFFTAmplitudeSquareDerivative(
                        self.max_frequency_real_component,
                        real_sensitivity[i],
                        self.max_frequency_imag_component,
                        imag_sensitivity[i])

            self.gradients[node_id] = frequency_amplitude_square_sensitivity

        if self.response_settings["clean_primal_solution"].GetBool() and Path(self.response_settings["primal_solution_folder_name"].GetString()).is_dir():
            shutil.rmtree(self.response_settings["primal_solution_folder_name"].GetString())

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the total adjoint analysis = ", round(timer.time() - start_time,2),"s")

    def GetValue(self):
        if self.max_frequency_bin_index == -1:
            raise Exception("Please execute CalculateValue method first.")

        return self.max_frequency_amplitude

    def GetNodalGradient(self, variable):
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        return self.gradients

    def IsEvaluatedInFolder(self):
        return True

    def _CalculateSensitivities(self, component_type, settings):
        start_adjoint_simulation = timer.time()

        # update the list of files with place holders
        place_holder_settings = settings["place_holder_settings"]
        place_holder_dict = {
            place_holder_settings["frequency_bin_index_place_holder"].GetString(): str(self.max_frequency_bin_index),
            place_holder_settings["frequency_component_place_holder"].GetString(): component_type
        }

        for place_holder_file_settings in place_holder_settings["list_of_files_to_update"]:
            UpdateFilesWithPlaceHolders(place_holder_file_settings, place_holder_dict)

        adjoint_project_parameters_file = settings["adjoint_project_parameters_file"].GetString()
        adjoint_base_name = adjoint_project_parameters_file[:-5]
        with open(adjoint_project_parameters_file, "r") as file_input:
            input_adjoint_project_parameters_lines = file_input.read()

        input_adjoint_project_parameters_lines = UpdateStringWithPlaceHolders(input_adjoint_project_parameters_lines, place_holder_dict)
        adjoint_parameters = Kratos.Parameters(input_adjoint_project_parameters_lines)

        adjoint_filled_project_parameters_file_name = adjoint_base_name + "_final_bin_index_{:d}_{:s}.json".format(self.max_frequency_bin_index, component_type)
        adjoint_log_file_name = adjoint_base_name + "_bin_index_{:d}_{:s}.log".format(self.max_frequency_bin_index, component_type)

        DragFrequencyMaxAmplitude._WriteKratosParameters(adjoint_filled_project_parameters_file_name, adjoint_parameters)

        model = Kratos.Model()
        sensitivities = CalculateShapeSensitivity(model, self.drag_model_part_name, adjoint_filled_project_parameters_file_name, adjoint_log_file_name)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint for frequency bin {:d} {:s} component = {:f} s using {:s}".format(self.max_frequency_bin_index, component_type, round(timer.time() - start_adjoint_simulation,2), adjoint_project_parameters_file))

        return sensitivities

    @staticmethod
    def _GetLabel():
        return "DragFrequencyMaxAmplitude"

    @staticmethod
    def _WriteKratosParameters(file_name, parameters):
        with open(file_name, "w") as file_output:
            file_output.write(parameters.PrettyPrintJsonString())

