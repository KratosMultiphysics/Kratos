import KratosMultiphysics as KM
from KratosMultiphysics.time_based_ascii_file_writer_utility import (
    TimeBasedAsciiFileWriterUtility,
)
import KratosMultiphysics.HDF5Application as KratosHDF5


import numpy as np
import operator as op


def Factory(settings: KM.Parameters, model: KM.Model) -> KM.Process:
    if not isinstance(settings, KM.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if not isinstance(model, KM.Model):
        raise Exception("expected input shall be a Model object.")
    return DataGeneratorProcess(model, settings["Parameters"])


## All the processes python should be derived from "Process"
class DataGeneratorProcess(KM.Process):
    """This class generates a dataset for Neural Network training and testing

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KM.Process.__init__(self)  # calling the baseclass constructor

        default_settings = KM.Parameters(
            """{
            "help"                     : "This process generates the data series for the neural network training and testing",
            "interval"                 : [0, 1e30],
            "write_output_file"        : true,
            "input_variables"          : [""],
            "input_sources"            : [""],
            "input_order"              : "variables_first",
            "output_variables"         : [""],
            "output_sources"           : [""],
            "output_format"            : "ascii",
            "output_order"            : "variables_first",
            "perturbate_variables"     : [""],
            "random_distribution"      : [""],
            "random_parameters"        : [[]],
            "print_format"             : ".8f",
            "training_file_settings": {
                    "input_file_name"  : "training_in",
                    "output_file_name" : "training_out",
                    "file_format"      : "dat",
                    "output_path"      : "data/"
            }
        }"""
        )

        # getting the ModelPart from the Model
        self.output_model_parts = []
        if not settings.Has("model_part_name"):
            raise Exception("Output model part not present")
        if settings["model_part_name"].IsArray():
            output_is_vector = True
            self.output_model_parts_names = settings["model_part_name"].GetStringArray()
        else:
            output_is_vector = False
            self.output_model_parts_names = [settings["model_part_name"].GetString()]
        if not self.output_model_parts_names:
            raise Exception('No "model_part_name" was specified!')
        for name in self.output_model_parts_names:
            self.output_model_parts.append(model[name])
        settings.RemoveValue("model_part_name")

        # getting the input ModelPart from the Model
        self.input_model_parts = []
        if not settings.Has("input_model_part"):
            raise Exception("Input model part not present")
        if settings["input_model_part"].IsArray():
            input_is_vector = True
            self.input_model_parts_names = settings["input_model_part"].GetStringArray()
        else:
            input_is_vector = False
            self.input_model_parts_names = [settings["input_model_part"].GetString()]
        if not self.input_model_parts_names:
            raise Exception('No "input_model_part" was specified!')
        for name in self.input_model_parts_names:
            self.input_model_parts.append(model[name])
        settings.RemoveValue("input_model_part")

        settings.ValidateAndAssignDefaults(default_settings)

        # Detect 'End' as a tag and replace it by a large number
        if settings.Has("interval"):
            if settings["interval"][1].IsString():
                if settings["interval"][1].GetString() == "End":
        self.random_distribution = settings["random_distribution"].GetStringArray()
            KM.KratosGlobals.GetVariable(var) for var in perturbate_variable_names
        ]
        if not self.perturbate_variables:
            self.perturbate_variables.append(None)
        input_dist_names = settings["random_distribution"]
        self.random_distribution = [
            input_dist_names[i].GetString() for i in range(input_dist_names.size())
        ]

        # getting the list of distribution parameters for each variable
        random_parameters_matrix = settings["random_parameters"].GetMatrix()
        self.random_parameters = []
        for distribution in range(random_parameters_matrix.Size1()):
            distribution_parameters = []
            for parameter in range(random_parameters_matrix.Size2()):
                distribution_parameters.append(
                    random_parameters_matrix[distribution, parameter]
                )
            self.random_parameters.append(distribution_parameters)
        self.interval = settings["interval"].GetVector()
        self.format = settings["print_format"].GetString()
        self.output_format = settings["output_format"].GetString()

        # retrieving the input variables
        input_var_names = settings["input_variables"]
        variable_names = [
            input_var_names[i].GetString() for i in range(input_var_names.size())
        ]
        input_sources_names = settings["input_sources"]
        self.input_sources = [
            input_sources_names[i].GetString()
            for i in range(input_sources_names.size())
        ]
        self.input_variables = [
            KM.KratosGlobals.GetVariable(var) for var in variable_names
        ]
        if len(self.input_variables) == 0:
            raise Exception("No variables specified for input!")
        if not (len(self.input_variables) == len(self.input_sources)):
            raise Exception("The number of input variables and sources are different.")
        if not (len(self.input_variables) == len(self.input_model_parts)):
            if input_is_vector:
                raise Exception(
                    "The input_model_parts are given as vector -- The number of input variables and model parts must be the same."
                )
            else:
                while not (len(self.input_variables) == len(self.input_model_parts)):
                    self.input_model_parts.append(self.input_model_parts[-1])
        self.dict_input = list(
            zip(self.input_model_parts, self.input_variables, self.input_sources)
        )
        # getting input order
        self.input_order = settings["input_order"].GetString()

        # retrieving the output variables
        output_var_names = settings["output_variables"]
        variable_names = [
            output_var_names[i].GetString() for i in range(output_var_names.size())
        ]
        output_sources_names = settings["output_sources"]
        self.output_sources = [
            output_sources_names[i].GetString()
            for i in range(output_sources_names.size())
        ]
        self.output_variables = [
            KM.KratosGlobals.GetVariable(var) for var in variable_names
        ]
        if len(self.output_variables) == 0:
            raise Exception("No variables specified for output!")
        if not (len(self.output_variables) == len(self.output_sources)):
            raise Exception("The number of output variables and sources are different.")
        if not (len(self.output_variables) == len(self.output_model_parts)):
            if output_is_vector:
                raise Exception(
                    "The output_model_parts are given as vector -- The number of output variables and model parts must be the same."
                )
            else:
                while not (len(self.output_variables) == len(self.output_model_parts)):
                    self.output_model_parts.append(self.output_model_parts[-1])
        self.dict_output = list(
            zip(self.output_model_parts, self.output_variables, self.output_sources)
        )
        # getting output order
        self.output_order = settings["output_order"].GetString()

        # Output file names handling
        if self.write_output_file:
            model_part_name = self.output_model_parts_names[0]
            if model_part_name == "":
                raise Exception('No "model_part_name" was specified!')
            self.model_part = model[model_part_name]

            training_input_file_name = (
                settings["training_file_settings"]["input_file_name"].GetString()
                + "."
                + settings["training_file_settings"]["file_format"].GetString()
            )
            training_output_file_name = (
                settings["training_file_settings"]["output_file_name"].GetString()
                + "."
                + settings["training_file_settings"]["file_format"].GetString()
            )

            training_input_file_handler_params = KM.Parameters()
            training_output_file_handler_params = KM.Parameters()

            if self.output_format == "ascii":
                training_input_file_handler_params.AddEmptyValue("file_name")
                training_input_file_handler_params["file_name"].SetString(
                    training_input_file_name
                )
                training_input_file_handler_params.AddEmptyValue("output_path")
                training_input_file_handler_params["output_path"].SetString(
                    settings["training_file_settings"]["output_path"].GetString()
                )
                training_output_file_handler_params.AddEmptyValue("file_name")
                training_output_file_handler_params["file_name"].SetString(
                    training_output_file_name
                )
                training_output_file_handler_params.AddEmptyValue("output_path")
                training_output_file_handler_params["output_path"].SetString(
                    settings["training_file_settings"]["output_path"].GetString()
                )

                # Training input
                file_header = ""
                self.training_input_file = TimeBasedAsciiFileWriterUtility(
                    self.model_part, training_input_file_handler_params, file_header
                ).file

                # Training output
                self.training_output_file = TimeBasedAsciiFileWriterUtility(
                    self.model_part, training_output_file_handler_params, file_header
                ).file
            elif self.output_format == "h5":
                # Training input
                training_input_file_name = (
                    settings["training_file_settings"]["output_path"].GetString()
                    + training_input_file_name
                )
                training_input_file_handler_params.AddEmptyValue("file_access_mode")
                training_input_file_handler_params["file_access_mode"].SetString(
                    "truncate"
                )
                training_input_file_handler_params.AddEmptyValue("file_name")
                training_input_file_handler_params["file_name"].SetString(
                    training_input_file_name
                )
                self.hdf5_file_training_input = KratosHDF5.HDF5FileSerial(
                    training_input_file_handler_params
                )

                self.hdf5_input_parameters = KM.Parameters()
                self.hdf5_input_parameters.AddEmptyArray("list_of_variables")
                self.hdf5_input_parameters["list_of_variables"].SetStringArray(
                    input_var_names.GetStringArray()
                )
                self.hdf5_input_parameters.AddEmptyValue("prefix")
                # self.hdf5_input_parameters["prefix"].SetString("/<time>/<model_part_name>/InputData")

                # Training output
                training_output_file_name = (
                    settings["training_file_settings"]["output_path"].GetString()
                    + training_output_file_name
                )
                training_output_file_handler_params.AddEmptyValue("file_access_mode")
                training_output_file_handler_params["file_access_mode"].SetString(
                    "truncate"
                )
                training_output_file_handler_params.AddEmptyValue("file_name")
                training_output_file_handler_params["file_name"].SetString(
                    training_output_file_name
                )
                self.hdf5_file_training_output = KratosHDF5.HDF5FileSerial(
                    training_output_file_handler_params
                )

                self.hdf5_output_parameters = KM.Parameters()
                self.hdf5_output_parameters.AddEmptyArray("list_of_variables")
                self.hdf5_output_parameters["list_of_variables"].SetStringArray(
                    output_var_names.GetStringArray()
                )
                self.hdf5_output_parameters.AddEmptyValue("prefix")
                # self.hdf5_output_parameters["prefix"].SetString("/<time>/<model_part_name>/ResultsData")

            else:
                raise Exception(
                    "Format not supported. Supported formats are: ascii, h5."
                )

    def ExecuteInitializeSolutionStep(self):
        current_time = self.input_model_parts[0].ProcessInfo[KM.TIME]
        current_step = self.input_model_parts[0].ProcessInfo[KM.STEP]
        if (current_time >= self.interval[0]) and (current_time < self.interval[1]):
            index = 0
            input_value_list = []
            for model_part, variable, source in self.dict_input:
                model_input_value_list = []
                if (
                    index < len(self.perturbate_variables)
                    and variable == self.perturbate_variables[index]
                ):
                    dist = self.random_distribution[index]
                    rnd_parameters = self.random_parameters[index]
                    factor = getattr(np.random, dist)(*rnd_parameters)
                    # Process related variables (e.g. TIME, STEP)
                    if source == "process":
                        raise Exception("Process variables cannot be perturbated.")
                    # Node properties (e.g. position)
                    elif source == "node":
                        print("Warning: The node properties are being perturbated.")
                        for node in model_part.Nodes:
                            input_value = op.imul(getattr(node, variable.Name()))
                            node.SetValue(variable, input_value)
                            model_input_value_list.append(
                                getattr(node, variable.Name())
                            )
                    # Node step values (e.g. variables like displacement)
                    elif source == "solution_step":
                        for node in model_part.Nodes:
                            input_value = op.imul(
                                node.GetSolutionStepValue(variable), (1.0 + factor)
                            )
                            node.SetSolutionStepValue(variable, input_value)
                            model_input_value_list.append(
                                node.GetSolutionStepValue(variable, 0)
                            )
                    elif source == "previous_solution_step":
                        for node in model_part.Nodes:
                            input_value = op.imul(
                                node.GetSolutionStepValue(variable), (1.0 + factor)
                            )
                            node.SetSolutionStepValue(variable, input_value)
                            model_input_value_list.append(
                                node.GetSolutionStepValue(variable, 1)
                            )
                    # Condition values
                    elif source == "condition":
                        for condition in model_part.GetConditions():
                            input_value = op.imul(
                                condition.GetValue(variable), (1.0 + factor)
                            )
                            condition.SetValue(variable, input_value)
                            model_input_value_list.append(condition.GetValue(variable))
                    index = index + 1
                else:
                    # Process related variables (e.g. TIME, STEP)
                    if source == "process":
                        model_input_value_list.append(model_part.ProcessInfo[variable])
                    # Node properties (e.g. position)
                    elif source == "node":
                        for node in model_part.Nodes:
                            model_input_value_list.append(
                                getattr(node, variable.Name())
                            )
                    # Node step values (e.g. variables like displacement)
                    elif source == "solution_step":
                        for node in model_part.Nodes:
                            model_input_value_list.append(
                                node.GetSolutionStepValue(variable, 0)
                            )
                    elif source == "previous_solution_step":
                        for node in model_part.Nodes:
                            model_input_value_list.append(
                                node.GetSolutionStepValue(variable, 1)
                            )
                    # Condition values
                    elif source == "condition":
                        for condition in model_part.GetConditions():
                            model_input_value_list.append(condition.GetValue(variable))
                input_value_list.extend(model_input_value_list)
            # Reorder if indicated
            if self.input_order == "sources_first":
                try:
                    input_value_list = self._OrderSourcesFirst(
                        input_value_list, self.dict_input
                    )
                except IndexError:
                    pass

            # Writing input file
            if self.write_output_file:
                if self.output_format == "ascii":
                    self.training_input_file.write(
                        " ".join(str(v) for v in input_value_list) + "\n"
                    )
                if self.output_format == "h5":
                    self.hdf5_input_parameters["prefix"].SetString(
                        "/" + str(current_step) + "/InputData"
                    )
                    nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(
                        self.hdf5_input_parameters, self.hdf5_file_training_input
                    )
                    nodal_io.WriteNodalResults(self.input_model_parts[0], 1)

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.output_model_parts[0].ProcessInfo[KM.TIME]
        current_step = self.output_model_parts[0].ProcessInfo[KM.STEP]
        if (current_time >= self.interval[0]) and (current_time < self.interval[1]):
            output_value_list = []
            for model_part, variable, source in self.dict_output:
                model_output_value_list = []
                # Process related variables (e.g. TIME, STEP)
                if source == "process":
                    model_output_value_list.append(model_part.ProcessInfo[variable])
                # Node properties (e.g. position)
                elif source == "node":
                    for node in model_part.Nodes:
                        model_output_value_list.append(getattr(node, variable.Name()))
                # Node step values (e.g. variables like displacement)
                elif source == "solution_step":
                    for node in model_part.Nodes:
                        model_output_value_list.append(
                            node.GetSolutionStepValue(variable, 0)
                        )
                # Condition values
                elif source == "condition":
                    for condition in model_part.GetConditions():
                        model_output_value_list.append(condition.GetValue(variable))
                output_value_list.extend(model_output_value_list)
            # Reorder if indicated
            if self.output_order == "sources_first":
                try:
                    output_value_list = self._OrderSourcesFirst(
                        output_value_list, self.dict_output
                    )
                except IndexError:
                    pass

            # Writing output file
            if self.write_output_file:
                if self.output_format == "ascii":
                    self.training_output_file.write(
                        " ".join(str(v) for v in output_value_list) + "\n"
                    )
                if self.output_format == "h5":
                    self.hdf5_output_parameters["prefix"].SetString(
                        "/" + str(current_step) + "/OutputData"
                    )
                    nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(
                        self.hdf5_output_parameters, self.hdf5_file_training_output
                    )
                    nodal_io.WriteNodalResults(self.output_model_parts[0], 1)

    def IsOutputStep(self):
        return False

    @staticmethod
    def _OrderSourcesFirst(values_list, variables_dictionary):
        """Reorders the values list by sources instead of by variables (e.g. node by node)."""
        ordered_list = []
        k, m = divmod(len(values_list), len(variables_dictionary))
        # Split the lists
        split_list = list(
            values_list[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)]
            for i in range(len(variables_dictionary))
        )
        # Reorder the lists by source entry
        for index in range(len(split_list[0])):
            for variable_list in split_list:
                ordered_list.append(variable_list[index])

        return ordered_list
