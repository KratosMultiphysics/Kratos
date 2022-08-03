# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

if CheckIfApplicationsAvailable("StatisticsApplication"):
    from KratosMultiphysics.StatisticsApplication import SpatialMethods as spatial_methods
else:
    msg = "YPlusOutputProcess requires StatisticsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)

from math import sqrt

def Factory(settings, model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return YPlusOutputProcess(model, settings["Parameters"])


class YPlusOutputProcess(Kratos.OutputProcess):
    """
    A class responsible for the Y+ output, which is an condition value in Kratos.
    """

    def __init__(self, model, params):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""
            {
                "model_part_name"                  : "",
                "interval"                         : [0.0, 1e30],
                "y_plus_output_limit"              : 1000.0,
                "print_to_screen"                  : false,
                "print_format"                     : ".8f",
                "write_output_file"                : true,
                "output_step"                      : 8,
                "output_to_elements"               : false,
                "calculate_normals_every_time_step": false,
                "echo_level"                       : 0,
                "output_file_settings"             : {}
            }
            """)

        self.interval_utility = Kratos.IntervalUtility(params)
        params.ValidateAndAssignDefaults(default_settings)

        # getting the ModelPart from the Model
        self.model_part_name = params["model_part_name"].GetString()
        if self.model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.model_part_name]

        # getting output limit for summarization
        self.y_plus_output_limit = params["y_plus_output_limit"].GetDouble()

        self.format = params["print_format"].GetString()
        self.output_step = params["output_step"].GetInt()
        self.print_to_screen = params["print_to_screen"].GetBool()
        self.write_output_file = params["write_output_file"].GetBool()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):

                output_file_name = params["model_part_name"].GetString() + "_y_plus.dat"

                file_handler_params = Kratos.Parameters(
                    params["output_file_settings"])

                if file_handler_params.Has("file_name"):
                    warn_msg  = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                    warn_msg += '"' + file_handler_params["file_name"].GetString() + '"}\n'
                    warn_msg += 'Using this specififed file name instead of the default "' + output_file_name + '"'
                    Kratos.Logger.PrintWarning("YPlusOutputProcess", warn_msg)
                else:
                    file_handler_params.AddEmptyValue("file_name")
                    file_handler_params["file_name"].SetString(output_file_name)

                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                   file_handler_params, file_header).file

        self.distribution_params = Kratos.Parameters('''{
            "number_of_value_groups" : 1,
            "min_value"              : 0.0,
            "max_value"              : "max"
        }''')
        self.distribution_params["max_value"].SetDouble(self.y_plus_output_limit)


        compute_y_plus_process_settings = Kratos.Parameters("""{
            "model_part_name"                  : "",
            "output_variable_name"             : "Y_PLUS",
            "output_to_elements"               : false,
            "calculate_normals_every_time_step": false,
            "echo_level"                       : 0
        }""")
        compute_y_plus_process_settings["model_part_name"].SetString(params["model_part_name"].GetString())
        compute_y_plus_process_settings["output_to_elements"].SetBool(params["output_to_elements"].GetBool())
        compute_y_plus_process_settings["calculate_normals_every_time_step"].SetBool(params["calculate_normals_every_time_step"].GetBool())
        compute_y_plus_process_settings["echo_level"].SetInt(params["echo_level"].GetInt())
        self.compute_y_plus_process = KratosCFD.ComputeYPlusProcess(model, compute_y_plus_process_settings)

    def Check(self):
        return self.compute_y_plus_process.Check()

    def ExecuteInitializeSolutionStep(self):
        if self.IsOutputStep():
            self.compute_y_plus_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        if self.IsOutputStep():
            self.compute_y_plus_process.ExecuteFinalizeSolutionStep()

    def PrintOutput(self):
        output = self._CalculateWithRespectToThreshold()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            current_time = self.model_part.ProcessInfo[Kratos.TIME]
            output_vals = [format(val, self.format) for val in output]

            # not formatting time in order to not lead to problems with time recognition
            # in the file writer when restarting
            output_vals.insert(0, str(current_time))

            res_labels = ["time: ", "mean: ", "std: ", "max: ", "y_plus > " +
                            "{:.1f}".format(self.y_plus_output_limit) + ": "]

            if (self.print_to_screen):
                result_msg = 'y+ evaluation for model part ' + \
                    self.model_part_name + '\n'
                result_msg += ', '.join([a+b for a,
                                            b in zip(res_labels, output_vals)])
                self._PrintToScreen(result_msg + "%")

            if (self.write_output_file):
                self.output_file.write(' '.join(output_vals) + "\n")

    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

    def _GetFileHeader(self):
        header = '# Y+ for model part ' + self.model_part_name + \
            '| Y+_threshold: ' + str(self.y_plus_output_limit) + '\n'
        header += '# Time Mean Std Max HowMany >' + \
            "{:.1f}".format(self.y_plus_output_limit) + ' [%]\n'
        return header

    def _PrintToScreen(self, result_msg):
        Kratos.Logger.PrintInfo("YPlusOutputProcess", "Y+ VALUE RESULTS:")
        Kratos.Logger.PrintInfo("YPlusOutputProcess", "Current time: " + result_msg)

    def _CalculateWithRespectToThreshold(self):
        current_container = spatial_methods.NonHistorical.Conditions.NormMethods

        _, _, _, group_histogram, group_percentage_distribution, group_means, group_variances = current_container.Distribution(
            self.model_part, KratosCFD.Y_PLUS, "value", self.distribution_params)

        # % of element with y+ above threshold
        how_many = group_percentage_distribution[-1]*100.0

        # quantifying the mean and std for values below the threshold
        total_elements_in_threshold_range = group_histogram[0] + group_histogram[1]
        if (total_elements_in_threshold_range > 0):
            y_mean = (group_means[0] * group_histogram[0] + group_means[1] * group_histogram[1]) / total_elements_in_threshold_range

            threshold_sum_squared = (group_variances[0] + pow(group_means[0], 2)) * group_histogram[0] + (group_variances[1] + pow(group_means[1], 2)) * group_histogram[1]
            y_std = sqrt((threshold_sum_squared  - total_elements_in_threshold_range * pow(y_mean, 2)) / (total_elements_in_threshold_range - 1.0))
        else:
            y_mean = 0.0
            y_std = 0.0

        x_max, _ = current_container.Max(self.model_part, KratosCFD.Y_PLUS, "value")
        return [y_mean, y_std, x_max, how_many]

    def IsOutputStep(self):
        current_time = self.model_part.ProcessInfo[Kratos.TIME]
        current_step = self.model_part.ProcessInfo[Kratos.STEP]
        return (self.interval_utility.IsInInterval(current_time)) and (current_step % self.output_step == 0)