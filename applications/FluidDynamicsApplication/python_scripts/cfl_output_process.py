# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

if CheckIfApplicationsAvailable("StatisticsApplication"):
    from KratosMultiphysics.StatisticsApplication import SpatialMethods as spatial_methods
else:
    msg = "CFLOutputProcess requires StatisticsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)

from math import sqrt


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return CFLOutputProcess(model, settings["Parameters"])


class CFLOutputProcess(KratosMultiphysics.Process):
    """
    A class responsible for the CFL output, which is an element value in Kratos.
    """

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "",
                "interval"             : [0.0, 1e30],
                "cfl_output_limit"     : 2.5,
                "print_to_screen"      : false,
                "print_format"         : ".8f",
                "write_output_file"    : true,
                "output_step"          : 8,
                "output_file_settings": {}
            }
            """)


        # Detect "End" as a tag and replace it by a large number
        if(params.Has("interval")):
            if(params["interval"][1].IsString()):
                if(params["interval"][1].GetString() == "End"):
                    params["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:" +
                                    params["interval"].PrettyPrintJsonString())

        params.ValidateAndAssignDefaults(default_settings)

        # getting the ModelPart from the Model
        self.model_part_name = params["model_part_name"].GetString()
        if self.model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.model_part_name]

        self.interval = params["interval"].GetVector()

        # getting output limit for summarization
        self.cfl_output_limit = params["cfl_output_limit"].GetDouble()

        # TODO: Is it ok to do this check? If not, distribution calculation is going to be messy with if conditions for
        #       case with cfl_output_limit <= 1.0
        if (self.cfl_output_limit <= 1.0):
            raise Exception("Please provide cfl_output_limit greater than 1.0")

        self.format = params["print_format"].GetString()
        self.output_step = params["output_step"].GetInt()
        self.print_to_screen = params["print_to_screen"].GetBool()
        self.write_output_file = params["write_output_file"].GetBool()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):
                file_handler_params = KratosMultiphysics.Parameters(
                    params["output_file_settings"])

                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                   file_handler_params, file_header).file

        self.distribution_params = KratosMultiphysics.Parameters('''{
            "number_of_value_groups" : 1,
            "min_value"              : "min",
            "max_value"              : "max"
        }''')
        self.distribution_params["min_value"].SetDouble(min(self.cfl_output_limit, 1.0))
        self.distribution_params["max_value"].SetDouble(max(self.cfl_output_limit, 1.0))

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])) and (current_step % self.output_step == 0):
            self._EvaluateCFL()
            output = self._CalculateWithRespectToThreshold()

            if (self.model_part.GetCommunicator().MyPID() == 0):
                output_vals = [format(val, self.format) for val in output]

                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                output_vals.insert(0, str(current_time))

                res_labels = ["time: ", "mean: ", "std: ", "max: ", "cfl" +
                              "{:.1f}".format(self.cfl_output_limit) + ": ", "cfl1.0: "]

                if (self.print_to_screen):

                    result_msg = 'CFL evaluation for model part ' + \
                        self.model_part_name + '\n'
                    result_msg += ', '.join([a+b for a,
                                             b in zip(res_labels, output_vals)])
                    self._PrintToScreen(result_msg)

                if (self.write_output_file):
                    self.output_file.write(' '.join(output_vals) + "\n")

    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

    def _GetFileHeader(self):
        header = '# CFL for model part ' + self.model_part_name + \
            '| CFL_threshold: ' + str(self.cfl_output_limit) + '\n'
        header += '# Time Mean Std Max HowMany>' + \
            "{:.1f}".format(self.cfl_output_limit) + ' [%] HowMany>1.0 [%]\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            "CFLOutputProcess", "CFL VALUE RESULTS:")
        KratosMultiphysics.Logger.PrintInfo(
            "CFLOutputProcess", "Current time: " + result_msg)

    def _CalculateWithRespectToThreshold(self):
        current_container = spatial_methods.NonHistorical.Elements.NormMethods

        _, _, _, group_histogram, group_percentage_distribution, group_means, group_variances = current_container.Distribution(
            self.model_part, KratosMultiphysics.CFL_NUMBER, "value", self.distribution_params)

        # % of element with cfl above threshold
        how_many = group_percentage_distribution[-1]*100.0
        # % of element with cfl above 1
        how_many1 = (1.0 - group_percentage_distribution[0])*100.0

        # quantifying the mean and std for values below the threshold
        total_elements_in_threshold_range = group_histogram[0] + group_histogram[1]
        if (total_elements_in_threshold_range > 0):
            y_mean = (group_means[0] * group_histogram[0] + group_means[1] * group_histogram[1]) / total_elements_in_threshold_range

            threshold_sum_squared = (group_variances[0] + pow(group_means[0], 2)) * group_histogram[0] + (group_variances[1] + pow(group_means[1], 2)) * group_histogram[1]
            y_std = sqrt((threshold_sum_squared  - total_elements_in_threshold_range * pow(y_mean, 2)) / (total_elements_in_threshold_range - 1.0))
        else:
            y_mean = 0.0
            y_std = 0.0

        # qunatifying the global max
        # TODO: @Mate, where should we put the id of the element, where max is (second bland output argument is max_id)?
        x_max, _ = current_container.Max(self.model_part, KratosMultiphysics.CFL_NUMBER, "value")
        return [y_mean, y_std, x_max, how_many, how_many1]

    def _EvaluateCFL(self):
        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            KratosCFD.EstimateDtUtility2D.CalculateLocalCFL(self.model_part)
        else:
            KratosCFD.EstimateDtUtility3D.CalculateLocalCFL(self.model_part)
