# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

from statistics import mean, stdev


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

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])) and (current_step % self.output_step == 0):
            cfl_value = self._EvaluateCFL()

            if (self.model_part.GetCommunicator().MyPID() == 0):
                output = self._SummarizeCFL(cfl_value)
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

    def _CalculateWithRespectToThreshold(self, x):

        y = [val for val in x if val < self.cfl_output_limit]
        y1 = [val for val in x if val < 1.0]
        # % of element with cfl above threshold
        how_many = ((len(x)-len(y))/len(x))*100
        # % of element with cfl above 1
        how_many1 = ((len(x)-len(y1))/len(x))*100

        # quantifying the mean and std for values below the threshold
        y_mean = mean(y)
        y_std = stdev(y)

        # qunatifying the global max
        x_max = max(x)
        return [y_mean, y_std, x_max, how_many, how_many1]

    def _EvaluateCFL(self):

        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            KratosCFD.EstimateDtUtility2D.CalculateLocalCFL(self.model_part)
        else:
            KratosCFD.EstimateDtUtility3D.CalculateLocalCFL(self.model_part)

        local_cfl = []
        for elem in self.model_part.Elements:
            local_cfl.append(elem.GetValue(KratosMultiphysics.CFL_NUMBER))

        local_cfl = self.model_part.GetCommunicator().GetDataCommunicator().GathervDoubles(local_cfl, 0)

        return local_cfl

    def _SummarizeCFL(self, local_cfl):

        global_cfl = []
        for k in local_cfl:
            global_cfl.extend(k)

        cfl_mean, cfl_std, cfl_max, cfl_how_many, cfl_how_many1 = self._CalculateWithRespectToThreshold(
            global_cfl)

        return [cfl_mean, cfl_std, cfl_max, cfl_how_many, cfl_how_many1]
