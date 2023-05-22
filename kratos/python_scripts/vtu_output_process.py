from pathlib import Path
import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils

def Factory(parameters: Kratos.Parameters, model: Kratos.Model):
    if not isinstance(parameters, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, Kratos.Model):
        raise Exception("expected input shall be a model object")

    return VtuOutputProcess(model, parameters["Parameters"])


class VtuOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""
        {
            "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "file_format"                                 : "binary",
            "output_precision"                            : 7,
            "output_control_type"                         : "step",
            "output_interval"                             : 1.0,
            "output_sub_model_parts"                      : false,
            "output_path"                                 : "VTU_Output",
            "custom_name_prefix"                          : "",
            "save_output_files_in_folder"                 : true,
            "write_deformed_configuration"                : false,
            "nodal_solution_step_data_variables"          : [],
            "nodal_data_value_variables"                  : [],
            "nodal_flags"                                 : [],
            "element_data_value_variables"                : [],
            "element_flags"                               : [],
            "condition_data_value_variables"              : [],
            "condition_flags"                             : []

        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__()

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[parameters["model_part_name"].GetString()]
        self.write_deformed_configuration = parameters["write_deformed_configuration"].GetBool()
        self.output_precision = parameters["output_precision"].GetInt()
        self.output_sub_model_parts = parameters["output_sub_model_parts"].GetBool()
        self.custom_name_prefix = parameters["custom_name_prefix"].GetString()

        file_format = parameters["file_format"].GetString()
        if file_format == "ascii":
            self.writer_format = Kratos.VtuOutput.ASCII
        elif file_format == "binary":
            self.writer_format = Kratos.VtuOutput.BINARY
        else:
            raise RuntimeError(f"Unsupported file format requested [ requested format = {file_format} ]. Supported file formats:\n\tascii\n\tbinary")

        if parameters["save_output_files_in_folder"].GetBool():
            output_path = parameters["output_path"].GetString()
            if not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
                kratos_utils.DeleteDirectoryIfExisting(output_path)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
            # now create the output path
            p_output_path = Path(output_path)
            p_output_path.mkdir(parents=True, exist_ok=True)
            self.p_output_file_name_prefix = p_output_path / self.custom_name_prefix
        else:
            self.p_output_file_name_prefix = Path(self.custom_name_prefix)


        self.output_interval = parameters["output_interval"].GetDouble()
        self.output_control = parameters["output_control_type"].GetString()
        self.next_output = 0.0

        self.vtu_output_ios: 'list[Kratos.VtuOutput]' = []

        self.__ScheduleNextOutput() # required here esp for restart

    def ExecuteInitialize(self) -> None:
        # check and create all the vtu outputs
        self.vtu_output_ios.append(Kratos.VtuOutput(self.model_part, not self.write_deformed_configuration, self.writer_format, self.output_precision))


        if self.output_sub_model_parts:
            for sub_model_part in self.model_part.SubModelParts():
                self.vtu_output_ios.append(Kratos.VtuOutput(sub_model_part, not self.write_deformed_configuration, self.writer_format, self.output_precision))

    def PrintOutput(self):
        for vtu_output in self.vtu_output_ios:
            vtu_output.PrintOutput(str(self.p_output_file_name_prefix) + vtu_output.GetModelPart().FullName())

        self.__ScheduleNextOutput()

    def IsOutputStep(self):
        if self.output_control == "time":
            return self.__GetTime() >= self.next_output
        else:
            return self.model_part.ProcessInfo[Kratos.STEP] >= self.next_output

    def __ScheduleNextOutput(self):
        if self.output_interval > 0.0: # Note: if == 0, we'll just always print
            if self.output_control == "time":
                while self.next_output <= self.__GetTime():
                    self.next_output += self.output_interval
            else:
                while self.next_output <= self.model_part.ProcessInfo[Kratos.STEP]:
                    self.next_output += self.output_interval

    def __GetTime(self):
        # remove rounding errors that mess with the comparison
        # e.g. 1.99999999999999999 => 2.0
        return float("{0:.12g}".format(self.model_part.ProcessInfo[Kratos.TIME]))
