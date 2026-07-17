from pathlib import Path

import KratosMultiphysics


def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MeshioOutputProcess(model, settings["Parameters"])


class MeshioOutputProcess(KratosMultiphysics.OutputProcess):
    """Output process writing results in any meshio++-supported format.

    It wraps KratosMultiphysics.MeshioPlusPlusIO: every output step calls
    WriteModelPart, which extends the current output instead of overwriting it
    (an XDMF temporal collection in a single file, or a file series
    <output_name>_<label>.<ext> for the other formats).

    The format is taken from the "format" setting, or resolved from the
    extension of "output_name" when it is "auto". Query the formats available
    in this build with KratosMultiphysics.MeshioPlusPlusIO.GetSupportedWriteFormats().
    """

    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        super().__init__()

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[settings["model_part_name"].GetString()]

        output_name = settings["output_name"].GetString()
        if not output_name:
            output_name = self.model_part.Name
        if settings["format"].GetString() in ("auto", "automatic", "") and not Path(output_name).suffix:
            raise Exception('"output_name" needs a file extension when "format" is "auto" '
                            '(e.g. "results.vtu"), or set "format" explicitly.')

        output_path = Path(settings["output_path"].GetString())
        if settings["save_output_files_in_folder"].GetBool():
            output_path.mkdir(parents=True, exist_ok=True)
            file_name = output_path / output_name
        else:
            file_name = Path(output_name)

        # Forward the IO-specific settings; defaults are defined in
        # MeshioPlusPlusIO::GetDefaultParameters (meshioplusplus_io.cpp)
        io_settings = KratosMultiphysics.Parameters("""{}""")
        for key in ("format", "file_format", "time_series", "output_control_type",
                    "output_precision", "label_precision",
                    "custom_name_prefix", "custom_name_postfix",
                    "entity_type", "output_sub_model_parts",
                    "write_deformed_configuration", "write_ids", "xdmf_data_format",
                    "nodal_solution_step_data_variables", "nodal_data_value_variables",
                    "nodal_flags",
                    "element_data_value_variables", "element_flags",
                    "condition_data_value_variables", "condition_flags",
                    "gauss_point_variables_extrapolated_to_nodes",
                    "gauss_point_variables_in_elements"):
            io_settings.AddValue(key, settings[key])

        self.meshio_io = KratosMultiphysics.MeshioPlusPlusIO(str(file_name), io_settings)

        self.__controller = KratosMultiphysics.OutputController(model, settings)

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_name"                                 : "",
            "output_path"                                 : "meshio_output",
            "save_output_files_in_folder"                 : true,
            "format"                                      : "auto",
            "file_format"                                 : "default",
            "time_series"                                 : "automatic",
            "output_control_type"                         : "step",
            "output_interval"                             : 1.0,
            "output_precision"                            : 7,
            "label_precision"                             : 4,
            "custom_name_prefix"                          : "",
            "custom_name_postfix"                         : "",
            "entity_type"                                 : "automatic",
            "output_sub_model_parts"                      : false,
            "write_deformed_configuration"                : false,
            "write_ids"                                   : false,
            "xdmf_data_format"                            : "auto",
            "nodal_solution_step_data_variables"          : [],
            "nodal_data_value_variables"                  : [],
            "nodal_flags"                                 : [],
            "element_data_value_variables"                : [],
            "element_flags"                               : [],
            "condition_data_value_variables"              : [],
            "condition_flags"                             : [],
            "gauss_point_variables_extrapolated_to_nodes" : [],
            "gauss_point_variables_in_elements"           : []
        }""")

    def Check(self) -> int:
        return self.__controller.Check()

    def IsOutputStep(self) -> bool:
        return self.__controller.Evaluate()

    def PrintOutput(self) -> None:
        # The IO keeps the transient state: this extends the existing output
        self.meshio_io.WriteModelPart(self.model_part)
        self.__controller.Update()
