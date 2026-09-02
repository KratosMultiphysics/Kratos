from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.pyvista_utilities as pyvista_utilities


def Factory(parameters: Kratos.Parameters, model: Kratos.Model) -> Kratos.OutputProcess:
    if not isinstance(parameters, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, Kratos.Model):
        raise Exception("expected input shall be a model object")

    return PyVistaAnimationOutputProcess(model, parameters["Parameters"])


class PyVistaAnimationOutputProcess(Kratos.OutputProcess):
    """Record the solution loop as an animation using PyVista.

    One frame is rendered every time the OutputController fires, and the animation
    is encoded in ExecuteFinalize. The output format follows the "file_name"
    suffix: ".gif" needs the optional 'imageio' package, movie suffixes such as
    ".mp4" additionally need 'imageio-ffmpeg', and ".png" writes a zero-padded
    frame sequence that needs nothing beyond PyVista.

    Only rank 0 renders; the other ranks are silent no-ops.
    """

    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""
        {
            "model_part_name"                   : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_control_type"               : "step",
            "output_interval"                   : 1.0,
            "output_path"                       : "PyVista_Animation",
            "file_name"                         : "animation.gif",
            "echo_level"                        : 0,
            "fps"                               : 10.0,
            "quality"                           : 5,
            "loop"                              : 0,
            "window_size"                       : [1024, 768],
            "variable_name"                     : "",
            "component"                         : 0,
            "label"                             : "",
            "cmap"                              : "turbo",
            "clim"                              : [],
            "lock_clim"                         : false,
            "warp_by_vector"                    : "",
            "warp_factor"                       : 1.0,
            "show_undeformed"                   : true,
            "show_edges"                        : true,
            "use_deformed_configuration"        : false,
            "view"                              : "default",
            "theme"                             : "",
            "time_annotation"                   : true,
            "time_format"                       : "Time = {:.4g}",
            "nodal_solution_step_data_variables": [],
            "nodal_data_value_variables"        : [],
            "element_data_value_variables"      : [],
            "condition_data_value_variables"    : []
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        super().__init__()

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[parameters["model_part_name"].GetString()]
        self.echo_level = parameters["echo_level"].GetInt()

        # Only rank 0 renders: VTK rendering is not distributed, and gathering the
        # whole distributed mesh onto one rank every output step is not worth it.
        self.is_writing_rank = self.model_part.GetCommunicator().MyPID() == 0

        self.output_path = Path(parameters["output_path"].GetString())
        if self.is_writing_rank and not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
            kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))
        self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
        if self.is_writing_rank:
            self.output_path.mkdir(parents=True, exist_ok=True)

        self.file_name = str(self.output_path / parameters["file_name"].GetString())

        self.recorder_settings = self.__CreateRecorderSettings(parameters)
        self.recorder = None

        self.__controller = Kratos.OutputController(model, parameters)

    def __CreateRecorderSettings(self, parameters: Kratos.Parameters) -> dict:
        """Translate the JSON parameters into TransientPlotter keyword arguments."""
        def optional_string(key):
            value = parameters[key].GetString()
            return value if value != "" else None

        clim = parameters["clim"].GetVector()
        nodal_variables = (parameters["nodal_solution_step_data_variables"].GetStringArray()
                           + parameters["nodal_data_value_variables"].GetStringArray())

        return {
            "fps"                     : parameters["fps"].GetDouble(),
            "quality"                 : parameters["quality"].GetInt(),
            "loop"                    : parameters["loop"].GetInt(),
            "windowSize"              : [int(v) for v in parameters["window_size"].GetVector()],
            "view"                    : parameters["view"].GetString(),
            "theme"                   : optional_string("theme"),
            "name"                    : optional_string("variable_name"),
            "component"               : parameters["component"].GetInt(),
            "label"                   : optional_string("label"),
            "cmap"                    : parameters["cmap"].GetString(),
            "clim"                    : (float(clim[0]), float(clim[1])) if len(clim) == 2 else None,
            "lockClim"                : parameters["lock_clim"].GetBool(),
            "warpByVector"            : optional_string("warp_by_vector"),
            "factor"                  : parameters["warp_factor"].GetDouble(),
            "showUndeformed"          : parameters["show_undeformed"].GetBool(),
            "showEdges"               : parameters["show_edges"].GetBool(),
            "useDeformedConfiguration": parameters["use_deformed_configuration"].GetBool(),
            "timeAnnotation"          : parameters["time_annotation"].GetBool(),
            "timeFormat"              : parameters["time_format"].GetString(),
            "nodalVariables"          : nodal_variables,
            "elementVariables"        : parameters["element_data_value_variables"].GetStringArray(),
            "conditionVariables"      : parameters["condition_data_value_variables"].GetStringArray(),
        }

    def ExecuteInitialize(self) -> None:
        if not self.is_writing_rank:
            return

        self.recorder = pyvista_utilities.TransientPlotter(
            self.file_name, **self.recorder_settings
        )

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(
                "PyVistaAnimationOutputProcess",
                f"Recording animation to '{self.file_name}'."
            )

    def Check(self) -> int:
        return self.__controller.Check()

    def IsOutputStep(self) -> bool:
        return self.__controller.Evaluate()

    def PrintOutput(self) -> None:
        if self.recorder is not None:
            self.recorder.AddFrame(self.model_part)
        self.__controller.Update()

    def ExecuteFinalize(self) -> None:
        if self.recorder is None:
            return

        number_of_frames = self.recorder.Close()
        self.recorder = None

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(
                "PyVistaAnimationOutputProcess",
                f"Wrote '{self.file_name}' with {number_of_frames} frames."
            )
