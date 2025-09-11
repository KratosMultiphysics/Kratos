# Importing the Kratos Library
import KratosMultiphysics
# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeReactionProcess(model, settings["Parameters"])

class ComputeReactionProcess(KratosMultiphysics.OutputProcess):
    """
    Auxiliary base class to output total reaction of over a line or surface.
    A derived class needs to be implemented to be able to use this functionality
    (in both conforming and non-conforming boundary conditions),
    as calling the base class alone is not enough.
    """
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:

        super().__init__()

        self.params = settings
        self.params.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.model_part_name = self.params["model_part_name"].GetString()
        if self.model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        elif self.model_part_name.startswith('Background_Grid'):
            self.model_part_name = self.model_part_name.replace('Background_Grid','MPM_Material')

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "",
                "interval"                  : [0.0, 1e30],
                "output_interval"           : 1.0,
                "output_control_type"       : "step",
                "print_format"              : ".8f",
                "output_file_settings"      : {}
            }
            """)

    def ExecuteBeforeSolutionLoop(self):
        self.model_part = self.model[self.model_part_name]

        self.format = self.params["print_format"].GetString()
        self.interval = KratosMultiphysics.IntervalUtility(self.params)
        self.controller = KratosMultiphysics.OutputController(self.model, self.params)

        file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])
        output_file_name = self.params["model_part_name"].GetString() + "_reaction.dat"
        if file_handler_params.Has("file_name"):
            warn_msg  = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
            warn_msg += '"' + file_handler_params["file_name"].GetString() + '"}\n'
            warn_msg += 'Using this specified file name instead of the default "' + output_file_name + '"'
            KratosMultiphysics.Logger.PrintWarning("ComputeReactionProcess", warn_msg)
        else:
            file_handler_params.AddEmptyValue("file_name")
            file_handler_params["file_name"].SetString(output_file_name)
        self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_params, self._GetFileHeader()).file
        self.output_file.flush()

    def PrintOutput(self) -> None:
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        reaction_force = self._GetReaction() # Compute reaction force
        self.output_file.write(f"{str(current_time)} {reaction_force[0]:{self.format}} {reaction_force[1]:{self.format}} {reaction_force[2]:{self.format}}\n")
        self.output_file.flush()
        self.controller.Update()

    def IsOutputStep(self) -> bool:
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        return self.interval.IsInInterval(current_time) and self.controller.Evaluate()

    def ExecuteFinalize(self):
        self.output_file.close()

    def _GetFileHeader(self):
        err_msg  = 'ComputeReactionProcess: _GetFileHeader called in base class\n'
        err_msg += 'this needs to be implemented and called from derived class'
        raise Exception(err_msg)

    def _GetReaction(self):
        err_msg  = 'ComputeReactionProcess: _GetReaction called in base class\n'
        err_msg += 'this needs to be implemented and called from derived class'
        raise Exception(err_msg)
