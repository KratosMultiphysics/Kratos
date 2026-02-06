import KratosMultiphysics
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMReactionOutputProcess(model, settings["Parameters"])

class MPMReactionOutputProcess(KratosMultiphysics.OutputProcess):
    """
    Auxiliary base class to output total reaction of over a line or surface.
    A derived class needs to be implemented to be able to use this functionality
    (in both conforming and non-conforming boundary conditions),
    as calling the base class alone is not enough.
    """
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        super().__init__()

        self.model = model
        self.params = settings
        self.params.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.output_file = None
        self.controller = KratosMultiphysics.OutputController(self.model, self.params)
        self.format = self.params["print_format"].GetString()

        self.model_part_name = self.params["model_part_name"].GetString()
        if not self.model_part_name:
            raise Exception('No "model_part_name" was specified!')
        elif self.model_part_name.startswith('Background_Grid'):
            self.model_part_name = self.model_part_name.replace('Background_Grid','MPM_Material')

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "",
                "output_interval"           : 1.0,
                "output_control_type"       : "step",
                "print_format"              : ".8f",
                "output_file_settings"      : {}
            }
            """)

    def ExecuteBeforeSolutionLoop(self) -> None:
        self.model_part = self.model[self.model_part_name]
        file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])

        # If output file name is not specified, use default value: <model_part_name>_reaction.dat
        if not file_handler_params.Has("file_name"):
            output_file_name = f"{self.model_part_name}_reaction"
            wrn_msg = f"File name not specified, using the default name: {output_file_name}"
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__,wrn_msg)
            file_handler_params.AddString("file_name",output_file_name)

        self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_params, self._GetFileHeader()).file
        self.PrintOutput()

    def PrintOutput(self) -> None:
        if self.output_file:
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            reaction_force = self._GetReaction() # Compute reaction force
            out = f"{str(time)} {reaction_force[0]:{self.format}} {reaction_force[1]:{self.format}} {reaction_force[2]:{self.format}}\n"
            self.output_file.write(out)
            self.output_file.flush()
            self.controller.Update()

    def IsOutputStep(self) -> bool:
        return self.controller.Evaluate()

    def ExecuteFinalize(self) -> None:
        if self.output_file:
            self.output_file.close()

    def _GetFileHeader(self):
        err_msg  = '{self.__class__.__name__}: _GetFileHeader called in base class\n'
        err_msg += 'this needs to be implemented and called from derived class'
        raise Exception(err_msg)

    def _GetReaction(self):
        err_msg  = '{self.__class__.__name__}: _GetReaction called in base class\n'
        err_msg += 'this needs to be implemented and called from derived class'
        raise Exception(err_msg)
