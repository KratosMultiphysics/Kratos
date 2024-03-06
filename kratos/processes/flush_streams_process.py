# --- Kratos Imports ---
import KratosMultiphysics

# --- STD Imports ---
import sys
import time


class FlushStreamsProcess(KratosMultiphysics.Process):
    """ @brief Flush output streams at a given frequency.
        @param interval interval between forced flushes in seconds.
        @param stdout controls whether stdout is flushed.
        @param stderr controls whether stderr is flushed.
        @param stages list of process stages to flush at. Options are
                      - Initialize
                      - BeforeSolutionLoop
                      - InitializeSolutionStep
                      - FinalizeSolutionStep
                      - BeforeOutputStep
                      - AfterOutputStep
                      - Finalize
        @ingroup KratosCore
    """

    def __init__(self,
                 _: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters):
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__buffer_interval: float = parameters["interval"].GetDouble()
        self.__flush_stdout: bool = parameters["stdout"].GetBool()
        self.__flush_stderr: bool = parameters["stderr"].GetBool()
        self.__stages: "set[str]" = set()
        self.__last_flush: float = time.time()

        # Parse trigger stages
        admissible_stages: "set[str]" = set(["Initialize",
                                             "BeforeSolutionLoop",
                                             "InitializeSolutionStep",
                                             "FinalizeSolutionStep",
                                             "BeforeOutputStep",
                                             "AfterOutputStep",
                                             "Finalize"])

        for stage_parameter in parameters["stages"].values():
            stage_name = stage_parameter.GetString()
            if stage_name not in admissible_stages:
                raise ValueError(f"invalid stage {stage_name}. Options: {' '.join(admissible_stages)}")
            self.__stages.add(stage_name)

    def Execute(self) -> None:
        current_time = time.time()
        if self.__buffer_interval <= current_time - self.__last_flush:
            self.__Flush()
            self.__last_flush = current_time

    def ExecuteInitialize(self) -> None:
        if "Initialize" in self.__stages:
            self.Execute()

    def ExecuteBeforeSolutionLoop(self) -> None:
        if "BeforeSolutionLoop" in self.__stages:
            self.Execute()

    def ExecuteInitializeSolutionStep(self) -> None:
        if "InitializeSolutionStep" in self.__stages:
            self.Execute()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if "FinalizeSolutionStep" in self.__stages:
            self.Execute()

    def ExecuteBeforeOutputStep(self) -> None:
        if "BeforeOutputStep" in self.__stages:
            self.Execute()

    def ExecuteAfterOutputStep(self) -> None:
        if "AfterOutputStep" in self.__stages:
            self.Execute()

    def ExecuteFinalize(self) -> None:
        if "Finalize" in self.__stages:
            self.Execute()

    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "interval" : 0.0,
            "stdout" : true,
            "stderr" : true,
            "stages" : ["FinalizeSolutionStep"]
        }""")

    def __Flush(self) -> None:
        if self.__flush_stdout:
            sys.stdout.flush()
        if self.__flush_stderr:
            sys.stderr.flush()


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> FlushStreamsProcess:
    return FlushStreamsProcess(model, parameters["Parameters"])
