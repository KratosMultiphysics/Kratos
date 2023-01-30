# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# Additional imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Create(*args):
    return PrintIterationNumberOperation(*args)

class PrintIterationNumberOperation(CoSimulationCouplingOperation):
    """This operation is used to print the number of iterations on the strong coupling schemes
    TODO:
    - add tests
    - more cleanup
    """
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = self.model[self.model_part_name]

        self.interval = KM.IntervalUtility(settings)

        if self.model_part.GetCommunicator().MyPID() == 0:
            output_file_name = self.model_part_name + "_number_iterations.dat"
            file_handler_settings = KM.Parameters(self.settings["output_file_settings"])
            if file_handler_settings.Has("file_name"):
                warn_msg = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                warn_msg += '"' + file_handler_settings["file_name"].GetString() + '"}\n'
                warn_msg += 'Using this specififed file name instead of the default "' + output_file_name + '"'
                cs_tools.cs_print_info(self._ClassName(), warn_msg)
            else:
                file_handler_settings.AddEmptyValue("file_name")
                file_handler_settings["file_name"].SetString(output_file_name)
            file_header = self._GetFileHeader()
            self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_settings, file_header).file

    def InitializeSolutionStep(self):
        self.iteration_number = 0

    def InitializeCouplingIteration(self):
        self.iteration_number += 1

    def FinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KM.TIME]
        if self.interval.IsInInterval(current_time):
            if self.model_part.GetCommunicator().MyPID() == 0:
                self.output_file.write(str(current_time) + "\t" + str(self.iteration_number) + "\n")

    def Finalize(self):
        if self.model_part.GetCommunicator().MyPID() == 0:
            self.output_file.close()

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    def _GetFileHeader(self):
        header = "#TIME[s]" + "\t" + "ITERATION_NUMBER\n"
        return header

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : "",
            "interval"              : [0.0, 1e30],
            "output_file_settings"  : {}
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
