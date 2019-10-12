from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import subprocess

def Create(settings, solver_name):
    return DummySolverWrapper(settings, solver_name)

class DummySolverWrapper(CoSimulationSolverWrapper):
    """This class serves as wrapper for the cpp ping and pong solvers
    """
    def __init__(self, settings, solver_name):
        super(DummySolverWrapper, self).__init__(settings, solver_name)

        wrapper_settings = self.settings["solver_wrapper_settings"]

        start_external_solver = wrapper_settings["start_external_solver"].GetBool()
        if start_external_solver:
            command_txt = wrapper_settings["external_solver_start_command"].GetString()
            print("Running : ", command_txt)
            self.rv = subprocess.Popen(command_txt, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True, start_new_session=True) # TODO check this, esp it should also block the caller if it hangs (basically propagate all errors etc)

        self.controlling_external_solver = wrapper_settings["controlling_external_solver"].GetBool()

    def AdvanceInTime(self, current_time):
        if self.controlling_external_solver:
            self.__SendControlSignal("AdvanceInTime")
            # TODO this requires more!
        else:
            return 0.0

    def InitializeSolutionStep(self):
        super(DummySolverWrapper, self).InitializeSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal("InitializeSolutionStep")

    def SolveSolutionStep(self):
        super(DummySolverWrapper, self).SolveSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal("SolveSolutionStep")

    def FinalizeSolutionStep(self):
        super(DummySolverWrapper, self).FinalizeSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal("FinalizeSolutionStep")

    def Finalize(self):
        super(DummySolverWrapper, self).Finalize()
        if self.controlling_external_solver:
            self.__SendControlSignal("Finalize")

    def ImportCouplingInterface(self, interface_config):
        if self.controlling_external_solver:
            self.__SendControlSignal("ExportMesh") # TODO this can also be geometry at some point
        super(DummySolverWrapper, self).ImportCouplingInterface(interface_config)

    def ExportCouplingInterface(self, interface_config):
        if self.controlling_external_solver:
            self.__SendControlSignal("ImportMesh") # TODO this can also be geometry at some point
        super(DummySolverWrapper, self).ExportCouplingInterface(interface_config)

    def ImportData(self, data_config):
        if self.controlling_external_solver:
            self.__SendControlSignal("ExportData")
        super(DummySolverWrapper, self).ImportData(data_config)

    def ExportData(self, data_config):
        if self.controlling_external_solver:
            self.__SendControlSignal("ImportData")
        super(DummySolverWrapper, self).ExportData(data_config)

    def PrintInfo(self):
        cs_tools.cs_print_info(self._ClassName(), "printing info...")

    def _GetIOType(self):
        return "dummy_solver_io"

    def __SendControlSignal(self, signal):
        data_config = {
            "type" : "control_signal",
            "signal" : signal
        }
        self.ExportData(data_config)
