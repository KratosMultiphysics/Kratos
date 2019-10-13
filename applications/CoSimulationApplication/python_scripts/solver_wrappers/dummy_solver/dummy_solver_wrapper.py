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
    """This class serves as wrapper for the cpp/fortran dummy solver
    """
    def __init__(self, settings, solver_name):
        super(DummySolverWrapper, self).__init__(settings, solver_name)

        wrapper_settings = self.settings["solver_wrapper_settings"]

        start_external_solver = wrapper_settings["start_external_solver"].GetBool()
        if start_external_solver:
            command_txt = wrapper_settings["external_solver_start_command"].GetString()
            command_args = wrapper_settings["external_solver_arguments"].GetStringArray()
            cs_tools.cs_print_info(self._ClassName(), 'Running external solver with command: "{}" | arguments: "{}"'.format(command_txt, command_args))

            full_command = [command_txt]
            full_command.extend(command_args)
            self.external_solver_process = subprocess.Popen(full_command, stderr=subprocess.PIPE, start_new_session=True)

        self.controlling_external_solver = wrapper_settings["controlling_external_solver"].GetBool()

        self.model_part_name = self.settings["solver_wrapper_settings"]["main_model_part_name"].GetString()
        cs_tools.CreateMainModelPartsFromCouplingData(self.data_dict.values(), self.model, self.name)
        cs_tools.AllocateHistoricalVariablesFromCouplingData(self.data_dict.values(), self.model, self.name)

    def Initialize(self):
        super(DummySolverWrapper, self).Initialize()

        interface_config = { "model_part_name" : self.model_part_name }

        self.ImportCouplingInterface(interface_config)

    def AdvanceInTime(self, current_time):
        self.__CheckExternalSolverProcess()
        if self.controlling_external_solver:
            self.__SendControlSignal("AdvanceInTime")
            # TODO this requires more, then delete the next line!
            return 0.0
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
        if self.controlling_external_solver:
            self.__SendControlSignal("Finalize")
        super(DummySolverWrapper, self).Finalize() # this also does the disconnect

    def ImportCouplingInterface(self, interface_config):
        if self.controlling_external_solver:
            self.__SendControlSignal("ExportMesh", interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super(DummySolverWrapper, self).ImportCouplingInterface(interface_config)

    def ExportCouplingInterface(self, interface_config):
        if self.controlling_external_solver:
            self.__SendControlSignal("ImportMesh", interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super(DummySolverWrapper, self).ExportCouplingInterface(interface_config)

    def ImportData(self, data_config):
        if self.controlling_external_solver:
            self.__SendControlSignal("ExportData", data_config["interface_data"].name)
        super(DummySolverWrapper, self).ImportData(data_config)

    def ExportData(self, data_config):
        if not data_config["type"] == "control_signal" and self.controlling_external_solver:
            self.__SendControlSignal("ImportData", data_config["interface_data"].name)
        super(DummySolverWrapper, self).ExportData(data_config)

    def PrintInfo(self):
        cs_tools.cs_print_info(self._ClassName(), "printing info...")

    def _GetIOType(self):
        return "dummy_solver.dummy_solver_io"

    def __SendControlSignal(self, signal, identifier=""):
        data_config = {
            "type"       : "control_signal",
            "signal"     : signal,
            "identifier" : identifier
        }
        self.ExportData(data_config)

    def __CheckExternalSolverProcess(self):
        if hasattr(self, 'external_solver_process') and self.external_solver_process.poll() is None:
            _, process_stderr = self.external_solver_process.communicate()
            if process_stderr:
                raise Exception("{} terminated with the following error:\n{}".format(self._ClassName(), process_stderr.decode('ascii')))
