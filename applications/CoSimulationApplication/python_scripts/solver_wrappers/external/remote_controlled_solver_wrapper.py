# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import subprocess, os

def Create(settings, solver_name):
    return RemoteControlledSolverWrapper(settings, solver_name)

class RemoteControlledSolverWrapper(CoSimulationSolverWrapper):
    """Interface for the CFD-Solver TAU
    """
    def __init__(self, settings, solver_name):
        super().__init__(settings, solver_name)

        wrapper_settings = self.settings["solver_wrapper_settings"]

        start_external_solver = wrapper_settings["start_external_solver"].GetBool()
        if start_external_solver:
            command_txt = wrapper_settings["external_solver_start_command"].GetString()
            path_to_tau = wrapper_settings["path_to_tau"].GetString()
            # QUESTION: shall we use this:
            parent_path = os.path.join(os.path.dirname(__file__), '../..')
            # OR this:
            # parent_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
            path_to_tau_solver = parent_path + '/helpers/tau_solver.py'
            path_to_tau_solver = parent_path + '/helpers/TauSolver.py'
            tau_input_file = wrapper_settings["tau_input_file"].GetString()
            tau_log_file = wrapper_settings["tau_log_file"].GetString()
            # command_args = wrapper_settings["external_solver_arguments"].GetStringArray()
            command_args = [path_to_tau, path_to_tau_solver, tau_input_file, tau_log_file]
            cs_tools.cs_print_info(self._ClassName(), 'Running external solver with command: "{}" | arguments: "{}"'.format(command_txt, command_args))

            full_command = [command_txt]
            full_command.extend(command_args)
            # self.external_solver_process = subprocess.Popen(full_command, stderr=subprocess.PIPE, start_new_session=True) # TODO check what to use here
            # self.external_solver_process = subprocess.Popen(full_command)

        self.controlling_external_solver = wrapper_settings["controlling_external_solver"].GetBool()

        self.model_part_name = wrapper_settings["main_model_part_name"].GetString()
        cs_tools.CreateMainModelPartsFromCouplingData(self.data_dict.values(), self.model, self.name)
        cs_tools.AllocateHistoricalVariablesFromCouplingData(self.data_dict.values(), self.model, self.name)

    def Initialize(self):
        super().Initialize()

        interface_config = { "model_part_name" : self.model_part_name }

        # self.ImportCouplingInterface(interface_config)

    def AdvanceInTime(self, current_time):
        # self.__CheckExternalSolverProcess() # TODO check why this is blocking
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.AdvanceInTime)
            data_config = {
                "type"       : "time",
                "time"     : current_time
            }
            self.ExportData(data_config)
            self.ImportData(data_config)
            return 100.0 #data_config["time"]
        else:
            return 100.0

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.InitializeSolutionStep)

    def SolveSolutionStep(self):
        super().SolveSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.SolveSolutionStep)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.FinalizeSolutionStep)

    def Finalize(self):
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.BreakSolutionLoop)
        super().Finalize() # this also does the disconnect

    def ImportCouplingInterface(self, interface_config):
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ExportMesh, interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super().ImportCouplingInterface(interface_config)

    def ExportCouplingInterface(self, interface_config):
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ImportMesh, interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super().ExportCouplingInterface(interface_config)

    def ImportData(self, data_config):
        if self.controlling_external_solver and data_config["type"] == "coupling_interface_data":
            # CoSim imports, the external solver exports
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ExportData, data_config["interface_data"].name)
        super().ImportData(data_config)

    def ExportData(self, data_config):
        if self.controlling_external_solver:
            if data_config["type"] == "coupling_interface_data":
                # CoSim exports, the external solver imports
                self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ImportData, data_config["interface_data"].name)
            elif data_config["type"] == "convergence_signal":
                return # we control the ext solver, no need for sending a convergence signal
        super().ExportData(data_config)

    def PrintInfo(self):
        cs_tools.cs_print_info(self._ClassName(), "printing info...")

    def _GetIOType(self):
        return "kratos_co_sim_io"

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
