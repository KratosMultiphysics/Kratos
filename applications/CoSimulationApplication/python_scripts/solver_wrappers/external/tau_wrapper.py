from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import subprocess, os

def Create(settings, solver_name):
    return TAUWrapper(settings, solver_name)

class TAUWrapper(CoSimulationSolverWrapper):
    """Interface for the CFD-Solver TAU
    """
    def __init__(self, settings, solver_name):
        super(TAUWrapper, self).__init__(settings, solver_name)
        print("Hellow world wrapper")

        wrapper_settings = self.settings["solver_wrapper_settings"]
        self.coupling_interface_imported = False

        start_external_solver = wrapper_settings["start_external_solver"].GetBool()
        print("start_external_solver = ", start_external_solver)
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
            self.external_solver_process = subprocess.Popen(full_command, stderr=subprocess.PIPE, start_new_session=True) # TODO check what to use here
            # self.external_solver_process = subprocess.Popen(full_command)

        self.controlling_external_solver = wrapper_settings["controlling_external_solver"].GetBool()

        self.model_part_name = wrapper_settings["main_model_part_name"].GetString()
        # cs_tools.CreateMainModelPartsFromCouplingData(self.data_dict.values(), self.model, self.name)
        cs_tools.CreateModelPartsFromCouplingData(self.data_dict.values(), self.model, self.name)
        cs_tools.AllocateHistoricalVariablesFromCouplingData(self.data_dict.values(), self.model, self.name)

        self.time = wrapper_settings["start_time"].GetDouble()
        self.time_step = wrapper_settings["time_step"].GetDouble()

    def Initialize(self):
        print('TAUWrapper Initialize')
        super(TAUWrapper, self).Initialize()

        print('TAUWrapper Initialize')

    def AdvanceInTime(self, current_time):
        # self.__CheckExternalSolverProcess() # TODO check why this is blocking
        if self.controlling_external_solver:
            print('tau_wrapper advanceintime')
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.AdvanceInTime)
            data_config = {
                "type"       : "time",
                "time"     : current_time
            }
            print('tau_wrapper ExportData')
            self.ExportData(data_config)
            # print('tau_wrapper advanceintime ImportData')
            self.ImportData(data_config)
            print('tau_wrapper advanceintime Finish')
            print(data_config["time"])
            return data_config["time"]
        else:
            self.time += self.time_step
            return 0.0

    def InitializeSolutionStep(self):
        print('TAUWrapper InitializeSolutionStep')
        super(TAUWrapper, self).InitializeSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.InitializeSolutionStep)
            # self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.BreakSolutionLoop)
        print('TAUWrapper InitializeSolutionStep End')

    def SolveSolutionStep(self):
        print('TAUWrapper SolveSolutionStep')
        super(TAUWrapper, self).SolveSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.SolveSolutionStep)

        if not self.coupling_interface_imported:
            for model_part_name, comm_name in self.settings["solver_wrapper_settings"]["model_parts_recv"].items():
                interface_config = {
                    "comm_name" : comm_name.GetString(),
                    "model_part_name" : model_part_name
                }

                self.ImportCouplingInterface(interface_config)
                self.coupling_interface_imported = True
        print('TAUWrapper SolveSolutionStep End')

    def FinalizeSolutionStep(self):
        print('TAUWrapper FinalizeSolutionStep')
        super(TAUWrapper, self).FinalizeSolutionStep()
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.FinalizeSolutionStep)
        print('TAUWrapper FinalizeSolutionStep End')

    def Finalize(self):
        print('TAUWrapper Finalize')
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.BreakSolutionLoop)
        super(TAUWrapper, self).Finalize() # this also does the disconnect
        print('TAUWrapper Finalize End')

    def ImportCouplingInterface(self, interface_config):
        print('TAUWrapper ImportCouplingInterface')
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ExportMesh, interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super(TAUWrapper, self).ImportCouplingInterface(interface_config)
        print('TAUWrapper ImportCouplingInterface End')

    def ExportCouplingInterface(self, interface_config):
        print('TAUWrapper ExportCouplingInterface')
        if self.controlling_external_solver:
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ImportMesh, interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super(TAUWrapper, self).ExportCouplingInterface(interface_config)
        print('TAUWrapper ExportCouplingInterface End')

    def ImportData(self, data_config):
        print('TAUWrapper ImportData')
        if self.controlling_external_solver and data_config["type"] == "coupling_interface_data":
            # CoSim imports, the external solver exports
            self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ExportData, data_config["interface_data"].name)
        super(TAUWrapper, self).ImportData(data_config)
        print('TAUWrapper ImportData End')

    def ExportData(self, data_config):
        print('TAUWrapper ExportData')
        if self.controlling_external_solver:
            if data_config["type"] == "coupling_interface_data":
                # CoSim exports, the external solver imports
                self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ImportData, data_config["interface_data"].name)
            elif data_config["type"] == "convergence_signal":
                return # we control the ext solver, no need for sending a convergence signal
        super(TAUWrapper, self).ExportData(data_config)
        print('TAUWrapper ExportData End')

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
        print("tau_wrapper __SendControlSignal", data_config)
        self.ExportData(data_config)

    def __CheckExternalSolverProcess(self):
        if hasattr(self, 'external_solver_process') and self.external_solver_process.poll() is None:
            _, process_stderr = self.external_solver_process.communicate()
            if process_stderr:
                raise Exception("{} terminated with the following error:\n{}".format(self._ClassName(), process_stderr.decode('ascii')))
