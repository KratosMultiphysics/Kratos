# CoSimulation imports
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.StructuralMechanicsApplication # needed for some variables

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities
import subprocess, os

def Create(settings, model, solver_name):
    return RemoteControlledSolverWrapper(settings, model, solver_name)

class RemoteControlledSolverWrapper(CoSimulationSolverWrapper):
    """Interface for the CFD-Solver TAU
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        wrapper_settings = self.settings["solver_wrapper_settings"]


        model_part_utilities.CreateModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

        # start_external_solver = wrapper_settings["start_external_solver"].GetBool()
        # if start_external_solver:
        #     command_txt = wrapper_settings["external_solver_start_command"].GetString()
        #     path_to_tau = wrapper_settings["path_to_tau"].GetString()
        #     # QUESTION: shall we use this:
        #     parent_path = os.path.join(os.path.dirname(__file__), '../..')
        #     # OR this:
        #     # parent_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
        #     path_to_tau_solver = parent_path + '/helpers/tau_solver.py'
        #     path_to_tau_solver = parent_path + '/helpers/TauSolver.py'
        #     tau_input_file = wrapper_settings["tau_input_file"].GetString()
        #     tau_log_file = wrapper_settings["tau_log_file"].GetString()
        #     # command_args = wrapper_settings["external_solver_arguments"].GetStringArray()
        #     command_args = [path_to_tau, path_to_tau_solver, tau_input_file, tau_log_file]
        #     cs_tools.cs_print_info(self._ClassName(), 'Running external solver with command: "{}" | arguments: "{}"'.format(command_txt, command_args))

        #     full_command = [command_txt]
        #     full_command.extend(command_args)
        #     # self.external_solver_process = subprocess.Popen(full_command, stderr=subprocess.PIPE, start_new_session=True) # TODO check what to use here
        #     # self.external_solver_process = subprocess.Popen(full_command)

        # self.controlling_external_solver = wrapper_settings["controlling_external_solver"].GetBool()

    def Initialize(self):
        super().Initialize()

        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = { "model_part_name" : model_part_name }
            self.ImportCouplingInterface(interface_config)

    def AdvanceInTime(self, current_time):
        # self.__CheckExternalSolverProcess() # TODO check why this is blocking
        self.__SendControlSignal("AdvanceInTime")
        return 0.0
        # data_config = {
        #     "type"       : "time",
        #     "time"     : current_time
        # }
        # self.ExportData(data_config)
        # self.ImportData(data_config)
        # return 100.0 #data_config["time"]

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.__SendControlSignal("InitializeSolutionStep")

    def SolveSolutionStep(self):
        super().SolveSolutionStep()
        self.__SendControlSignal("InitializeSolutionStep")

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.__SendControlSignal("InitializeSolutionStep")

    def Finalize(self):
        self.__SendControlSignal("exit")
        super().Finalize() # this also does the disconnect

    def ImportCouplingInterface(self, interface_config):
        self.__SendControlSignal("ExportMesh", interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super().ImportCouplingInterface(interface_config)

    def ExportCouplingInterface(self, interface_config):
        self.__SendControlSignal("ImportMesh", interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super().ExportCouplingInterface(interface_config)

    def ImportData(self, data_config):
        # CoSim imports, the external solver exports
        self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ExportData, data_config["interface_data"].name)
        super().ImportData(data_config)

    # def ExportData(self, data_config):
    #     if self.controlling_external_solver:
    #         if data_config["type"] == "coupling_interface_data":
    #             # CoSim exports, the external solver imports
    #             self.__SendControlSignal(KratosCoSim.CoSimIO.ControlSignal.ImportData, data_config["interface_data"].name)
    #         elif data_config["type"] == "convergence_signal":
    #             return # we control the ext solver, no need for sending a convergence signal
    #     super().ExportData(data_config)

    def _GetIOType(self):
        return "kratos_co_sim_io"

    def __SendControlSignal(self, signal, identifier=""):
        data_config = {
            "type"           : "control_signal",
            "control_signal" : signal,
            "identifier"     : identifier
        }
        self.ExportData(data_config)

    def __CheckExternalSolverProcess(self):
        if hasattr(self, 'external_solver_process') and self.external_solver_process.poll() is None:
            _, process_stderr = self.external_solver_process.communicate()
            if process_stderr:
                raise Exception("{} terminated with the following error:\n{}".format(self._ClassName(), process_stderr.decode('ascii')))
