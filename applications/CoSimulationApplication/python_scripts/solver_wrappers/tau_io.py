from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# Other imports
import os, subprocess, time

def Create(settings, model, solver_name):
    print("Creating TAU-IO")
    return TAUIO(settings, model, solver_name)

class TAUIO(CoSimulationIO):
    """IO for the legacy EMPIRE_API
    """
    def __init__(self, settings, model, solver_name):
        super(TAUIO, self).__init__(settings, model, solver_name)

        tau_path = "/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq/tau.py"
        # QUESTION: shall we use this:
        parent_path = os.path.join(os.path.dirname(__file__), '..')
        # OR this:
        # parent_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
        tau_solver_path = parent_path + '/helpers/tau_solver.py'
        tau_input_file = "airfoil_Structured.cntl"
        tau_log_file = "log_TAU.out"
        self.tau_subprocess = subprocess.Popen(
            ["python", tau_path, tau_solver_path, tau_input_file, tau_log_file], stdin=subprocess.PIPE)

    def ImportCouplingInterface(self, interface_config):
        print('ImportCouplingInterface in TAUIO not implemented yet!!!')

    def ExportCouplingInterface(self, interface_config):
        print('ExportCouplingInterface in TAUIO not implemented yet!!!')

    def ImportData(self, data_config):
        print('ImportData in TAUIO not implemented yet!!!')

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "control_signal":
            control_signal_key = data_config["signal"] + '\n'
            self.tau_subprocess.stdin.write(control_signal_key.encode())
            self.tau_subprocess.stdin.flush()
            time.sleep(10)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def Finalize(self):
        print("TAU IO Finalize")
        self.tau_subprocess.kill()

    def PrintInfo(self):
        print("This is the TAU-IO")
