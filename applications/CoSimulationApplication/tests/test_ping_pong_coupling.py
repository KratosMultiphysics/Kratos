import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

import os
class TestPingPong(KratosUnittest.TestCase):
    '''TODO add description
    '''
    def _createTest(self, problem_dir_name, parameter_file_name):
        self.problem_dir_name = problem_dir_name

        full_parameter_file_name = os.path.join(problem_dir_name, parameter_file_name + '.json')

        with open(full_parameter_file_name, 'r') as parameter_file:
            self.cosim_parameters = KM.Parameters(parameter_file.read())

        # To avoid many prints
        echo_level = self.cosim_parameters["problem_data"]["echo_level"].GetInt()
        if (echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
        else:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

    def test_ping_pong(self):
        self._createTest("ping_pong_test", "cosim_ping_pong_parameters")
        CoSimulationAnalysis(self.cosim_parameters).Run()

        DeleteFileIfExisting('ping.log')
        DeleteFileIfExisting('pong.log')
        DeleteFileIfExisting('EMPIRE_datafield_pong_send_data.dat')
        DeleteFileIfExisting('EMPIRE_datafield_ping_send_data.dat')
        DeleteFileIfExisting('EMPIRE_datafield_pong_recv_data.dat')
        DeleteFileIfExisting('EMPIRE_datafield_ping_recv_data.dat')



if __name__ == '__main__':
    KratosUnittest.main()
