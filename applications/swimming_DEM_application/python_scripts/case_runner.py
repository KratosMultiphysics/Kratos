import sys
import os
from DEM_procedures import KratosPrint as Say
import KratosSwimmingDEM as script


class CaseRunner:
    def __init__(self,
                 main_path,
                 algorithm,
                 total_number_of_simulations=None):
        self.main_path = main_path
        self.algorithm = algorithm
        self.n_simulations = total_number_of_simulations

    def FinishAndSayMessage(self, message_start):
        if self.n_simulations:
            message_start += ' out of ' + str(self.n_simulations) + '...\n\n'
        else:
            message_start += '...\n\n'
        Say(message_start)

    def RunCase(self, parameters, simulation_id):
        message_start = '\n\nRunning simulation number ' + str(simulation_id)
        self.FinishAndSayMessage(message_start)

        try:
            with script.Solution(self.algorithm, parameters) as test:
                test.Run()
            error = None
            message_start = '\n\nSuccessfully finished running simulation number ' + str(simulation_id)
            self.FinishAndSayMessage(message_start)
        except:
            error = sys.exc_info()
            message_start = '\n\nFinished running simulation number ' + str(simulation_id)
            self.FinishAndSayMessage(message_start)
            Say('The simulation crashed.\n\n')

        os.chdir(self.main_path)

        return error