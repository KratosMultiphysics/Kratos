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
        self.print_line_length = 85
        self.number_of_trailing_dots = 10

    def SayStartMessage(self):
        message_start = '\n' + '=' * self.print_line_length
        message_start += '\n\nRunning simulation number ' + str(self.simulation_id)
        if self.n_simulations:
            message_start += ' out of ' + str(self.n_simulations) + '...\n'
        else:
            message_start += '...\n'

        if self.identification_text:
            message_start += self.identification_text

        message_start += '_' * self.print_line_length + '\n'

        Say(message_start)

    def FinishAndSayEndMessage(self, message_start):
        if self.n_simulations:
            message_start += ' out of ' + str(self.n_simulations) + '...\n'
        else:
            message_start += '...\n'

        message_start += '=' * self.print_line_length + '\n'
        message_start += ('\n' + ' ' * int(self.print_line_length / 2) + '*') * self.number_of_trailing_dots
        Say(message_start)

    def RunCase(self, parameters, simulation_id, identification_text=''):
        self.simulation_id = simulation_id
        self.identification_text = identification_text

        self.SayStartMessage()

        try:
            with script.Solution(self.algorithm, parameters) as test:
                test.Run()
            error = None
            message_start = 'Successfully finished running simulation number ' + str(simulation_id)
            self.FinishAndSayEndMessage(message_start)
        except Exception:
            error = sys.exc_info()
            message_start = 'Finished running simulation number ' + str(simulation_id)
            self.FinishAndSayEndMessage(message_start)
            Say('The simulation crashed.')

        os.chdir(self.main_path)

        return error