import os

# This test................

# Importing the Kratos Library
import KratosMultiphysics as Kratos
import numpy as np
# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis
# This utility will control the execution scope

debug_mode=False
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class GentleInjectionTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with open(self.file_parameters,'r') as parameter_file:
            parameters_harsh = Kratos.Parameters(parameter_file.read())

        fluid_solver_settings = parameters_harsh['fluid_parameters']['solver_settings']
        self.max_nonlinear_iterations = fluid_solver_settings['maximum_iterations'].GetInt()

        if not debug_mode:
            parameters_harsh['fluid_parameters']['output_processes'] = Kratos.Parameters('''{}''')

        parameters_gentle = Kratos.Parameters(parameters_harsh.WriteJsonString())

        self.gentle_initiation_interval = 0.1
        parameters_gentle['coupling']['gentle_coupling_initiation']['initiation_interval'].SetDouble(self.gentle_initiation_interval)

        # Create Model
        model_harsh = Kratos.Model()
        model_gentle = Kratos.Model()

        self.test_harsh_injection = GentleInjectionAnalysis(model_harsh, parameters_harsh)
        self.test_gentle_injection = GentleInjectionAnalysis(model_gentle, parameters_gentle)
        self.test_harsh_injection.name = 'harsh'
        self.test_gentle_injection.name = 'gentle'

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test_harsh_injection.Run()
            self.test_gentle_injection.Run()
            times = np.array(self.test_harsh_injection.times)
            n_iterations_harsh = np.array(self.test_harsh_injection.n_iterations)
            n_iterations_gentle = np.array(self.test_gentle_injection.n_iterations)

            # check that the number of iterations needed with gentle injection is always less
            assert((n_iterations_harsh >= n_iterations_gentle).all())
            assert((n_iterations_harsh > n_iterations_gentle).any())

            if debug_mode:
                import matplotlib.pyplot as plt
                plt.plot(times, n_iterations_harsh, label='initiation_interval=' + str(0.0))
                plt.plot(times, n_iterations_gentle, label='initiation_interval=' + str(round(self.gentle_initiation_interval, 1)))
                plt.xlabel('time (s)')
                plt.ylabel('nonlinear iterations')
                ax = plt.gca()
                ax.set_ylim([0, self.max_nonlinear_iterations])
                plt.legend()
                plt.savefig('nonlinear_iterations.pdf')

class GentleInjectionAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, parameters=Kratos.Parameters("{}")):
        super().__init__(model, parameters)
        self.n_iterations = []
        self.times = []
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'gentle_injection_tests/')
        self.problem_name = parameters

    def FinalizeSolutionStep(self):
        time = self.fluid_model_part.ProcessInfo[Kratos.TIME]
        n_iterations = self.fluid_model_part.ProcessInfo[Kratos.NL_ITERATION_NUMBER]
        self.times.append(time)
        self.n_iterations.append(n_iterations)
        super(GentleInjectionAnalysis, self).FinalizeSolutionStep()