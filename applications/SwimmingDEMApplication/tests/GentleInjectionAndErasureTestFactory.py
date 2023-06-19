import os

# This factory encodes to different tests consisting in a cubic domain of fluid
# with a few particles submerged in it.

# The first test involves a DEM inlet and it checks that the 'gentle' injection of
# particles improves the convergence of the fluid solver. In this case the bounding
# box (which is responsible for deleating particles) is inactive.

# The second test only involves particles initially in the domain (no DEM inlet),
# while the bounding box is active and set to eliminate particles that cross a certain
# height within the fluid domain. Here the improvement of the fluid convergence due to
# gentle elimination of particles is checked.

# In both cases the problem is run twice: one time without the 'gentle' technique and
# the other with it. It is checked that the number of nonlinear iterations is, on average,
# less for the latter case.

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

class GentleInjectionAndErasureTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with open(self.file_parameters_harsh, 'r') as parameter_file:
            self.parameters_harsh = Kratos.Parameters(parameter_file.read())

        with open(self.file_parameters_gentle,'r') as parameter_file:
            self.parameters_gentle = Kratos.Parameters(parameter_file.read())

        fluid_solver_settings = self.parameters_harsh['fluid_parameters']['solver_settings']
        self.max_nonlinear_iterations = fluid_solver_settings['maximum_iterations'].GetInt()

        if not debug_mode:
            for parameters in {self.parameters_harsh, self.parameters_gentle}:
                parameters['do_print_results_option'].SetBool(False)
                parameters['fluid_parameters']['output_processes'] = Kratos.Parameters('''{}''')

        # Create Model
        model_harsh = Kratos.Model()
        model_gentle = Kratos.Model()

        self.test_harsh_injection = GentleAnalysis(model_harsh, self.parameters_harsh)
        self.test_gentle_injection = GentleAnalysis(model_gentle, self.parameters_gentle)
        self.test_harsh_injection.name = 'harsh'
        self.test_gentle_injection.name = 'gentle'

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test_harsh_injection.Run()
            self.test_gentle_injection.Run()
            times = np.array(self.test_harsh_injection.times)
            n_iterations_harsh = np.array(self.test_harsh_injection.n_iterations)
            n_iterations_gentle = np.array(self.test_gentle_injection.n_iterations)

            # Check that the number of iterations needed with gentle injection is on average less
            assert(sum(n_iterations_harsh) > sum(n_iterations_gentle))

            if debug_mode:
                import matplotlib.pyplot as plt
                parameter_value, parameter_name = self.GetGentleParameterValueAndName(self.parameters_gentle)
                plt.plot(times, n_iterations_harsh, label='harsh (' + parameter_name + '=' + str(0.0) + ')')
                plt.plot(times, n_iterations_gentle, label='gentle (' +parameter_name + '=' + str(round(parameter_value, 2)) + ')')
                plt.xlabel('time (s)')
                plt.ylabel('nonlinear iterations')
                ax = plt.gca()
                ax.set_ylim([0, self.max_nonlinear_iterations])
                plt.legend()
                plt.savefig('nonlinear_iterations' + type(self).__name__ + '.pdf')
                plt.close()

class GentleAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, parameters=Kratos.Parameters("{}")):
        super().__init__(model, parameters)
        self.n_iterations = []
        self.times = []
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'fluid_convergence_tests/')
        self.problem_name = parameters

    def FinalizeSolutionStep(self):
        time = self.fluid_model_part.ProcessInfo[Kratos.TIME]
        n_iterations = self.fluid_model_part.ProcessInfo[Kratos.NL_ITERATION_NUMBER]
        self.times.append(time)
        self.n_iterations.append(n_iterations)
        super(GentleAnalysis, self).FinalizeSolutionStep()