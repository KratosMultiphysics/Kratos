from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.IgaApplication.iga_analysis import IgaAnalysis

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Other imports
import os

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class IgaTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            # The mechanical solver selects automatically the fastest linear-solver available
            # this might not be appropriate for a test, therefore in case nothing is specified,
            # the previous default linear-solver is set
            if not ProjectParameters["solver_settings"].Has("linear_solver_settings"):
                default_lin_solver_settings = KratosMultiphysics.Parameters("""{
                "solver_type": "ExternalSolversApplication.super_lu",
                    "max_iteration": 500,
                    "tolerance": 1e-9,
                    "scaling": false,
                    "symmetric_scaling": true,
                    "verbosity": 0
                }""")
                ProjectParameters["solver_settings"].AddValue("linear_solver_settings", default_lin_solver_settings)

            self.modify_parameters(ProjectParameters)

            # To avoid many prints
            if (ProjectParameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = IgaAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Finalize()

class BeamLinear5pTest(IgaTestFactory):
    file_name = "beam_linear_5p/beam_linear_5p"

    
if __name__ == '__main__':
    KratosUnittest.main()