from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics import SeismicApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class SeismicTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            # The mechanical solver selects automatically the fastest linear-solver available
            # this might not be appropriate for a test, therefore in case nothing is specified,
            # the previous default linear-solver is set
            if not ProjectParameters["solver_settings"].Has("linear_solver_settings"):
                # check if running in MPI because there we use a different default linear solver
                if IsDistributedRun():
                    default_lin_solver_settings = KratosMultiphysics.Parameters("""{
                        "solver_type" : "amesos",
                        "amesos_solver_type" : "Amesos_Klu"
                    }""")

                else:
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
            if ProjectParameters["problem_data"]["echo_level"].GetInt() == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
            else:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.INFO)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = StructuralMechanicsAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Finalize()



class FiberBeamElementTest(SeismicTestFactory):
    file_name = "fiber_beam_element_test/reinforced_concrete_hysteritic_test"



if __name__ == '__main__':
    KratosUnittest.main()
