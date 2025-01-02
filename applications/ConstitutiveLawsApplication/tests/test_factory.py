# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
if kratos_utils.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def SelectAndVerifyLinearSolver(settings, skiptest):
    # The mechanical solver selects automatically the fastest linear-solver available
    # this might not be appropriate for a test, therefore in case nothing is specified,
    # the previous default linear-solver is set
    if not settings["solver_settings"].Has("linear_solver_settings"):
        default_lin_solver_settings = KratosMultiphysics.Parameters("""{
                "solver_type": "LinearSolversApplication.sparse_lu"
            }""")
        settings["solver_settings"].AddValue("linear_solver_settings", default_lin_solver_settings)

    solver_type = settings["solver_settings"]["linear_solver_settings"]["solver_type"].GetString()
    solver_type_splitted = solver_type.split(".")
    if len(solver_type_splitted) == 2:
        # this means that we use a solver from an application
        # hence we have to check if it exists, otherwise skip the test
        app_name = solver_type_splitted[0]
        solver_name = solver_type_splitted[1]
        if not kratos_utils.CheckIfApplicationsAvailable(app_name):
            skiptest('Application "{}" is needed for the specified solver "{}" but is not available'.format(app_name, solver_name))


@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            SelectAndVerifyLinearSolver(ProjectParameters, self.skipTest)

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


class SimpleSmallDeformationPlasticityMCTest(TestFactory):
    file_name = "SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_MC_test"

class SimpleSmallDeformationPlasticityVMTest(TestFactory):
    file_name = "SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_VM_test"

class SimpleSmallDeformationPlasticityDPTest(TestFactory):
    file_name = "SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_DP_test"

class SimpleSmallDeformationPlasticityTTest(TestFactory):
    file_name = "SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_T_test"

class BigCubeSmallDeformationPlasticityMCTest(TestFactory):
    file_name = "BigCubeSmallDeformationPlasticity/bigcube_small_deformation_plasticity_MC_test"

class BigCubeSmallDeformationPlasticityVMTest(TestFactory):
    file_name = "BigCubeSmallDeformationPlasticity/bigcube_small_deformation_plasticity_VM_test"

class BigCubeSmallDeformationPlasticityDPTest(TestFactory):
    file_name = "BigCubeSmallDeformationPlasticity/bigcube_small_deformation_plasticity_DP_test"

class BigCubeSmallDeformationPlasticityTTest(TestFactory):
    file_name = "BigCubeSmallDeformationPlasticity/bigcube_small_deformation_plasticity_T_test"

class SerialParallelRuleOfMixturesCubeDamageTest(TestFactory):
    file_name = "SerialParallelRuleOfMixturesCube/serial_parallel_damage_test"

class PlasticDamageTest(TestFactory):
    file_name = "PlasticDamageModel/plastic_damage_test"

class AnisotropyTest(TestFactory):
    file_name = "AnisotropyCube/anisotropy_test"

class Anisotropy2DTest(TestFactory):
    file_name = "AnisotropyCube/anisotropy_2d_test"

class InitialStateInelasticityTest(TestFactory):
    file_name = "InitialStateInelasticity/initial_state2_test"

class InitialStateInelasticity2Test(TestFactory):
    file_name = "InitialStateInelasticity/initial_state3_test"

class SmallDeformationPlasticityTest(TestFactory):
    file_name = "SmallDeformationPlasticity/small_deformation_plasticity_test"

class SimpleJ2PlasticityTest(TestFactory):
    file_name = "SimpleSmallDeformationPlasticity/plasticity_j2_cube_test"

class TensileTestStructuralTest(TestFactory):
    file_name = "TensileTestStructural/TensileTestStructural"

class HighCycleFatigueTest(TestFactory):
    file_name = "HighCycleFatigue/high_cycle_fatigue_test"

class AutomatedInitialDamageTest(TestFactory):
    file_name = "AutomatedInitialDamageProcess/automated_initial_damage_process_test"

class TractionSeparationLawTest(TestFactory):
    file_name = "TractionSeparationLaw/traction_separation_law_test"

class CurveByPointsPlasticityTest(TestFactory):
    file_name = "CurveByPointsPlasticity/plastic_test"



if __name__ == '__main__':
    KratosUnittest.main()
