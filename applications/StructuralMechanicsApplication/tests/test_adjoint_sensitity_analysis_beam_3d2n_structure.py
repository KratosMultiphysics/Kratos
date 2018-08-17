from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

try:
    from KratosMultiphysics.HDF5Application import *
    has_hdf5_application = True
except ImportError:
    has_hdf5_application = False

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

def solve_primal_problem():
    with open("beam_test_parameters.json",'r') as parameter_file:
        ProjectParametersPrimal = Parameters( parameter_file.read())

    # To avoid many prints
    if (ProjectParametersPrimal["problem_data"]["echo_level"].GetInt() == 0):
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    model_primal = Model()
    primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
    primal_analysis.Run()

def _get_test_working_dir():
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(this_file_dir, "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n")

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisBeamStructure(KratosUnittest.TestCase):

    # called only once for this class, opposed of setUp()
    @classmethod
    def setUpClass(cls):
        with controlledExecutionScope(_get_test_working_dir()):
            solve_primal_problem()

    def test_local_stress_response(self):
        #Create the adjoint solver
        with controlledExecutionScope(_get_test_working_dir()):
            with open("beam_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = Parameters( parameter_file.read())

            model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
            model_adjoint = Model()

            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
            adjoint_analysis.Run()

            # Check sensitivities for the parameter I22
            reference_values = [-87622.77099397512, 38125.18144970003, 625.0029074349038]
            sensitivities_to_check = []
            element_list = [1,6,10]
            for element_id in element_list:
                sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(I22_SENSITIVITY))

            self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
            self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
            self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_nodal_displacement_response(self):
        # Create the adjoint solver
        with controlledExecutionScope(_get_test_working_dir()):
            with open("beam_test_nodal_disp_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = Parameters( parameter_file.read())

            model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
            model_adjoint = Model()

            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)

            adjoint_analysis.Run()

            # Check sensitivities for the parameter I22
            reference_values = [-454.1027959305903, -378.2187594016309, -6.200311358415619]
            sensitivities_to_check = []
            element_list = [1,6,10]
            for element_id in element_list:
                sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(I22_SENSITIVITY))

            self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
            self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
            self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_strain_energy_response(self):
        # Create the adjoint solver
        with controlledExecutionScope(_get_test_working_dir()):
            with open("beam_test_strain_energy_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = Parameters( parameter_file.read())

            model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
            model_adjoint = Model()

            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)

            adjoint_analysis.Run()

            # Check sensitivities for the parameter I22
            reference_values = [-9082.055913430777, -7564.375188032661, -124.00616062621255]
            sensitivities_to_check = []
            element_list = [1,6,10]
            for element_id in element_list:
                sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(I22_SENSITIVITY))

            self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
            self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
            self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        with controlledExecutionScope(_get_test_working_dir()):
            kratos_utilities.DeleteFileIfExisting("Beam_structure.time")
            for file_name in os.listdir():
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

if __name__ == '__main__':
    KratosUnittest.main()
