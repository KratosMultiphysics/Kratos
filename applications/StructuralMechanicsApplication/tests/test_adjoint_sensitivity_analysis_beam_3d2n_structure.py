from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

has_hdf5_application = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")

def solve_primal_problem():
    with open("beam_test_parameters.json",'r') as parameter_file:
        ProjectParametersPrimal = KratosMultiphysics.Parameters( parameter_file.read())

    # To avoid many prints
    if (ProjectParametersPrimal["problem_data"]["echo_level"].GetInt() == 0):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    model_primal = KratosMultiphysics.Model()
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
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            solve_primal_problem()

    def test_local_stress_response(self):
        #Create the adjoint solver
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            with open("beam_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = KratosMultiphysics.Parameters( parameter_file.read())

            model_part_name = ProjectParametersAdjoint["solver_settings"]["model_part_name"].GetString()
            model_adjoint = KratosMultiphysics.Model()

            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
            adjoint_analysis.Run()

            # Check sensitivities for the parameter I22
            reference_values = [-87.62277093392399, 9.497391494932984, 38.125186783868, 0.6250049974719261, 0.15624887499699122]
            sensitivities_to_check = []
            element_list = [1,2,3,4,5,6,10]
            for element_id in element_list:
                sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 3)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[4], reference_values[1], 3)
        self.assertAlmostEqual(sensitivities_to_check[5], reference_values[2], 3)
        self.assertAlmostEqual(sensitivities_to_check[6], reference_values[3], 3)
        self.assertAlmostEqual(sensitivities_to_check[7], reference_values[4], 5)

    def test_nodal_displacement_response(self):
        # Create the adjoint solver
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            with open("beam_test_nodal_disp_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = KratosMultiphysics.Parameters( parameter_file.read())

            model_part_name = ProjectParametersAdjoint["solver_settings"]["model_part_name"].GetString()

            model_adjoint = KratosMultiphysics.Model()

            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)

            adjoint_analysis.Run()

            # Check sensitivities for the parameter I22
            reference_values = [-0.45410279537614157, -0.37821875982596204, -0.006200296058668847, 0.0004340210813670321]
            sensitivities_to_check = []
            element_list = [1,6,10]
            for element_id in element_list:
                sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 4)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 4)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 4)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[3], 5)

    def test_strain_energy_response(self):
        # Create the adjoint solver
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            with open("beam_test_strain_energy_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = KratosMultiphysics.Parameters( parameter_file.read())

            model_part_name = ProjectParametersAdjoint["solver_settings"]["model_part_name"].GetString()

            model_adjoint = KratosMultiphysics.Model()
            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
            adjoint_analysis.Run()

            # Check sensitivities for the parameter I22
            reference_values = [-9.082055907522943, -7.5643751965193164, -0.12400592117339182, 0.017360843254681547]
            sensitivities_to_check = []
            element_list = [1,6,10]
            for element_id in element_list:
                sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 4)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 4)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 4)
        self.assertAlmostEqual(sensitivities_to_check[3], reference_values[3], 5)


    def test_nodal_reaction_response(self):
        # Create the adjoint solver
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            with open("beam_test_nodal_reaction_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = KratosMultiphysics.Parameters( parameter_file.read())

            model_part_name = ProjectParametersAdjoint["solver_settings"]["model_part_name"].GetString()
            model_adjoint = KratosMultiphysics.Model()
            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
            adjoint_analysis.Run()

        # Check sensitivity for the parameter POINT_LOAD
        reference_values = [-0.31249774999397384, -1.249]
        sensitivities_to_check = []
        sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Conditions[1].GetValue(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)[2])
        sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[10].GetValue(StructuralMechanicsApplication.I22_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 4)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 2)

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            kratos_utilities.DeleteFileIfExisting("Beam_structure.time")
            for file_name in os.listdir():
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

if __name__ == '__main__':
    KratosUnittest.main()
