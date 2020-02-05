from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

has_hdf5_application = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")

def solve_primal_problem(file_name):
    with open(file_name,'r') as parameter_file:
        ProjectParametersPrimal = KratosMultiphysics.Parameters( parameter_file.read())

    # To avoid many prints
    if (ProjectParametersPrimal["problem_data"]["echo_level"].GetInt() == 0):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    else:
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.INFO)

    model_primal = KratosMultiphysics.Model()
    primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
    primal_analysis.Run()

def _get_test_working_dir():
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(this_file_dir, "adjoint_sensitivity_analysis_tests/adjoint_spring_damper_element_3d2n/")

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisSpringDamperStructure(KratosUnittest.TestCase):

    # called only once for this class, opposed of setUp()
    @classmethod
    def setUpClass(cls):
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            solve_primal_problem("ProjectParameters.json")

    def test_displacement_response(self):
        #Create the adjoint solver
        with KratosUnittest.WorkFolderScope(_get_test_working_dir() , __file__):
            with open("AdjointParameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = KratosMultiphysics.Parameters( parameter_file.read())

            model_adjoint = KratosMultiphysics.Model()
            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
            adjoint_analysis.Run()

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        with KratosUnittest.WorkFolderScope(_get_test_working_dir(), __file__):
            kratos_utilities.DeleteFileIfExisting("Spring_structure.time")
            for file_name in os.listdir():
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

if __name__ == '__main__':
    KratosUnittest.main()
