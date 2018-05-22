import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication as kfd
have_required_applications = True
missing_applications_message = ["Missing required application(s): ",]
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    have_required_applications = False
    missing_applications_message.append("ExternalSolversApplication")

try:
    import KratosMultiphysics.AdjointFluidApplication as kaf
    have_adjoint_fluid = True
except ImportError:
    have_required_applications = False
    missing_applications_message.append("AdjointFluidApplication")

try:
    import KratosMultiphysics.HDF5Application as kh5
except ImportError:
    have_required_applications = False
    missing_applications_message.append("HDF5Application")

from adjoint_fluid_analysis import AdjointFluidAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_required_applications," ".join(missing_application_message))
class AdjointFluidTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def testCylinder(self):
        work_folder = "CylinderTest"
        primal_settings_file_name = "cylinder_fluid_parameters.json"
        adjoint_settings_file_name = "cylinder_adjoint_parameters.json"

        with WorkFolderScope(work_folder):
            self._run_test(primal_settings_file_name,adjoint_settings_file_name)

            kratos_utilities.DeleteFileIfExisting("cylinder_2d.time")

    def _run_test(self,primal_settings_file_name,adjoint_settings_file_name):
        model = km.Model()
        settings = km.Parameters(r'''{}''')

        with open(primal_parameter_file_name,'r') as primal_parameter_file:
            settings.AddValue("primal settings", km.Parameters(primal_parameter_file.read()))

        with open(adjoint_parameter_file_name,'r') as adjoint_parameter_file:
            settings.AddValue("adjoint settings", km.Parameters(adjoint_parameter_file.read()))


        # to check the results: add output settings block if needed
        if self.print_output:
            settings["adjoint settings"].AddValue("output_configuration", km.Parameters(r'''{
                "result_file_configuration" : {
                    "gidpost_flags"       : {
                        "GiDPostMode"           : "GiD_PostBinary",
                        "WriteDeformedMeshFlag" : "WriteDeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "file_label"          : "time",
                    "output_control_type" : "step",
                    "output_frequency"    : 1,
                    "body_output"         : true,
                    "node_output"         : false,
                    "skin_output"         : false,
                    "plane_output"        : [],
                    "nodal_results"       : ["VELOCITY","PRESSURE"],
                    "gauss_point_results" : []
                },
                "point_data_configuration"  : []
            }'''))

        primal_analysis = FluidDynamicsAnalysis(model,settings)
        primal_analysis.Run()
        adjoint_analysis = AdjointFluidAnalysis(model,settings)
        adjoint_analysis.Run()

if __name__ == '__main__':
    test_case = AdjointFluidTest()
    test_case.setUp()
    test_case.testAdjointFluidAnalysis()
    test_case.tearDown()

