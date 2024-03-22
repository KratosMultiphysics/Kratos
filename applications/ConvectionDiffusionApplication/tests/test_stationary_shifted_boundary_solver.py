import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

import math

class ShiftedBoundarySolverAnnulusAnalysis(ConvectionDiffusionAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)

    def ModifyInitialGeometry(self):
        # Set the levelset function
        rad = 1.0
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            dist = math.sqrt(node.X**2 + node.Y**2) - rad
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, dist)

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        # Get variables
        conv_diff_settings = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS]
        unknown_variable = conv_diff_settings.GetUnknownVariable()
        source_variable = conv_diff_settings.GetVolumeSourceVariable()

        # Set BCs in active nodes
        for node in self.model.GetModelPart("ThermalModelPart.ImposedTemperature2D_OuterWall").Nodes:
            node.Fix(unknown_variable)
            temp = self.__EvaluateSolutionField(node.X, node.Y)
            node.SetSolutionStepValue(unknown_variable, 0, temp)

        # Set Dirichlet BCs in the immersed boundaries
        sbm_model_part = self.model.GetModelPart("ThermalModelPart.shifted_boundary")
        for condition in sbm_model_part.Conditions:
            coords = condition.GetValue(KratosMultiphysics.INTEGRATION_COORDINATES)
            temp = self.__EvaluateSolutionField(coords[0], coords[1])
            condition.SetValue(unknown_variable, temp)

        # Set source term
        rad = 1
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            if math.sqrt(node.X**2+node.Y**2) > rad:
                heat_flux = self.__EvaluateSourceTerm(node.X, node.Y)
                node.SetSolutionStepValue(source_variable, 0, heat_flux)

        # Set penalty coefficient
        sbm_model_part.ProcessInfo[KratosMultiphysics.PENALTY_COEFFICIENT] = 1e1

    def __EvaluateSolutionField(self, x, y):
        return 0.25*(9.0 - x**2 - y**2 - 2.0*math.log(3) + math.log(x**2+y**2)) + 0.25*math.sin(x)*math.sinh(y)

    def __EvaluateSourceTerm(self, x, y):
        return -0.25*((-2.0 + (2*(x**2+y**2)-4*y**2)/(x**2+y**2)**2) - math.sin(x)*math.sinh(y) + (-2.0 + (2*(x**2+y**2)-4*x**2)/(x**2+y**2)**2) + math.sin(x)*math.sinh(y))

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestShiftedBoundarySolver(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "test_stationary_shifted_boundary_solver"
        self.print_output = False
        self.print_reference_values = False

    def __ReadAndCustomizeTestSettings(self):
        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(self.file_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # If required, add the output process to the test settings
        if self.print_output:
            self.__AddOutput()

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self.__AddReferenceValuesOutput()
        else:
            self.__AddReferenceValuesCheck()

    def __AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "ThermalModelPart",
                "output_name"            : "",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "time",
                        "output_control_type"         : "time",
                        "output_interval"             : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["TEMPERATURE","DISTANCE","HEAT_FLUX"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : [],
                        "nodal_flags_results"         : ["BOUNDARY","INTERFACE"]
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        gid_output_settings["Parameters"]["output_name"].SetString(self.problem_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def __AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["TEMPERATURE"],
                "output_file_name" : "TO_BE_SET",
                "model_part_name"  : "ThermalModelPart",
                "time_frequency"   : 1.0
            }
        }""")
        json_output_settings["Parameters"]["output_file_name"].SetString(self.problem_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def __AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["TEMPERATURE"],
                "input_file_name"      : "TO_BE_SET",
                "model_part_name"      : "ThermalModelPart",
                "tolerance"            : 1.0e-6,
                "relative_tolerance"   : 1.0e-8,
                "time_frequency"       : 1.0
            }
        }""")
        json_check_settings["Parameters"]["input_file_name"].SetString(self.problem_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

    def tearDown(self):
        if not self.print_output:
            KratosUtilities.DeleteFilesEndingWith(self.work_folder, ".bin")
            KratosUtilities.DeleteFilesEndingWith(self.work_folder, ".post.lst")

    def testMovingLeastSquaresExtension(self):
        self.file_name = "ProjectParameters.json"
        self.problem_name = "annulus_3_unstr_ref_0_MLS"
        self.__ReadAndCustomizeTestSettings()
        self.parameters["problem_data"]["problem_name"].SetString(self.problem_name)
        self.parameters["solver_settings"]["extension_operator_type"].SetString("MLS")

        # Test solver with MLS extension
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            sbm_simulation = ShiftedBoundarySolverAnnulusAnalysis(self.model, self.parameters)
            sbm_simulation.Run()

    def testGradientBasedExtension(self):
        self.file_name = "ProjectParameters.json"
        self.problem_name = "annulus_3_unstr_ref_0_gradient_based"
        self.__ReadAndCustomizeTestSettings()
        self.parameters["problem_data"]["problem_name"].SetString(self.problem_name)
        self.parameters["solver_settings"]["extension_operator_type"].SetString("gradient_based")

        # Test solver with gradient-based extension
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            sbm_simulation = ShiftedBoundarySolverAnnulusAnalysis(self.model, self.parameters)
            sbm_simulation.Run()

    def testRadialBasisFunctionExtension(self):
        self.file_name = "ProjectParameters.json"
        self.problem_name = "annulus_3_unstr_ref_0_RBF"
        self.__ReadAndCustomizeTestSettings()
        self.parameters["problem_data"]["problem_name"].SetString(self.problem_name)
        self.parameters["solver_settings"]["extension_operator_type"].SetString("RBF")

        # Test solver with gradient-based extension
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            sbm_simulation = ShiftedBoundarySolverAnnulusAnalysis(self.model, self.parameters)
            sbm_simulation.Run()

if __name__ == "__main__":
    KratosUnittest.main()
