import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

import math


class EmbeddedSolverDirichletCircleAnalysis(ConvectionDiffusionAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)

    def ModifyInitialGeometry(self):
        # Set the level set function (circle of radius 1.0)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            dist = 1.0 - math.sqrt((node.X)**2 + (node.Y)**2 + (node.Z)**2)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, dist)

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        # Set source term
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX, 0, 1.0)

        # Set penalty coefficient
        self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.PENALTY_COEFFICIENT] = 1.0e0
        # Set boundary value for Nitsche imposition of DBC
        for elem in self._GetSolver().GetComputingModelPart().Elements:
            elem.SetValue(KratosMultiphysics.ConvectionDiffusionApplication.EMBEDDED_SCALAR, 0.0)


@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestEmbeddedSolver(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "test_convection_diffusion_embedded_solver"
        self.print_output = False
        self.print_reference_values = False

    def _readAndCustomizeTestSettings(self):
        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(self.file_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # If required, add the output process to the test settings
        if self.print_output:
            self._addOutput()

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self._addReferenceValuesOutput()
        else:
            self._addReferenceValuesCheck()

    def _readAndCustomizeTestSettingsDistanceModification(self):
        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(self.file_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
        self.parameters["solver_settings"]["use_distance_modification"].SetBool(True)

        # If required, add the output process to the test settings
        if self.print_output:
            self._addOutput()
            self.parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString("square_distance_modification")

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self._addReferenceValuesOutput()
        else:
            self._addReferenceValuesCheck()

    def _readAndCustomizeTestSettingsMLSConstraints(self):
        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(self.file_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
        self.parameters["solver_settings"]["use_mls_constraints"].SetBool(True)

        # If required, add the output process to the test settings
        if self.print_output:
            self._addOutput()
            self.parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString("square_mls_constraints")

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self._addReferenceValuesOutput()
        else:
            self._addReferenceValuesCheck()

    def _addOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "ThermalModelPart",
                "output_name"            : "square",
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
                        "output_interval"             : 1,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["DISTANCE","TEMPERATURE"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _addReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["TEMPERATURE"],
                "output_file_name" : "reference_embedded_dirichlet_circle",
                "model_part_name"  : "ThermalModelPart",
                "time_frequency"   : 1.0
            }
        }""")
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _addReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["TEMPERATURE"],
                "input_file_name"      : "reference_embedded_dirichlet_circle",
                "model_part_name"      : "ThermalModelPart",
                "tolerance"            : 1.0e-6,
                "relative_tolerance"   : 1.0e-8,
                "time_frequency"       : 1.0
            }
        }""")
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

    def tearDown(self):
        if not self.print_output:
            with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
                KratosUtilities.DeleteFileIfExisting("square.post.bin")
                KratosUtilities.DeleteFileIfExisting("square_distance_modification.post.bin")
                KratosUtilities.DeleteFileIfExisting("square_mls_constraints.post.bin")
                KratosUtilities.DeleteFileIfExisting("test_convection_diffusion_embedded_solver.post.lst")
                KratosUtilities.DeleteFileIfExisting("test_convection_diffusion_embedded_solver_distance_modification.post.lst")
                KratosUtilities.DeleteFileIfExisting("test_convection_diffusion_embedded_solver_mls_constraints.post.lst")

    def testEmbeddedSolverDirichletCircle(self):
        self.file_name = "ProjectParameters.json"

        self._readAndCustomizeTestSettings()

        # test solver without distance modification and without MLS constraints
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            embedded_simulation = EmbeddedSolverDirichletCircleAnalysis(self.model, self.parameters)
            embedded_simulation.Run()

    def testEmbeddedSolverDirichletCircleDistanceModification(self):
        self.file_name = "ProjectParameters.json"

        self._readAndCustomizeTestSettingsDistanceModification()

        # test solver without distance modification and without MLS constraints
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            embedded_simulation = EmbeddedSolverDirichletCircleAnalysis(self.model, self.parameters)
            embedded_simulation.Run()

    def testEmbeddedSolverDirichletCircleMLSConstraints(self):
        self.file_name = "ProjectParameters.json"

        self._readAndCustomizeTestSettingsMLSConstraints()

        # test solver without distance modification and without MLS constraints
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            embedded_simulation = EmbeddedSolverDirichletCircleAnalysis(self.model, self.parameters)
            embedded_simulation.Run()


if __name__ == "__main__":
    KratosUnittest.main()
