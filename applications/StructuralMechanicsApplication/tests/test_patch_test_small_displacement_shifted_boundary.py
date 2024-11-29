import KratosMultiphysics

from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestPatchTestSmallDisplacementShiftedBoundary(KratosUnittest.TestCase):

    def setUp(self):
        self.print_output = False

    def testPatchTestSmallDisplacementShiftedBoundary2D3N(self):
        self.work_folder = "patch_test/small_disp_shifted_boundary"
        self.file_name = "patch_test_2D_tension_tri_shifted_boundary"
        self.__RunTest()

    def __RunTest(self):
        class SmallDisplacementShiftedBoundaryPatchTestAnalysis(structural_mechanics_analysis.StructuralMechanicsAnalysis):

            def ModifyInitialGeometry(self):
                super().ModifyInitialGeometry()
                level_set_x_position = 0
                for node in self._GetSolver().GetComputingModelPart().Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, node.X - level_set_x_position)

            def ApplyBoundaryConditions(self):
                super().ApplyBoundaryConditions()

                penalty_factor = 1.0e3
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.PENALTY_COEFFICIENT] = penalty_factor
                for condition in self._GetSolver().GetComputingModelPart().GetSubModelPart("shifted_boundary").Conditions:
                    condition.SetValue(KratosMultiphysics.DISPLACEMENT, [0.0,0.0,0.0])

        with open(f"{self.work_folder}/{self.file_name}_parameters.json",'r') as parameter_file:
            # Read and customize settings
            self.settings = KratosMultiphysics.Parameters(parameter_file.read())
            if self.print_output:
                self.__AddOutput()

            # Creating the test
            model = KratosMultiphysics.Model()
            simulation = SmallDisplacementShiftedBoundaryPatchTestAnalysis(model, self.settings)
            simulation.Run()

    def __AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "Parameters"    : {
                "model_part_name"        : "Structure",
                "output_name"            : "TO_BE_DEFINED",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["DISPLACEMENT","DISTANCE"],
                        "nodal_flags_results"         : ["INTERFACE","BOUNDARY","ACTIVE"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        gid_output_settings["Parameters"]["output_name"].SetString(f"{self.work_folder}/{self.file_name}")
        self.settings["output_processes"]["gid_output"].Append(gid_output_settings)

if __name__ == '__main__':
    # KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()