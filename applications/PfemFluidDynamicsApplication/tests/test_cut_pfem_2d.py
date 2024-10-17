import sys
import time
import importlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis import PfemFluidDynamicsAnalysis


class TestCutPfem(KratosUnittest.TestCase):
    def testDamBreak2D(self):

        # Auxiliary class defining the test levelset
        class DamBreak2DAnalysisStage(PfemFluidDynamicsAnalysis):
            def __init__(self, model, project_parameters):
                super().__init__(model,project_parameters)

            def InitializeSolutionStep(self):
                # This does the remeshing and the "json" BC imposition
                super().InitializeSolutionStep()

                # Set levelset in the new mesh
                epsilon = 1.0e-5
                for node in self._GetSolver().main_model_part.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, node.Y-epsilon)
                    if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.001:
                        if node.Is(KratosMultiphysics.FREE_SURFACE):
                            node.Set(KratosMultiphysics.FREE_SURFACE, False)
                    # Fix the negative nodes to avoid them to move when doing the convection
                    if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:

                        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                        node.Fix(KratosMultiphysics.DISPLACEMENT_Y)

                # Deactivate negative elements -> To be moved to a PFEM-CutFEM solver
                for element in self._GetSolver().main_model_part.Elements:
                    n_pos = 0
                    n_neg = 0
                    for node in element.GetGeometry():
                        if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0:
                            n_pos += 1
                        else:
                            n_neg += 1
                    if n_neg == element.GetGeometry().PointsNumber():
                        element.Set(KratosMultiphysics.ACTIVE, False)
                    else:
                        element.Set(KratosMultiphysics.ACTIVE, True)

        # Test specifications
        self.case_name ="dam_break_2d"
        self.work_folder = f"cut_pfem_tests/{self.case_name}"
        self._ReadAndCustomizeTestSettings()

        # Test run
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            simulation = DamBreak2DAnalysisStage(self.model, self.parameters)
            simulation.Run()

        # Check some results
        model_part = self.model.GetModelPart("PfemFluidModelPart")
        n_nodes = model_part.NumberOfNodes()
        n_elems = model_part.NumberOfElements()
        self.assertEqual(n_nodes, 73)
        self.assertEqual(n_elems, 121)

        node_10 = model_part.GetNode(10)
        self.assertVectorAlmostEqual(node_10, [0.24,0.0,0.0], self.check_places)
        self.assertAlmostEqual(node_10.GetSolutionStepValue(KratosMultiphysics.DISTANCE), -1.0e-5, self.check_places)
        self.assertAlmostEqual(node_10.GetSolutionStepValue(KratosMultiphysics.PRESSURE), -2713.2575470526017, self.check_places)
        self.assertVectorAlmostEqual(node_10.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.002965895169715622,-0.0023034592967681178,0], self.check_places)

        node_20 = model_part.GetNode(20)
        self.assertVectorAlmostEqual(node_20, [0.25073523965843825,0.1525073863658415,0.0], self.check_places)
        self.assertAlmostEqual(node_20.GetSolutionStepValue(KratosMultiphysics.DISTANCE), 0.1525614286, self.check_places)
        self.assertAlmostEqual(node_20.GetSolutionStepValue(KratosMultiphysics.PRESSURE), -699.1205536766438, self.check_places)
        self.assertVectorAlmostEqual(node_20.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.005757319479420786,-0.02134741138617036,0], self.check_places)

    def setUp(self):
        self.print_output = False
        self.check_places = 10
        self.check_absolute_tolerance = 1.0e-8
        self.check_relative_tolerance = 1.0e-10

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('dam_break_2d.time')

    def _ReadAndCustomizeTestSettings(self):
        # Read the simulation settings
        settings_filename = "ProjectParameters.json"
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(settings_filename,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # If required, add the output process to the test settings
        if self.print_output:
            self._AddOutput()

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "Parameters"    : {
                "model_part_name"        : "PfemFluidModelPart",
                "output_name"            : "TO_BE_DEFINED",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "MultipleFiles"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["VELOCITY","PRESSURE","DENSITY"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        gid_output_settings["Parameters"]["output_name"].SetString(self.case_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

if __name__ == '__main__':
    KratosUnittest.main()
