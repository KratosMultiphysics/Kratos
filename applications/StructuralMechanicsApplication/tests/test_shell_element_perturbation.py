from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import time as timer
from gid_output_process import GiDOutputProcess

class TestShellElementPerturbation(KratosUnittest.TestCase):

    def execute_thick_shell_test(self, mp, model_part_primal):
        with open("response_function_tests/one_element_test_shell_surface_load.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(mp, parameters)
        primal_analysis.Initialize()

        startTime = timer.time()
        if not primal_analysis.time < primal_analysis.end_time:
            primal_analysis.end_time += 1

        primal_analysis.RunSolutionLoop()

        gid_output = GiDOutputProcess(model_part_primal,
                             "response_function_tests/one_element_test",
                            KratosMultiphysics.Parameters("""
                                {
                                    "result_file_configuration" : {
                                        "gidpost_flags": {
                                            "GiDPostMode": "GiD_PostBinary",
                                            "WriteDeformedMeshFlag": "WriteUndeformed",
                                            "WriteConditionsFlag": "WriteConditions",
                                            "MultiFileFlag": "SingleFile"
                                        },
                                        "file_label"          : "step",
                                        "output_control_type" : "step",
                                         "output_frequency"    : 1,
                                        "nodal_results"       : ["DISPLACEMENT" , "ROTATION"]
                                    }
                                }
                                """))

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()

    def force_perturbation_calculations(self, mp, model_part_primal):
        LHS_1 = KratosMultiphysics.Matrix(6,6)
        RHS_1 = KratosMultiphysics.Vector(12)

        RHS_unperturbed = []
        for condition in model_part_primal.Conditions:
            condition.CalculateLocalSystem(LHS_1,RHS_1,model_part_primal.ProcessInfo)
            RHS_unperturbed.append(RHS_1)

        i = 1
        perturbation_vector_X = []
        perturbed_RHS_X = []
        perturbed_RHS_Y = []
        perturbed_RHS_Z = []
        self.RHS_difference_X = []
        self.RHS_difference_Y = []
        self.RHS_difference_Z = []

        #calculation of the difference in RHS due to displacement perturbation
        for node in model_part_primal.Nodes:
            RHS_difference_X_node = []
            RHS_difference_Y_node = []
            RHS_difference_Z_node = []
            LHS = KratosMultiphysics.Matrix(6,6)
            RHS = KratosMultiphysics.Vector(12)
            displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0)
            
            delta = 0.001
            node.X += delta
            for condition in model_part_primal.Conditions:
                condition.CalculateLocalSystem(LHS,RHS,model_part_primal.ProcessInfo)
                perturbed_RHS_X.append(RHS)
                RHS_difference_X_node.append(RHS - RHS_unperturbed[condition.Id - 1])
            node.X -= delta

            LHS_2 = KratosMultiphysics.Matrix(6,6)
            RHS_2 = KratosMultiphysics.Vector(12)
            node.Y += delta
            for condition in model_part_primal.Conditions:
                condition.CalculateLocalSystem(LHS_2,RHS_2,model_part_primal.ProcessInfo)
                perturbed_RHS_Y.append(RHS_2)
                RHS_difference_Y_node.append(RHS_2 - RHS_unperturbed[condition.Id - 1])
            node.Y -= delta

            LHS_3 = KratosMultiphysics.Matrix(6,6)
            RHS_3 = KratosMultiphysics.Vector(12)
            node.Z += delta
            for condition in model_part_primal.Conditions:
                condition.CalculateLocalSystem(LHS_3,RHS_3,model_part_primal.ProcessInfo)
                perturbed_RHS_Z.append(RHS_3)
                RHS_difference_Z_node.append(RHS_3 - RHS_unperturbed[condition.Id - 1])
            node.Z -= delta

            self.RHS_difference_X.append(RHS_difference_X_node)
            self.RHS_difference_Y.append(RHS_difference_Y_node)
            self.RHS_difference_Z.append(RHS_difference_Z_node)

            i+=1


    # TODO def _check_outputs()
    def _check_outputs_thin_shell_element(self, RHS_difference_X, RHS_difference_Y, RHS_difference_Z):

        for i in range(0,4):
            for j in range(0,12):
                self.assertAlmostEqual(RHS_difference_X[i][j], self.RHS_difference_X[i][0][j], 4)
                self.assertAlmostEqual(RHS_difference_Y[i][j], self.RHS_difference_Y[i][0][j], 4)
                self.assertAlmostEqual(RHS_difference_Z[i][j], self.RHS_difference_Z[i][0][j], 4)

    def test_thin_shell_element(self):
        # these values are not correct, it is just a trial
        RHS_difference_node_1_X = [-7.27596e-12,-24.449,83.3333,0,-12.2245,83.3333,-7.27596e-12,-12.2245,41.6667,0,-24.449,41.6667]
        RHS_difference_node_1_Y = [24.449,0,79.6661,12.2245,0,39.8331,12.2245,0,39.8331,24.449,0,79.6661]
        RHS_difference_node_1_Z = [-83.3333,-79.6661,0,-83.3333,-39.8331,-1.45519e-11,-41.6667,-39.8331,0,-41.6667,-79.6661,0]
        RHS_difference_node_2_X = [-7.27596e-12,24.449,41.6667,0,12.2245,41.6667,0,12.2245,83.3333,0,24.449,83.3333]
        RHS_difference_node_2_Y = [-24.449,0,-79.6661,-12.2245,0,-39.8331,-12.2245,0,-39.8331,-24.449,0,-79.6661]
        RHS_difference_node_2_Z = [-41.6667,79.6661,1.45519e-11,-41.6667,39.8331,0,-83.3333,39.8331,-1.45519e-11,-83.3333,79.6661,1.45519e-11]
        RHS_difference_node_3_X = [7.27596e-12,-12.2245,-83.3333,0,-24.449,-83.3333,-7.27596e-12,-24.449,-41.6667,7.27596e-12,-12.2245,-41.6667]
        RHS_difference_node_3_Y = [12.2245,0,39.8331,24.449,0,79.6661,24.449,0,79.6661,12.2245,0,39.8331]
        RHS_difference_node_3_Z = [83.3333,-39.8331,1.45519e-11,83.3333,-79.6661,0,41.6667,-79.6661,-2.91038e-11,41.6667,-39.8331,1.45519e-11]
        RHS_difference_node_4_X = [7.27596e-12,12.2245,-41.6667,0,24.449,-41.6667,0,24.449,-83.3333,7.27596e-12,12.2245,-83.3333]
        RHS_difference_node_4_Y = [-12.2245,0,-39.8331,-24.449,0,-79.6661,-24.449,0,-79.6661,-12.2245,0,-39.8331]
        RHS_difference_node_4_Z = [41.6667,39.8331,1.45519e-11,41.6667,79.6661,0,83.3333,79.6661,0,83.3333,39.8331,1.45519e-11]

        RHS_difference_X = [RHS_difference_node_1_X, RHS_difference_node_2_X, RHS_difference_node_3_X, RHS_difference_node_4_X]
        RHS_difference_Y = [RHS_difference_node_1_Y, RHS_difference_node_2_Y, RHS_difference_node_3_Y, RHS_difference_node_4_Y]
        RHS_difference_Z = [RHS_difference_node_1_Z, RHS_difference_node_2_Z, RHS_difference_node_3_Z, RHS_difference_node_4_Z]

        model = KratosMultiphysics.Model()
        model_part_primal = model.CreateModelPart("Structure" , 2)

        self.execute_thick_shell_test(model, model_part_primal)
        self.force_perturbation_calculations(model, model_part_primal)
        self._check_outputs_thin_shell_element(RHS_difference_X, RHS_difference_Y, RHS_difference_Z)
        

if __name__ == '__main__':
    KratosUnittest.main()