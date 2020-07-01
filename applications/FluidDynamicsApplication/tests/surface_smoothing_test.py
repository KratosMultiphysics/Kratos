from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
import math
import os
from random import seed
from random import random

from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestSurfaceSmoothing(KratosUnittest.TestCase):

    def test_smoothing_3d_cube(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        KratosMultiphysics.ModelPartIO(GetFilePath("SurfaceSmoothingTest/three_dim_symmetrical_cube")).ReadModelPart(model_part)

        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.1)

        model_part.SetBufferSize(2)

        # seed random number generator
        seed(1)

        for node in model_part.Nodes:
            # generate random numbers between 0-1
            noise = (random()-0.5)*0.001
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,
                (noise+math.sqrt((node.X-0.0)**2+(node.Y-0.0)**2+(node.Z-0.0)**2) - 0.008) )
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(
                kratos_comm, model_part).Execute()

        dimensions = 3
        avg_num_elements = 4
        find_neighbouring_elements_process = KratosMultiphysics.FindElementalNeighboursProcess(
            model_part, dimensions, avg_num_elements).Execute()

        # Set IS_STRUCTURE to define contact with solid
        for node in model_part.Nodes:
            if node.X == 0.0 or node.Y == 0.0 or node.Z == 0.0:
                node.Set(KratosMultiphysics.CONTACT, True)
            else:
                node.Set(KratosMultiphysics.CONTACT, False)

        #verify the orientation of the skin
        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
        flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        KratosMultiphysics.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

        # Calculate boundary normals
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(
            model_part, 3)

        KratosMultiphysics.ComputeNodalGradientProcess(
            model_part,
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        linear_solver = linear_solver_factory.ConstructSolver(KratosMultiphysics.Parameters("""
            {
                "solver_type"         		: "amgcl",
                "max_iteration"       		: 10000,
                "tolerance"           		: 1e-9,
                "provide_coordinates" 		: false,
                "gmres_krylov_space_dimension"   	: 40,
                "smoother_type"       		: "gauss_seidel",
                "krylov_type"         		: "lgmres",
                "use_block_matrices_if_possible" 	: false,
                "coarsening_type"     		: "aggregation",
                "scaling"             		: true
            }
            """)
        )

        smoothing_process = KratosCFD.SurfaceSmoothingProcess(
           model_part,
           linear_solver)

        for _ in range(1): # surface smoothing can be called multiple times
            smoothing_process.Execute()

        node = (model_part.Nodes)[428]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), 0.00022061696483317178)
        node = (model_part.Nodes)[582]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), -0.0004197732005270892)
        node = (model_part.Nodes)[634]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), -0.00034132053885862427)

        # gid_output = GiDOutputProcess(model_part,
        #                            "smoothing_test_3D",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","DISTANCE_GRADIENT","NORMAL"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

    def test_smoothing_3d_cube_with_parameters(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("MainModelPart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        KratosMultiphysics.ModelPartIO(GetFilePath("SurfaceSmoothingTest/three_dim_symmetrical_cube")).ReadModelPart(model_part)

        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.1)

        model_part.SetBufferSize(2)

        # seed random number generator
        seed(1)

        for node in model_part.Nodes:
            # generate random numbers between 0-1
            noise = (random()-0.5)*0.001
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,
                (noise+math.sqrt((node.X-0.0)**2+(node.Y-0.0)**2+(node.Z-0.0)**2) - 0.008) )
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(
                kratos_comm, model_part).Execute()

        dimensions = 3
        avg_num_elements = 4
        find_neighbouring_elements_process = KratosMultiphysics.FindElementalNeighboursProcess(
            model_part, dimensions, avg_num_elements).Execute()

        # Set IS_STRUCTURE to define contact with solid
        for node in model_part.Nodes:
            if node.X == 0.0 or node.Y == 0.0 or node.Z == 0.0:
                node.Set(KratosMultiphysics.CONTACT, True)
            else:
                node.Set(KratosMultiphysics.CONTACT, False)

        #verify the orientation of the skin
        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
        flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        KratosMultiphysics.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

        # Calculate boundary normals
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(
            model_part, 3)

        KratosMultiphysics.ComputeNodalGradientProcess(
            model_part,
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name": "MainModelPart",
                "linear_solver_settings"       : {
                    "solver_type"         		: "amgcl",
                    "max_iteration"       		: 10000,
                    "tolerance"           		: 1e-9,
                    "provide_coordinates" 		: false,
                    "gmres_krylov_space_dimension"   	: 40,
                    "smoother_type"       		: "gauss_seidel",
                    "krylov_type"         		: "lgmres",
                    "use_block_matrices_if_possible" 	: false,
                    "coarsening_type"     		: "aggregation",
                    "scaling"             		: true
                }
            }
            """)

        smoothing_process_with_parameters = KratosCFD.SurfaceSmoothingProcess(
           current_model,
           settings)

        for _ in range(1): # surface smoothing can be called multiple times
            smoothing_process_with_parameters.Execute()

        node = (model_part.Nodes)[428]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), 0.00022061696483317178)
        node = (model_part.Nodes)[582]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), -0.0004197732005270892)
        node = (model_part.Nodes)[634]
        self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), -0.00034132053885862427)

        # gid_output = GiDOutputProcess(model_part,
        #                            "smoothing_test_3D",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","DISTANCE_GRADIENT","NORMAL"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()