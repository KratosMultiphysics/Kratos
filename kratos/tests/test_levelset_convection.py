import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

# from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseDistance(x, y, z):
    if (x <= 5.0):
        return -0.16*x**2 + 0.8*x
    else:
        return 0.0

def BaseJumpedDistance(x, y, z):
    if (x >= 5.0 and x <= 15.0):
        return 1.0
    else:
        return 0.0

def ConvectionVelocity(x, y, z):
    vel = KratosMultiphysics.Vector(3, 0.0)
    vel[0] = 1.0
    return vel

class TestLevelSetConvection(KratosUnittest.TestCase):

    def tearDown(self):
        # Remove the .time file
        try:
            os.remove('levelset_convection_process_mesh.time')
        except :
            pass

    def test_levelset_convection(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/levelset_convection_process_mesh")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, BaseDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, ConvectionVelocity(node.X,node.Y,node.Z))

        for node in model_part.Nodes:
            if node.X < 0.001:
                node.Fix(KratosMultiphysics.DISTANCE)

        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))

        model_part.CloneTimeStep(40.0)

        KratosMultiphysics.LevelSetConvectionProcess2D(
            KratosMultiphysics.DISTANCE,
            model_part,
            linear_solver).Execute()

        max_distance = -1.0
        min_distance = +1.0
        for node in model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        self.assertAlmostEqual(max_distance, 0.733304104543163)
        self.assertAlmostEqual(min_distance,-0.06371359024393097)

    def test_levelset_convection_BFECC(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/levelset_convection_process_mesh")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, BaseJumpedDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, ConvectionVelocity(node.X,node.Y,node.Z))

        for node in model_part.Nodes:
            if node.X < 0.001:
                node.Fix(KratosMultiphysics.DISTANCE)

        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))

        model_part.CloneTimeStep(30.0)

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(
                kratos_comm, model_part).Execute()

        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
            model_part,
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "levelset_variable_name" : "DISTANCE",
            "levelset_convection_variable_name" : "VELOCITY",
            "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
            "max_CFL" : 1.0,
            "max_substeps" : 0,
            "eulerian_error_compensation" : true,
            "cross_wind_stabilization_factor" : 0.7,
            "algebraic_stabilization" : false
        }""")
        KratosMultiphysics.LevelSetConvectionProcess2D(
            model_part,
            linear_solver,
            levelset_convection_settings).Execute()

        max_distance = -1.0
        min_distance = +1.0
        for node in model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        # gid_output = GiDOutputProcess(model_part,
        #                            "levelset_test_2D_supg",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","VELOCITY"]
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

        self.assertAlmostEqual(max_distance, 1.0617777301844604)
        self.assertAlmostEqual(min_distance, -0.061745786561321375)

    def test_levelset_convection_BFECC_algebraic(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/levelset_convection_process_mesh")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, BaseJumpedDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, ConvectionVelocity(node.X,node.Y,node.Z))

        for node in model_part.Nodes:
            if node.X < 0.001:
                node.Fix(KratosMultiphysics.DISTANCE)

        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))

        model_part.CloneTimeStep(1.0)

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(
                kratos_comm, model_part).Execute()

        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
            model_part,
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "levelset_variable_name" : "DISTANCE",
            "levelset_convection_variable_name" : "VELOCITY",
            "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
            "max_CFL" : 1.0,
            "max_substeps" : 0,
            "eulerian_error_compensation" : false,
            "cross_wind_stabilization_factor" : 0.7,
            "algebraic_stabilization" : true,
            "high_order_terms" : true
        }""")
        KratosMultiphysics.LevelSetConvectionProcess2D(
            model_part,
            linear_solver,
            levelset_convection_settings).Execute()

        max_distance = -1.0
        min_distance = +1.0
        for node in model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        # gid_output = GiDOutputProcess(model_part,
        #                            "levelset_test_2D_algebraic",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","VELOCITY"]
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

        self.assertAlmostEqual(max_distance, 1.0004068724542385)
        self.assertAlmostEqual(min_distance, -0.0004103252870070794)


if __name__ == '__main__':
    KratosUnittest.main()