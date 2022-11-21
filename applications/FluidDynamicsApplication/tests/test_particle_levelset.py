import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import math

from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseDistance(x, y, z):
    d = math.sqrt( (x-0.6)**2 + (y-0.6)**2 ) - 0.15
    return d

def ConvectionVelocity(x, y, z):
    vel = KratosMultiphysics.Vector(3, 0.0)
    vel[0] = -(y-0.5)
    vel[1] = (x-0.5)
    vel[2] = 0.0
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
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VOLUME)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("Cavity/square10")).ReadModelPart(model_part)
#        KratosMultiphysics.ModelPartIO(GetFilePath("../../../kratos/tests/auxiliar_files_for_python_unittest/mdpa_files/levelset_convection_process_mesh")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        particles_mp = current_model.CreateModelPart("Particles")
        particles_mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        particles_mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        particles_mp.SetBufferSize(2)

        for node in model_part.Nodes:
            node.GetValue(KratosMultiphysics.NODAL_VOLUME)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, BaseDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, ConvectionVelocity(node.X,node.Y,node.Z))

        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))


        gid_output = GiDOutputProcess(model_part,
                                "levelset_test_2D_algebraic_new",
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "nodal_results"       : ["DISTANCE","VELOCITY"]
                                        }
                                    }
                                    """)
                                )


        # gid_output = GiDOutputProcess(particles_mp,
        #                         "pls",
        #                         KratosMultiphysics.Parameters("""
        #                             {
        #                                 "result_file_configuration" : {
        #                                     "gidpost_flags": {
        #                                         "GiDPostMode": "GiD_PostBinary",
        #                                         "WriteDeformedMeshFlag": "WriteDeformed",
        #                                         "WriteConditionsFlag": "WriteConditions",
        #                                         "MultiFileFlag": "MultipleFiles"
        #                                     },
        #                                     "node_output": true,
        #                                     "nodal_results"       : ["DISTANCE"]
        #                                 }
        #                             }
        #                             """)
        #                         )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()

        levelset_convection_settings = KratosMultiphysics.Parameters("""{
                "max_CFL" : 1.0,
                "max_substeps" : 0,
                "eulerian_error_compensation" : false,
                "element_type" : "levelset_convection_supg"
            }""")
        convection_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                model_part,
                linear_solver,
                levelset_convection_settings)

        print(85)
        pls_utility = KratosMultiphysics.FluidDynamicsApplication.ParticleLevelsetUtility2D(model_part, particles_mp, KratosMultiphysics.DISTANCE, 0.3)
        locator = KratosMultiphysics.BinBasedFastPointLocator2D(model_part)
        locator.UpdateSearchDatabase()
        print(89)

        particle_convection_utility = KratosMultiphysics.ParticleConvectUtily2D(locator)

        redistance_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "Main",
            "nodal_area_variable": "NODAL_VOLUME",
            "max_levels" : 20,
            "max_distance" : 2.0,
            "calculate_exact_distances_to_plane": false
        }""")

        distance_calculator = KratosMultiphysics.ParallelDistanceCalculationProcess2D(current_model, redistance_settings)

        nsteps = 200
        dt = 2*math.pi/nsteps
        for i in range(1,nsteps):
            t = dt*i
            model_part.CloneTimeStep(t)
            particles_mp.CloneTimeStep(t)
            print("time=",model_part.ProcessInfo[KratosMultiphysics.TIME])
            print("delta_time=",model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME])

            for node in model_part.Nodes:
                if(node.X < 0.0001 or node.X > 0.999 or node.Y < 0.0001 or node.Y > 0.9999):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,[0,0,0])
                    node.Fix(KratosMultiphysics.DISTANCE)

            pls_utility.Prepare()

            convection_process.Execute()

            #move particles
            particle_convection_utility.MoveParticles_Substepping(particles_mp,1)
            #particle_convection_utility.MoveParticles_RK4(particles_mp)
            pls_utility.MoveAndCorrect(locator)
            #err

            distance_calculator.Execute()
            pls_utility.MoveAndCorrect(locator)

            max_distance = -1.0
            min_distance = +1.0
            for node in model_part.Nodes:
                d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                max_distance = max(max_distance, d)
                min_distance = min(min_distance, d)

            # self.assertAlmostEqual(max_distance, 0.733304104543163)
            # self.assertAlmostEqual(min_distance,-0.06371359024393097)

            gid_output.ExecuteInitializeSolutionStep()
            gid_output.PrintOutput()
            gid_output.ExecuteFinalizeSolutionStep()
            #err
        gid_output.ExecuteFinalize()




if __name__ == '__main__':
    KratosUnittest.main()