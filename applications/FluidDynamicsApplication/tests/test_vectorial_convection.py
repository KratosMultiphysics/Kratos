import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseDistance(x, y, z):
        return -0.16*x+ 0.8*x

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


    def test_levelset_convection_BDF(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosCFD.FRACTIONAL_VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        KratosMultiphysics.ModelPartIO(GetFilePath("levelset_convection_process_mesh")).ReadModelPart(model_part)
        model_part.SetBufferSize(3)


        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_Y, model_part)
        # KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_Z, model_part)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.1)
        # model_part.ProcessInfo.SetValue(KratosMultiphysics.CROSS_WIND_STABILIZATION_FACTOR, 0.0)

        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        self.time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        # model_part.CloneTimeStep(0.0)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0, 1.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 1, 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 2, 0.75)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 1, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 2, 0.0)
            node.SetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY, [0.0, 0.0, 0.0])



        # for node in model_part.Nodes:
        #     if node.X < 0.001 or node.X>49.99:
        #         node.Fix(KratosMultiphysics.DISTANCE)

        # for node in model_part.Nodes:
        #     if node.X < 0.001:
        #         node.Fix(KratosMultiphysics.DISTANCE)


        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))


        gid_output = GiDOutputProcess(model_part,
                                      "testing_uxue/test_ux",
                                      KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                    "file_label"          : "time",
                                                    "output_control_type" : "time",
                                                    "output_interval"     : 0.5,
                                                "nodal_results"       : ["FRACTIONAL_VELOCITY","VELOCITY"],
                                                "nodal_nonhistorical_results":[]
                                            }
                                        }
                                        """)
                                      )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()



        while self.time <20:
            self.time += 0.01
            model_part.CloneTimeStep(self.time)
            model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.time)


            KratosMultiphysics.TimeDiscretization.BDF(2).ComputeAndSaveBDFCoefficients(model_part.ProcessInfo)
            # bdf_vec = [1.0/0.01,-1.0/0.01,0.0]
            # model_part.ProcessInfo.SetValue(KratosMultiphysics.BDF_COEFFICIENTS, bdf_vec)
            KratosMultiphysics.FindGlobalNodalNeighboursProcess(model_part).Execute()

            levelset_convection_settings = KratosMultiphysics.Parameters("""{

            "element_type": "levelset_convection_bdf"
            }""")
            KratosCFD.VectorialConvectionProcess2D(
                model_part,
                linear_solver,
                levelset_convection_settings).Execute()
            print("hola")
            gid_output.ExecuteInitializeSolutionStep()
            gid_output.PrintOutput()
            gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()




if __name__ == '__main__':
    KratosUnittest.main()