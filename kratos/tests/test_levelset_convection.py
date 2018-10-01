from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseDistance(x, y, z):
    if (x <= 5.0):
        return -0.16*x**2 + 0.8*x
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
        except FileNotFoundError as e:
            pass

    def test_levelset_convection(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("levelset_convection_process_mesh")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, BaseDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, ConvectionVelocity(node.X,node.Y,node.Z))

        for node in model_part.Nodes:
            if node.X < 0.001:
                node.Fix(KratosMultiphysics.DISTANCE)

        import new_linear_solver_factory
        linear_solver = new_linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "SkylineLUFactorizationSolver"}"""))

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

        self.assertAlmostEqual(max_distance, 0.7904255118014996)
        self.assertAlmostEqual(min_distance,-0.11710292469993094)


if __name__ == '__main__':
    KratosUnittest.main()