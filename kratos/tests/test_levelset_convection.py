﻿import KratosMultiphysics
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


if __name__ == '__main__':
    KratosUnittest.main()