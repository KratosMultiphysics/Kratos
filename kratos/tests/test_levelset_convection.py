from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import math
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseDistance(x, y, z):
    return 16.0*x*(1.0-x)*y*(1.0 - y)

def ConvectionVelocity(x, y, z):
    norm = (x**2 + y**2)**0.5
    vel = KratosMultiphysics.Vector(3, 0.0)
    if (x + y <= 0.5):
        vel[0] = norm
        vel[1] = norm
    else:
        vel[0] = 2**0.5 - norm
        vel[1] = 2**0.5 - norm
    return vel

class TestLevelSetConvection(KratosUnittest.TestCase):

    def tearDown(self):
        # Remove the .time file
        try:
            os.remove('test_processes.time')
        except FileNotFoundError as e:
            pass

    def test_levelset_convection(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("test_processes")).ReadModelPart(model_part)
        model_part.SetBufferSize(1)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, BaseDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, ConvectionVelocity(node.X,node.Y,node.Z))

        import new_linear_solver_factory
        linear_solver = new_linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "SkylineLUFactorizationSolver"}"""))

        model_part.CloneTimeStep(1.0)

        import gid_output
        gid_io = gid_output.GiDOutput("before_convection")
        gid_io._write_mesh(0.0, model_part)
        gid_io._initialize_results(0.0, model_part)
        gid_io._write_nodal_results(0.0, model_part, KratosMultiphysics.DISTANCE)
        gid_io._write_nodal_results(0.0, model_part, KratosMultiphysics.VELOCITY)
        gid_io._finalize_results()

        max_iterations = 2
        KratosMultiphysics.LevelSetConvectionProcess2D(
            KratosMultiphysics.DISTANCE,
            model_part, 
            linear_solver).Execute()

        import gid_output
        gid_io = gid_output.GiDOutput("after_convection")
        gid_io._write_mesh(0.0, model_part)
        gid_io._initialize_results(0.0, model_part)
        gid_io._write_nodal_results(0.0, model_part, KratosMultiphysics.DISTANCE)
        gid_io._write_nodal_results(0.0, model_part, KratosMultiphysics.VELOCITY)
        gid_io._finalize_results()

        # max_distance = -1.0;
        # min_distance = +1.0
        # for node in model_part.Nodes:
        #     d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
        #     max_distance = max(max_distance, d)
        #     min_distance = min(min_distance, d)

        # self.assertAlmostEqual(max_distance, 0.44556526310761013)
        # self.assertAlmostEqual(min_distance,-0.504972246827639)
        
        
if __name__ == '__main__':
    KratosUnittest.main()
