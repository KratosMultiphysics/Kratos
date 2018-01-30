from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import math

def GetFilePath(fileName):
    import os
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestRedistance(KratosUnittest.TestCase):

    def _ExpectedDistance(self,x,y,z):
        d = x
        if( d > 0.2):
            d = 0.2
        if( d < -0.2):
            d = -0.2
        return x 
        #return -(math.sqrt(x**2+y**2+z**2) - 0.4)

    def test_model_part_sub_model_parts(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        KratosMultiphysics.ModelPartIO(GetFilePath("coarse_sphere")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)


        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0, self._ExpectedDistance(node.X,node.Y,node.Z)  )

        import new_linear_solver_factory
        linear_solver = new_linear_solver_factory.ConstructSolver( KratosMultiphysics.Parameters( """ { "solver_type" : "SkylineLUFactorizationSolver" } """ ) )
                
        model_part.CloneTimeStep(1.0)

        max_iterations = 2
        KratosMultiphysics.VariationalDistanceCalculationProcess3D(model_part, linear_solver, max_iterations).Execute()

        max_distance = -1.0;
        min_distance = +1.0
        for node in model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        self.assertAlmostEqual(max_distance, 0.4473365725705744)
        self.assertAlmostEqual(min_distance,-0.5049982310145651)
        
        
if __name__ == '__main__':
    KratosUnittest.main()
