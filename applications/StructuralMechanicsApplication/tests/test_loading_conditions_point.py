from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestLoadingConditionsPoint(KratosUnittest.TestCase):
    def _execute_point_condition_test(self, Dimension, LoadingType):
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)

        if LoadingType == "load":
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
            mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        elif LoadingType == "moment":
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
        else:
            raise RuntimeError("Wrong LoadingType")
        
        # create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)
        
        # ensure that the property 1 is created
        mp.GetProperties()[1]

        if LoadingType == "load":
            for node in mp.Nodes:
                node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
                node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
                node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)

            if Dimension == 2:
                cond = mp.CreateNewCondition("PointLoadCondition2D1N", 1, [1], mp.GetProperties()[1])
            elif Dimension == 3:
                cond = mp.CreateNewCondition("PointLoadCondition3D1N", 1, [1], mp.GetProperties()[1])
            else:
                raise RuntimeError("Wrong Dimension")

        elif LoadingType == "moment":
            for node in mp.Nodes:
                node.AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X)
                node.AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y)
                node.AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z)

            if Dimension == 2:
                cond = mp.CreateNewCondition("PointMomentCondition2D1N", 1, [1], mp.GetProperties()[1])
            elif Dimension == 3:
                cond = mp.CreateNewCondition("PointMomentCondition3D1N", 1, [1], mp.GetProperties()[1])
            else:
                raise RuntimeError("Wrong Dimension")
        else:
            raise RuntimeError("Wrong LoadingType")
        

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)
        
        # first we apply a load to the condition 
        load_on_cond = (1.8, 2.6,-11.47)

        if LoadingType == "load":
            cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)
        elif LoadingType == "moment":
            cond.SetValue(StructuralMechanicsApplication.POINT_MOMENT, load_on_cond)
        
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertEqual(rhs[0], load_on_cond[0])
        self.assertEqual(rhs[1], load_on_cond[1])
        if Dimension == 3:
            self.assertEqual(rhs[2], load_on_cond[2])

        # now we apply a load to the node
        nodal_load = (5.5, 1.2,-9,3)

        if LoadingType == "load":
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD, nodal_load)
        elif LoadingType == "moment":
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_MOMENT, nodal_load)
        
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertEqual(rhs[0], load_on_cond[0] + nodal_load[0])
        self.assertEqual(rhs[1], load_on_cond[1] + nodal_load[1])
        if Dimension == 3:
            self.assertEqual(rhs[2], load_on_cond[2] + nodal_load[2])
            

    def test_PointLoadCondition2D1N(self):
        self._execute_point_condition_test(Dimension=2, LoadingType="load")

    def test_PointLoadCondition3D1N(self):
        self._execute_point_condition_test(Dimension=2, LoadingType="load")

    def test_PointMomentCondition2D1N(self):
        self._execute_point_condition_test(Dimension=2, LoadingType="moment")

    def test_PointMomentCondition3D1N(self):
        self._execute_point_condition_test(Dimension=3, LoadingType="moment")
    

if __name__ == '__main__':
    KratosUnittest.main()
