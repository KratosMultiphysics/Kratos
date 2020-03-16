from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math


class TestAdjointLoadingConditionsSurface(KratosUnittest.TestCase):

    def __CheckSensitivityMatrix(self, sen_matrix, reference_res_list,digits_to_check=5):
        for i in range(sen_matrix.Size1()):
            for j in range(sen_matrix.Size2()):
                self.assertAlmostEqual(sen_matrix[i,j],reference_res_list[i][j],digits_to_check)

    def _SimplestSurfaceLoadCondition3D4N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,1.0,1.0,0.0)
        mp.CreateNewNode(4,0.0,1.0,0.0)

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition(prefix + "AdjointSemiAnalyticSurfaceLoadCondition3D4N", 1, [1,2,3,4], mp.GetProperties()[1])

        #first we apply a constant SURFACE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] =  0.0
        load_on_cond[1] =  0.0
        load_on_cond[2] =  1.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        sen_matrix = KratosMultiphysics.Matrix(0,0)
        #lhs = KratosMultiphysics.Matrix(0,0)
        #rhs = KratosMultiphysics.Vector(0)
        #cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        #print(rhs)

        mp.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True
        mp.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = 1e-5

        reference_res_1 = [[0.25,0,0,0.25,0,0,0.25,0,0,0.25,0,0],[0,0.25,0,0,0.25,0,0,0.25,0,0,0.25,0],[0,0,0.25,0,0,0.25,0,0,0.25,0,0,0.25]]
        cond.CalculateSensitivityMatrix(StructuralMechanicsApplication.SURFACE_LOAD, sen_matrix, mp.ProcessInfo)
        self.__CheckSensitivityMatrix(sen_matrix, reference_res_1)

        reference_res_2 = ((0,0,-0.166667,0,0,-0.166667,0,0,-0.0833333,0,0,-0.0833333),(0,0,-0.166667,0,0,-0.0833333,0,0,-0.0833333,0,0,-0.166667),
        (0,0,1.25e-06,0,0,8.33328e-07,0,0,4.16669e-07,0,0,8.33336e-07),(0,0,0.166667,0,0,0.166667,0,0,0.0833333,0,0,0.0833333),
        (0,0,-0.0833333,0,0,-0.166667,0,0,-0.166667,0,0,-0.0833333),(0,0,8.33336e-07,0,0,1.25e-06,0,0,8.33336e-07,0,0,4.16664e-07),
        (0,0,0.0833333,0,0,0.0833333,0,0,0.166667,0,0,0.166667),(0,0,0.0833333,0,0,0.166667,0,0,0.166667,0,0,0.0833333),(0,0,4.16664e-07,0,0,8.33333e-07,0,0,1.25e-06,0,0,8.33336e-07),
        (0,0,-0.0833333,0,0,-0.0833333,0,0,-0.166667,0,0,-0.166667),(0,0,0.166667,0,0,0.0833333,0,0,0.0833333,0,0,0.166667),(0,0,8.33336e-07,0,0,4.16667e-07,0,0,8.33336e-07,0,0,1.25e-06))
        cond.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, sen_matrix, mp.ProcessInfo)
        self.__CheckSensitivityMatrix(sen_matrix, reference_res_2)


    #def test_SDSimplestSurfaceLoadCondition3D4N(self):
    #    self._SimplestSurfaceLoadCondition3D4N("SmallDisplacement")

    def test_SimplesAdjointtSurfaceLoadCondition3D4N(self):
        self._SimplestSurfaceLoadCondition3D4N()



if __name__ == '__main__':
    KratosUnittest.main()
