from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestAdjointLoadingConditions(KratosUnittest.TestCase):

    def __CheckSensitivityMatrix(self, sen_matrix, reference_sen_matrix, digits_to_check=5):
        if ((sen_matrix.Size1() != len(reference_sen_matrix)) or (sen_matrix.Size2() != len(reference_sen_matrix[0]))):
            raise Exception('Matrix sizes does not fit!')
        for i in range(sen_matrix.Size1()):
            for j in range(sen_matrix.Size2()):
                self.assertAlmostEqual(sen_matrix[i,j], reference_sen_matrix[i][j], digits_to_check)

    def _SurfaceLoadCondition3D4N(self, prefix = ""):
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
        props = mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition("AdjointSemiAnalytic" + prefix + "SurfaceLoadCondition3D4N", 1, [1,2,3,4], props)

        # first we apply a constant SURFACE_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] =  0.0
        load_on_cond[1] =  0.0
        load_on_cond[2] =  1.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        sen_matrix = KratosMultiphysics.Matrix(0,0)

        # settings for finite differencing
        mp.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True
        mp.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = 1e-5

        # check w.r.t. to design variable SURFACE_LOAD
        reference_res_1 = [[0.25,0,0,0.25,0,0,0.25,0,0,0.25,0,0],[0,0.25,0,0,0.25,0,0,0.25,0,0,0.25,0],[0,0,0.25,0,0,0.25,0,0,0.25,0,0,0.25]]
        cond.CalculateSensitivityMatrix(StructuralMechanicsApplication.SURFACE_LOAD, sen_matrix, mp.ProcessInfo)
        self.__CheckSensitivityMatrix(sen_matrix, reference_res_1)

        # check w.r.t. to design variable SHAPE_SENSITIVITY
        reference_res_2 = ((0,0,-0.166667,0,0,-0.166667,0,0,-0.0833333,0,0,-0.0833333),(0,0,-0.166667,0,0,-0.0833333,0,0,-0.0833333,0,0,-0.166667),
        (0,0,1.25e-06,0,0,8.33328e-07,0,0,4.16669e-07,0,0,8.33336e-07),(0,0,0.166667,0,0,0.166667,0,0,0.0833333,0,0,0.0833333),
        (0,0,-0.0833333,0,0,-0.166667,0,0,-0.166667,0,0,-0.0833333),(0,0,8.33336e-07,0,0,1.25e-06,0,0,8.33336e-07,0,0,4.16664e-07),
        (0,0,0.0833333,0,0,0.0833333,0,0,0.166667,0,0,0.166667),(0,0,0.0833333,0,0,0.166667,0,0,0.166667,0,0,0.0833333),(0,0,4.16664e-07,0,0,8.33333e-07,0,0,1.25e-06,0,0,8.33336e-07),
        (0,0,-0.0833333,0,0,-0.0833333,0,0,-0.166667,0,0,-0.166667),(0,0,0.166667,0,0,0.0833333,0,0,0.0833333,0,0,0.166667),(0,0,8.33336e-07,0,0,4.16667e-07,0,0,8.33336e-07,0,0,1.25e-06))
        cond.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, sen_matrix, mp.ProcessInfo)
        self.__CheckSensitivityMatrix(sen_matrix, reference_res_2)

        # change loading
        load_on_cond[2] =  0.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        cond.SetValue(KratosMultiphysics.PRESSURE, 1.0)

        # check w.r.t. to design variable PRESSURE
        reference_res_3 = [0,0,0.25,0,0,0.25,0,0,0.25,0,0,0.25]
        cond.CalculateSensitivityMatrix(KratosMultiphysics.PRESSURE, sen_matrix, mp.ProcessInfo)
        if ((sen_matrix.Size2() != len(reference_res_3)) or sen_matrix.Size1() != 1):
            raise Exception('Matrix sizes does not fit!')
        for i in range(sen_matrix.Size2()):
            self.assertAlmostEqual(sen_matrix[0,i], reference_res_3[i], 5)


    def _LineLoadCondition3D2N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,1.0,0.0)

        #ensure that the property 1 is created
        props = mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition("AdjointSemiAnalytic" + prefix + "LineLoadCondition3D2N", 1, [1,2], props)
        cond.SetValue(KratosMultiphysics.LOCAL_AXIS_2, [-1.0, 1.0, 0.0])

        #first we apply a constant LINE_LOAD to theh condition
        Line_Load_i = 10000.00/math.sqrt(2) #apply a 45 degrees load

        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 0.0
        load_on_cond[1] = -Line_Load_i
        load_on_cond[2] = -Line_Load_i
        cond.SetValue(StructuralMechanicsApplication.LINE_LOAD,load_on_cond)
        sen_matrix = KratosMultiphysics.Matrix(0,0)

        # settings for finite differencing
        mp.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = False
        mp.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = 1e-6

        # check w.r.t. to design variable SURFACE_LOAD
        reference_res_1 = ((0.707107,0,0,0.707107,0,0),(0,0.707107,0,0,0.707107,0),(0,0,0.707107,0,0,0.707107))
        cond.CalculateSensitivityMatrix(StructuralMechanicsApplication.LINE_LOAD, sen_matrix, mp.ProcessInfo)
        self.__CheckSensitivityMatrix(sen_matrix, reference_res_1)

        # check w.r.t. to design variable SHAPE_SENSITIVITY
        reference_res_2 = ((0,2499.999376,2499.999376,0,2499.999376,2499.999376),(0,2499.999376,2499.999376,0,2499.999376,2499.999376),
        (0,-0.001250,-0.001250,0,-0.0012450,-0.001250),(0,-2500.000624,-2500.000624,0,-2500.000624,-2500.000624),
        (0,-2500.000624,-2500.000624,0,-2500.000624,-2500.000624),(0,-0.001250,-0.001250,0,-0.001250,-0.001250))
        cond.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, sen_matrix, mp.ProcessInfo)
        self.__CheckSensitivityMatrix(sen_matrix, reference_res_2)

    def test_SDSurfaceLoadCondition3D4N(self):
        self._SurfaceLoadCondition3D4N(prefix = "SmallDisplacement")

    def test_AdjointSurfaceLoadCondition3D4N(self):
        self._SurfaceLoadCondition3D4N()

    def test_SDLineLoadCondition3D2N(self):
        self._LineLoadCondition3D2N(prefix = "SmallDisplacement")

    def test_AdjointLineLoadCondition3D2N(self):
        self._LineLoadCondition3D2N()

if __name__ == '__main__':
    KratosUnittest.main()
