import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestAdjointLoadingConditions(KratosUnittest.TestCase):

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
        reference_res_1 = KratosMultiphysics.Matrix(3, 12)
        reference_res_1.fill(0.0)
        reference_res_1[0,0] = reference_res_1[0,3] = reference_res_1[0,6] = reference_res_1[0,9]  = 0.25
        reference_res_1[1,1] = reference_res_1[1,4] = reference_res_1[1,7] = reference_res_1[1,10] = 0.25
        reference_res_1[2,2] = reference_res_1[2,5] = reference_res_1[2,8] = reference_res_1[2,11] = 0.25
        cond.CalculateSensitivityMatrix(StructuralMechanicsApplication.SURFACE_LOAD, sen_matrix, mp.ProcessInfo)
        self.assertMatrixAlmostEqual(reference_res_1, sen_matrix)

        # check w.r.t. to design variable SHAPE_SENSITIVITY
        reference_res_2 = KratosMultiphysics.Matrix(12, 12)
        reference_res_2.fill(0.0)
        reference_res_2[0,2] = reference_res_2[0,5] = reference_res_2[1,2] = reference_res_2[1,11] = -0.166667
        reference_res_2[0,8] = reference_res_2[0,11] = reference_res_2[1,5] = reference_res_2[1,8] = -0.0833333
        reference_res_2[3,2] = reference_res_2[3,5] = reference_res_2[10,2] = reference_res_2[10,11]  = 0.166667
        reference_res_2[3,8] = reference_res_2[3,11] = reference_res_2[10,5] = reference_res_2[10,8] = 0.0833333
        reference_res_2[6,8] = reference_res_2[6,11] = reference_res_2[7,5] = reference_res_2[7,8]  = 0.166667
        reference_res_2[6,2] = reference_res_2[6,5] = reference_res_2[7,2] = reference_res_2[7,11] = 0.0833333
        reference_res_2[4,5] = reference_res_2[4,8] = reference_res_2[9,8] = reference_res_2[9,11] = -0.166667
        reference_res_2[4,2] = reference_res_2[4,11] = reference_res_2[9,2] = reference_res_2[9,5] = -0.0833333
        reference_res_2[2,2] = reference_res_2[5,5] = reference_res_2[8,8] = reference_res_2[11,11] = 1.25e-06
        reference_res_2[2,8] = reference_res_2[5,11] = reference_res_2[8,2] = reference_res_2[11,5] = 4.1667e-07
        reference_res_2[2,5] = reference_res_2[2,11] = reference_res_2[5,2] = reference_res_2[5,8] = 4.1667e-07
        reference_res_2[8,5] = reference_res_2[8,11] = reference_res_2[11,2] = reference_res_2[11,8] = 8.333e-07
        cond.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, sen_matrix, mp.ProcessInfo)
        self.assertMatrixAlmostEqual(reference_res_2, sen_matrix, 6)

        # change loading
        load_on_cond[2] =  0.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        cond.SetValue(KratosMultiphysics.PRESSURE, 1.0)

        # check w.r.t. to design variable PRESSURE
        reference_res_3 = KratosMultiphysics.Matrix(1, 12)
        reference_res_3.fill(0.0)
        reference_res_3[0,2] = reference_res_3[0,5] = reference_res_3[0,8] = reference_res_3[0,11]  = 0.25
        cond.CalculateSensitivityMatrix(KratosMultiphysics.PRESSURE, sen_matrix, mp.ProcessInfo)
        self.assertMatrixAlmostEqual(reference_res_3, sen_matrix)


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
        reference_res_1 = KratosMultiphysics.Matrix(3, 6)
        reference_res_1.fill(0.0)
        reference_res_1[0,0] = reference_res_1[0,3] = reference_res_1[1,1] = reference_res_1[1,4]  = reference_res_1[2,2] = reference_res_1[2,5] = 0.70710678
        cond.CalculateSensitivityMatrix(StructuralMechanicsApplication.LINE_LOAD, sen_matrix, mp.ProcessInfo)
        self.assertMatrixAlmostEqual(reference_res_1, sen_matrix, 6)

        # check w.r.t. to design variable SHAPE_SENSITIVITY
        reference_res_2 = KratosMultiphysics.Matrix(6, 6)
        reference_res_2.fill(0.0)
        reference_res_2[0,1] = reference_res_2[0,2] = reference_res_2[0,4] = reference_res_2[0,5] = 2499.999376
        reference_res_2[1,1] = reference_res_2[1,2] = reference_res_2[1,4] = reference_res_2[1,5] = 2499.999376
        reference_res_2[2,1] = reference_res_2[2,2] = reference_res_2[2,4] = reference_res_2[2,5] = -0.001250
        reference_res_2[3,1] = reference_res_2[3,2] = reference_res_2[3,4] = reference_res_2[3,5] = -2500.000624
        reference_res_2[4,1] = reference_res_2[4,2] = reference_res_2[4,4] = reference_res_2[4,5] = -2500.000624
        reference_res_2[5,1] = reference_res_2[5,2] = reference_res_2[5,4] = reference_res_2[5,5] = -0.001250
        cond.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, sen_matrix, mp.ProcessInfo)
        self.assertMatrixAlmostEqual(reference_res_2, sen_matrix, 6)

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
