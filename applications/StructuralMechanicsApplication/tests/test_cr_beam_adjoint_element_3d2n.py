from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import random
import math

class TestCrBeamAdjointElement(KratosUnittest.TestCase):

    def setUp(self):
        # create test model part
        dim=3
        self.model_part = KratosMultiphysics.ModelPart("test")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,dim)
        self._add_variables(self.model_part)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.1, 0.3)
        self._apply_material_properties(self.model_part,dim)
        prop = self.model_part.GetProperties()[0]

        self.model_part.CreateNewElement("CrLinearBeamElement3D2N", 1, [1, 2], prop)
        StructuralMechanicsApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(
            self.model_part).Execute()
        self.adjoint_beam_element = self.model_part.GetElement(1)

        self.model_part.CreateNewElement("CrLinearBeamElement3D2N", 2, [1, 2], prop)
        self.beam_element = self.model_part.GetElement(2)

        self._assign_solution_step_data(0)

    def _create_shape_disturbed_elements(self,mp,delta):
        dim=3
        self.model_part_1 = KratosMultiphysics.ModelPart("Shape_Disturbed_Elements")
        self._add_variables(self.model_part_1)

        x1 = mp.Nodes[1].X
        y1 = mp.Nodes[1].Y
        z1 = mp.Nodes[1].Z
        x2 = mp.Nodes[2].X
        y2 = mp.Nodes[2].Y
        z2 = mp.Nodes[2].Z
        self.model_part_1.CreateNewNode(1, x1, y1, z1)
        self.model_part_1.CreateNewNode(2, x1+delta, y1, z1)
        self.model_part_1.CreateNewNode(3, x1, y1+delta, z1)
        self.model_part_1.CreateNewNode(4, x1, y1, z1+delta)
        self.model_part_1.CreateNewNode(5, x2, y2, z2)
        self.model_part_1.CreateNewNode(6, x2+delta, y2, z2)
        self.model_part_1.CreateNewNode(7, x2, y2+delta, z2)
        self.model_part_1.CreateNewNode(8, x2, y2, z2+delta)

        self._apply_material_properties(self.model_part_1,dim)
        prop = self.model_part_1.GetProperties()[0]

        self.model_part_1.CreateNewElement("CrLinearBeamElement3D2N", 1, [2, 5], prop)
        self.model_part_1.CreateNewElement("CrLinearBeamElement3D2N", 2, [3, 5], prop)
        self.model_part_1.CreateNewElement("CrLinearBeamElement3D2N", 3, [4, 5], prop)
        self.model_part_1.CreateNewElement("CrLinearBeamElement3D2N", 4, [1, 6], prop)
        self.model_part_1.CreateNewElement("CrLinearBeamElement3D2N", 5, [1, 7], prop)
        self.model_part_1.CreateNewElement("CrLinearBeamElement3D2N", 6, [1, 8], prop)

    def _create_property_disturbed_elements(self,mp,delta):
        dim = 3
        self.model_part_2 = KratosMultiphysics.ModelPart("Property_Disturbed_Elements")
        self._add_variables(self.model_part_2)
        self.model_part_2.CreateNewNode(1, mp.Nodes[1].X, mp.Nodes[1].Y, mp.Nodes[1].Z)
        self.model_part_2.CreateNewNode(2, mp.Nodes[2].X, mp.Nodes[2].Y, mp.Nodes[2].Z)
        self._apply_material_properties(self.model_part_2,dim)

        I22_initial = mp.GetProperties()[0][StructuralMechanicsApplication.I22]
        self.model_part_2.GetProperties()[0].SetValue(StructuralMechanicsApplication.I22, I22_initial + delta)
        prop = self.model_part_2.GetProperties()[0]

        self.model_part_2.CreateNewElement("CrLinearBeamElement3D2N", 1, [1, 2], prop)
        self.property_disturbed_beam_element = self.model_part_2.GetElement(1)

    def _get_displacement_vector(self,mp,disp):
        index=0
        for node in mp.Nodes:
            disp[index] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0)
            index = index + 1
            disp[index] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0)
            index = index + 1
            disp[index] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
            index = index + 1
            disp[index] = node.GetSolutionStepValue(KratosMultiphysics.ROTATION_X,0)
            index = index + 1
            disp[index] = node.GetSolutionStepValue(KratosMultiphysics.ROTATION_Y,0)
            index = index + 1
            disp[index] = node.GetSolutionStepValue(KratosMultiphysics.ROTATION_Z,0)
            index = index + 1

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.30)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TORSIONAL_INERTIA,0.00001)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I22,0.00001)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I33,0.00001)

        g = [0,0,0]
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _assign_solution_step_data(self, step=0):
        # generate nodal solution step test data
        random.seed(1.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,step,random.random())
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,step,random.random())
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,step,random.random())
            node.SetSolutionStepValue(KratosMultiphysics.ROTATION_X,step,random.random())
            node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Y,step,random.random())
            node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Z,step,random.random())

    def _zero_vector(self,size):
        v = KratosMultiphysics.Vector(size)
        for i in range(size):
            v[i] = 0.0
        return v

    def _shape_disturbance_correction_factor(self):
        dx = self.model_part.Nodes[1].X - self.model_part.Nodes[2].X
        dy = self.model_part.Nodes[1].Y - self.model_part.Nodes[2].Y
        dz = self.model_part.Nodes[1].Z - self.model_part.Nodes[2].Z
        l = math.sqrt(dx*dx + dy*dy + dz*dz)
        return l

    def _assert_matrix_almost_equal(self, matrix1, matrix2, prec=7):
        self.assertEqual(matrix1.Size1(), matrix2.Size1())
        self.assertEqual(matrix1.Size2(), matrix2.Size2())
        for i in range(matrix1.Size1()):
            for j in range(matrix1.Size2()):
                self.assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

    def test_CalculateSensitivityMatrix_Shape(self):
        # unperturbed residual
        LHSUndisturbed = KratosMultiphysics.Matrix(12,12)
        RHSUndisturbed = KratosMultiphysics.Matrix(12,12)
        dummy_RHS = self._zero_vector(12)
        PrimalDisplacement = self._zero_vector(12)
        self._get_displacement_vector(self.model_part,PrimalDisplacement)

        self.beam_element.CalculateLocalSystem(LHSUndisturbed,dummy_RHS,self.model_part.ProcessInfo)
        RHSUndisturbed = LHSUndisturbed * PrimalDisplacement

        # pseudo-load by finite difference approximation
        h = 0.00001
        corr_factor = self._shape_disturbance_correction_factor()
        alpha = corr_factor * h

        FDPseudoLoadMatrix = KratosMultiphysics.Matrix(6,12)
        LHSDisturbed = KratosMultiphysics.Matrix(12,12)
        RHSDisturbed = KratosMultiphysics.Matrix(12,12)

        self._create_shape_disturbed_elements(self.model_part,alpha)

        row_index = 0
        for element in self.model_part_1.Elements:
            element.CalculateLocalSystem(LHSDisturbed,dummy_RHS,self.model_part_1.ProcessInfo)
            RHSDisturbed = LHSDisturbed * PrimalDisplacement
            for j in range(12):
                FDPseudoLoadMatrix[row_index,j] = -(RHSDisturbed[j] - RHSUndisturbed[j]) / alpha
            row_index = row_index + 1

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(6,12)
        self.adjoint_beam_element.SetValue(StructuralMechanicsApplication.DISTURBANCE_MEASURE, h)
        self.adjoint_beam_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.SHAPE,PseudoLoadMatrix,self.model_part.ProcessInfo)
        self._assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 7)

    def test_CalculateSensitivityMatrix_Property(self):
        # unperturbed residual
        LHSUndisturbed = KratosMultiphysics.Matrix(12,12)
        RHSUndisturbed = self._zero_vector(12)
        dummy_RHS = self._zero_vector(12)

        PrimalDisplacement = self._zero_vector(12)
        self._get_displacement_vector(self.model_part,PrimalDisplacement)

        self.beam_element.CalculateLocalSystem(LHSUndisturbed, dummy_RHS, self.model_part.ProcessInfo)
        RHSUndisturbed = LHSUndisturbed * PrimalDisplacement

        # pseudo-load by finite difference approximation
        h = 0.00001
        FDPseudoLoadMatrix = KratosMultiphysics.Matrix(1,12)
        LHSDisturbed = KratosMultiphysics.Matrix(12,12)
        RHSDisturbed = self._zero_vector(12)

        inital_property_value = self.model_part.GetProperties()[0][StructuralMechanicsApplication.I22]
        delta = h * inital_property_value
        self._create_property_disturbed_elements(self.model_part,delta)

        self.property_disturbed_beam_element.CalculateLocalSystem(LHSDisturbed, dummy_RHS, self.model_part_2.ProcessInfo)
        RHSDisturbed = LHSDisturbed * PrimalDisplacement

        for j in range(12):
            FDPseudoLoadMatrix[0,j] = -(RHSDisturbed[j] - RHSUndisturbed[j]) / delta

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(1,12)
        self.adjoint_beam_element.SetValue(StructuralMechanicsApplication.DISTURBANCE_MEASURE, h)
        self.adjoint_beam_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.I22, PseudoLoadMatrix, self.model_part.ProcessInfo)
        self._assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 7)

if __name__ == '__main__':
    KratosUnittest.main()
