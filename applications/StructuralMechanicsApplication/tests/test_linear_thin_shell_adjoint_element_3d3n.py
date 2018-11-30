from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

def zero_vector(size):
    v = KratosMultiphysics.Vector(size)
    for i in range(size):
        v[i] = 0.0
    return v

def add_variables(mp):
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

def apply_material_properties(mp,dim):
    #define properties
    mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,1000)
    mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
    mp.GetProperties()[0].SetValue(KratosMultiphysics.THICKNESS,1.0)
    mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,1.0)
    g = [0,0,0]
    mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)
    cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
    mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

def copy_solution_step_data_of_node(node, mp_old, node_id, step=0):
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,step,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0))
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0))
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0))
    node.SetSolutionStepValue(KratosMultiphysics.ROTATION_X,step,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.ROTATION_X,0))
    node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Y,step,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.ROTATION_Y,0))
    node.SetSolutionStepValue(KratosMultiphysics.ROTATION_Z,step,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.ROTATION_Z,0))

class TestShellThinAdjointElement3D3N(KratosUnittest.TestCase):
    def setUp(self):
        # create test model part
        dim=3
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("test")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,dim)
        add_variables(self.model_part)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 3.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 0.0, 4.0, 0.0)
        apply_material_properties(self.model_part,dim)
        prop = self.model_part.GetProperties()[0]

        self.model_part.CreateNewElement("ShellThinElement3D3N", 1, [1, 2, 3], prop)
        StructuralMechanicsApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(
            self.model_part).Execute()
        self.adjoint_shell_element = self.model_part.GetElement(1)

        self.model_part.CreateNewElement("ShellThinElement3D3N", 2, [1, 2, 3], prop)
        self.shell_element = self.model_part.GetElement(2)

        self._assign_solution_step_data(0)

        self.shell_element.Initialize()
        self.adjoint_shell_element.Initialize()

    def _create_shape_perturbed_elements(self,mp,delta):
        dim=3
        self.model_part_1 = mp.GetModel().CreateModelPart("Shape_Perturbed_Elements")
        add_variables(self.model_part_1)

        x1 = mp.Nodes[1].X
        y1 = mp.Nodes[1].Y
        z1 = mp.Nodes[1].Z
        x2 = mp.Nodes[2].X
        y2 = mp.Nodes[2].Y
        z2 = mp.Nodes[2].Z
        x3 = mp.Nodes[3].X
        y3 = mp.Nodes[3].Y
        z3 = mp.Nodes[3].Z
        self.model_part_1.CreateNewNode(1, x1, y1, z1)
        self.model_part_1.CreateNewNode(2, x1+delta, y1, z1)
        self.model_part_1.CreateNewNode(3, x1, y1+delta, z1)
        self.model_part_1.CreateNewNode(4, x1, y1, z1+delta)
        self.model_part_1.CreateNewNode(5, x2, y2, z2)
        self.model_part_1.CreateNewNode(6, x2+delta, y2, z2)
        self.model_part_1.CreateNewNode(7, x2, y2+delta, z2)
        self.model_part_1.CreateNewNode(8, x2, y2, z2+delta)
        self.model_part_1.CreateNewNode(9, x3, y3, z3)
        self.model_part_1.CreateNewNode(10, x3+delta, y3, z3)
        self.model_part_1.CreateNewNode(11, x3, y3+delta, z3)
        self.model_part_1.CreateNewNode(12, x3, y3, z3+delta)

        apply_material_properties(self.model_part_1,dim)
        prop = self.model_part_1.GetProperties()[0]

        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 1, [2, 5, 9], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 2, [3, 5, 9], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 3, [4, 5, 9], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 4, [1, 6, 9], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 5, [1, 7, 9], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 6, [1, 8, 9], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 7, [1, 5, 10], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 8, [1, 5, 11], prop)
        self.model_part_1.CreateNewElement("ShellThinElement3D3N", 9, [1, 5, 12], prop)

        for element in self.model_part_1.Elements:
            element.Initialize()

        index = 1
        for i in range(3):
            for _ in range(4):
                copy_solution_step_data_of_node(self.model_part_1.Nodes[index], self.model_part, i + 1, 0)
                index += 1

    def _create_property_perturbed_elements(self,mp,delta):
        dim = 3
        self.model_part_2 = mp.GetModel().CreateModelPart("Property_Perturbed_Elements")
        add_variables(self.model_part_2)
        self.model_part_2.CreateNewNode(1, mp.Nodes[1].X, mp.Nodes[1].Y, mp.Nodes[1].Z)
        self.model_part_2.CreateNewNode(2, mp.Nodes[2].X, mp.Nodes[2].Y, mp.Nodes[2].Z)
        self.model_part_2.CreateNewNode(3, mp.Nodes[3].X, mp.Nodes[3].Y, mp.Nodes[3].Z)
        apply_material_properties(self.model_part_2,dim)

        THICKNESS_initial = mp.GetProperties()[0][KratosMultiphysics.THICKNESS]
        self.model_part_2.GetProperties()[0].SetValue(KratosMultiphysics.THICKNESS, THICKNESS_initial + delta)
        prop = self.model_part_2.GetProperties()[0]

        self.model_part_2.CreateNewElement("ShellThinElement3D3N", 1, [1, 2, 3], prop)
        self.property_perturbed_shell_element = self.model_part_2.GetElement(1)

        for i in range(3):
            copy_solution_step_data_of_node(self.model_part_2.Nodes[i+1], self.model_part, i+1, 0)

        self.property_perturbed_shell_element.Initialize()

    def _assign_solution_step_data(self, step=0):
        # generate nodal solution step test data
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,step,0.014725)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,step,0.001200)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,step,0.0725715)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.ROTATION_X,step,0.00125)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.ROTATION_Y,step,-0.114905)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.ROTATION_Z,step,0.258032)

        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,step,0.019735)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,step,0.002400)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,step,0.377976)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.ROTATION_X,step,-0.00155)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.ROTATION_Y,step,-0.217714)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.ROTATION_Z,step,0.2544032)

        self.model_part.Nodes[3].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,step,0.027725)
        self.model_part.Nodes[3].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,step,0.05803)
        self.model_part.Nodes[3].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,step,0.0875715)
        self.model_part.Nodes[3].SetSolutionStepValue(KratosMultiphysics.ROTATION_X,step,0.035892)
        self.model_part.Nodes[3].SetSolutionStepValue(KratosMultiphysics.ROTATION_Y,step,-0.784905)
        self.model_part.Nodes[3].SetSolutionStepValue(KratosMultiphysics.ROTATION_Z,step,0.4568742)

    def _shape_perturbation_correction_factor(self):
        dx = self.model_part.Nodes[1].X - self.model_part.Nodes[2].X
        dy = self.model_part.Nodes[1].Y - self.model_part.Nodes[2].Y
        dz = self.model_part.Nodes[1].Z - self.model_part.Nodes[2].Z
        l = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx = self.model_part.Nodes[1].X - self.model_part.Nodes[3].X
        dy = self.model_part.Nodes[1].Y - self.model_part.Nodes[3].Y
        dz = self.model_part.Nodes[1].Z - self.model_part.Nodes[3].Z
        l += math.sqrt(dx*dx + dy*dy + dz*dz)
        dx = self.model_part.Nodes[2].X - self.model_part.Nodes[3].X
        dy = self.model_part.Nodes[2].Y - self.model_part.Nodes[3].Y
        dz = self.model_part.Nodes[2].Z - self.model_part.Nodes[3].Z
        l += math.sqrt(dx*dx + dy*dy + dz*dz)
        l = l/3
        return l

    def _assert_matrix_almost_equal(self, matrix1, matrix2, prec=4):
        self.assertEqual(matrix1.Size1(), matrix2.Size1())
        self.assertEqual(matrix1.Size2(), matrix2.Size2())
        for i in range(matrix1.Size1()):
            for j in range(matrix1.Size2()):
                self.assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

    def test_CalculateSensitivityMatrix_Shape(self):
        # unperturbed residual
        dummy_LHS = KratosMultiphysics.Matrix(18,18)
        RHSUnperturbed = zero_vector(18)

        self.shell_element.CalculateLocalSystem(dummy_LHS, RHSUnperturbed, self.model_part.ProcessInfo)

        # pseudo-load by finite difference approximation
        h = 0.000001
        corr_factor = self._shape_perturbation_correction_factor()
        alpha = corr_factor * h

        FDPseudoLoadMatrix = KratosMultiphysics.Matrix(9,18)
        dummy_LHS = KratosMultiphysics.Matrix(18,18)
        RHSPerturbed = zero_vector(18)

        self._create_shape_perturbed_elements(self.model_part,alpha)
        row_index = 0
        for element in self.model_part_1.Elements:
            element.CalculateLocalSystem(dummy_LHS, RHSPerturbed, self.model_part_1.ProcessInfo)
            for j in range(18):
                FDPseudoLoadMatrix[row_index,j] = (RHSPerturbed[j] - RHSUnperturbed[j]) / alpha
            row_index = row_index + 1

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(9,18)
        self.adjoint_shell_element.SetValue(StructuralMechanicsApplication.PERTURBATION_SIZE, h)
        self.adjoint_shell_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.SHAPE,PseudoLoadMatrix,self.model_part.ProcessInfo)
        self._assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 4)

    def test_CalculateSensitivityMatrix_Property(self):
        # unperturbed residual
        dummy_LHS = KratosMultiphysics.Matrix(18,18)
        RHSUnperturbed = zero_vector(18)

        self.shell_element.CalculateLocalSystem(dummy_LHS, RHSUnperturbed, self.model_part.ProcessInfo)

        # pseudo-load by finite difference approximation
        h = 0.00001
        FDPseudoLoadMatrix = KratosMultiphysics.Matrix(1,18)
        RHSPerturbed = zero_vector(18)

        inital_property_value = self.model_part.GetProperties()[0][KratosMultiphysics.THICKNESS]
        delta = h * inital_property_value
        self._create_property_perturbed_elements(self.model_part,delta)

        self.property_perturbed_shell_element.CalculateLocalSystem(dummy_LHS, RHSPerturbed, self.model_part_2.ProcessInfo)

        for j in range(18):
            FDPseudoLoadMatrix[0,j] = (RHSPerturbed[j] - RHSUnperturbed[j]) / delta

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(1,18)
        self.adjoint_shell_element.SetValue(StructuralMechanicsApplication.PERTURBATION_SIZE, h)
        self.adjoint_shell_element.CalculateSensitivityMatrix(KratosMultiphysics.THICKNESS, PseudoLoadMatrix, self.model_part.ProcessInfo)
        self._assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 4)

if __name__ == '__main__':
    KratosUnittest.main()
