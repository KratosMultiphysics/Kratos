from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

def add_variables(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)

def apply_material_properties(model_part, dim):
    #define properties
    model_part.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 1000)
    model_part.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 7850)
    model_part.GetProperties()[1].SetValue(StructuralMechanicsApplication.CROSS_AREA, 10)
    model_part.GetProperties()[1].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2, 5)
    g = [0,0,0]
    model_part.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION, g)
    cl = StructuralMechanicsApplication.TrussConstitutiveLaw()
    model_part.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cl)

def copy_solution_step_data_of_node(node, old_node, step=0, step_old=0):
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, step,
    old_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, step_old))

def shape_perturbation_correction_factor(node_1, node_2):
    dx = node_1.X - node_2.X
    dy = node_1.Y - node_2.Y
    dz = node_1.Z - node_2.Z
    l = math.sqrt(dx*dx + dy*dy + dz*dz)
    return l

def create_property_perturbed_elements(model_part, delta, new_element_name):
    dim = 3
    perturbed_model_part = model_part.GetModel().CreateModelPart("Property_Perturbed_Elements")
    perturbed_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dim)
    add_variables(perturbed_model_part)
    perturbed_model_part.CreateNewNode(1, model_part.Nodes[1].X, model_part.Nodes[1].Y, model_part.Nodes[1].Z)
    perturbed_model_part.CreateNewNode(2, model_part.Nodes[2].X, model_part.Nodes[2].Y, model_part.Nodes[2].Z)
    apply_material_properties(perturbed_model_part,dim)

    A_initial = model_part.GetProperties()[1][StructuralMechanicsApplication.CROSS_AREA]
    perturbed_model_part.GetProperties()[1].SetValue(StructuralMechanicsApplication.CROSS_AREA, A_initial + delta )
    prop = perturbed_model_part.GetProperties()[1]

    perturbed_model_part.CreateNewElement(new_element_name, 1, [1, 2], prop)

    for node, old_node in zip(perturbed_model_part.Nodes, model_part.Nodes):
        copy_solution_step_data_of_node(node, old_node, 0, 0)

    for element in perturbed_model_part.Elements:
        element.Initialize(perturbed_model_part.ProcessInfo)

    return perturbed_model_part

def create_shape_perturbed_elements(model_part, delta, new_element_name):
    dim=3
    perturbed_model_part = model_part.GetModel().CreateModelPart("Shape_Perturbed_Elements")
    add_variables(perturbed_model_part)

    x1 = model_part.Nodes[1].X
    y1 = model_part.Nodes[1].Y
    z1 = model_part.Nodes[1].Z
    x2 = model_part.Nodes[2].X
    y2 = model_part.Nodes[2].Y
    z2 = model_part.Nodes[2].Z
    perturbed_model_part.CreateNewNode(1, x1, y1, z1)
    perturbed_model_part.CreateNewNode(2, x1+delta, y1, z1)
    perturbed_model_part.CreateNewNode(3, x1, y1+delta, z1)
    perturbed_model_part.CreateNewNode(4, x1, y1, z1+delta)
    perturbed_model_part.CreateNewNode(5, x2, y2, z2)
    perturbed_model_part.CreateNewNode(6, x2+delta, y2, z2)
    perturbed_model_part.CreateNewNode(7, x2, y2+delta, z2)
    perturbed_model_part.CreateNewNode(8, x2, y2, z2+delta)

    apply_material_properties(perturbed_model_part, dim)
    prop = perturbed_model_part.GetProperties()[1]

    perturbed_model_part.CreateNewElement(new_element_name, 1, [2, 5], prop)
    perturbed_model_part.CreateNewElement(new_element_name, 2, [3, 5], prop)
    perturbed_model_part.CreateNewElement(new_element_name, 3, [4, 5], prop)
    perturbed_model_part.CreateNewElement(new_element_name, 4, [1, 6], prop)
    perturbed_model_part.CreateNewElement(new_element_name, 5, [1, 7], prop)
    perturbed_model_part.CreateNewElement(new_element_name, 6, [1, 8], prop)

    for i, old_node in enumerate(model_part.Nodes):
        copy_solution_step_data_of_node(perturbed_model_part.Nodes[i*4+1], old_node, 0, 0)
        copy_solution_step_data_of_node(perturbed_model_part.Nodes[i*4+2], old_node, 0, 0)
        copy_solution_step_data_of_node(perturbed_model_part.Nodes[i*4+3], old_node, 0, 0)
        copy_solution_step_data_of_node(perturbed_model_part.Nodes[i*4+4], old_node, 0, 0)

    for element in perturbed_model_part.Elements:
        element.Initialize(perturbed_model_part.ProcessInfo)

    return perturbed_model_part

def FD_calculate_sensitivity_matrix(primal_element, primal_model_part, perturbed_model_part, delta):
    # unperturbed residual
    LHS_dummy = KratosMultiphysics.Matrix()
    RHSUnperturbed = KratosMultiphysics.Vector()
    primal_element.CalculateLocalSystem(LHS_dummy, RHSUnperturbed, primal_model_part.ProcessInfo)

    # pseudo-load by finite difference approximation
    num_derivatives = perturbed_model_part.NumberOfElements()
    num_dofs = RHSUnperturbed.Size()
    FDPseudoLoadMatrix = KratosMultiphysics.Matrix(num_derivatives, num_dofs)
    RHSPerturbed = KratosMultiphysics.Vector()

    for i, element in enumerate(perturbed_model_part.Elements):
        element.CalculateLocalSystem(LHS_dummy, RHSPerturbed, perturbed_model_part.ProcessInfo)
        for j in range(num_dofs):
            FDPseudoLoadMatrix[i,j] = (RHSPerturbed[j] - RHSUnperturbed[j]) / delta

    return FDPseudoLoadMatrix

class TestTrussLinearAdjointElement(KratosUnittest.TestCase):

    def setUp(self):
        # create test model part
        dim=3
        self.current_model = KratosMultiphysics.Model()

        self.model_part = self.current_model.CreateModelPart("test")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dim)
        add_variables(self.model_part)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 1.0, 0.3)
        apply_material_properties(self.model_part, dim)
        prop = self.model_part.GetProperties()[1]

        self.model_part.CreateNewElement("AdjointFiniteDifferenceTrussLinearElement3D2N", 1, [1, 2], prop)
        self.adjoint_truss_element = self.model_part.GetElement(1)

        self.model_part.CreateNewElement("TrussLinearElement3D2N", 2, [1, 2], prop)
        self.truss_element = self.model_part.GetElement(2)

        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.014725)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.001200)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.0725715)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.019735)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.002400)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.377976)

        self.truss_element.Initialize(self.model_part.ProcessInfo)
        self.adjoint_truss_element.Initialize(self.model_part.ProcessInfo)


    def test_CalculateSensitivityMatrix_Property(self):
        # Pertubation measure
        h = 0.00001
        A_initial = self.model_part.GetProperties()[1][StructuralMechanicsApplication.CROSS_AREA]
        delta = A_initial * h

        # Create perturbed elements
        perturbed_model_part = create_property_perturbed_elements(self.model_part,delta, "TrussLinearElement3D2N")

        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_model_part, delta)

        # Pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix()
        self.model_part.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True
        self.model_part.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = h
        self.adjoint_truss_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.CROSS_AREA, PseudoLoadMatrix, self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)


    def test_CalculateSensitivityMatrix_Shape(self):
        # Pertubation measure
        h = 0.00001
        corr_factor = shape_perturbation_correction_factor(self.model_part.Nodes[1], self.model_part.Nodes[2])
        delta = corr_factor * h

        perturbed_model_part = create_shape_perturbed_elements(self.model_part,delta, "TrussLinearElement3D2N")

        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_model_part, delta)

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix()
        self.model_part.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True

        self.model_part.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = h
        self.adjoint_truss_element.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, PseudoLoadMatrix,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)


class TestTrussAdjointElement(KratosUnittest.TestCase):

    def setUp(self):
        # create test model part
        dim=3
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("test")

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dim)
        add_variables(self.model_part)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 1.0, 0.3)
        apply_material_properties(self.model_part, dim)
        prop = self.model_part.GetProperties()[1]

        self.model_part.CreateNewElement("AdjointFiniteDifferenceTrussElement3D2N", 1, [1, 2], prop)
        self.adjoint_truss_element = self.model_part.GetElement(1)

        self.model_part.CreateNewElement("TrussElement3D2N", 2, [1, 2], prop)
        self.truss_element = self.model_part.GetElement(2)

        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.14725)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.01200)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.725715)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 1.49735)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.2400)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.377976)

        self.truss_element.Initialize(self.model_part.ProcessInfo)
        self.adjoint_truss_element.Initialize(self.model_part.ProcessInfo)


    def test_CalculateSensitivityMatrix_Property(self):
        # Pertubation measure
        h = 0.00001
        A_initial = self.model_part.GetProperties()[1][StructuralMechanicsApplication.CROSS_AREA]
        delta = A_initial * h

        # Create perturbed elements
        perturbed_model_part = create_property_perturbed_elements(self.model_part,delta, "TrussElement3D2N")

        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_model_part, delta)

        # Pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix()
        self.model_part.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True

        self.model_part.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = h
        self.adjoint_truss_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.CROSS_AREA, PseudoLoadMatrix, self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)


    def test_CalculateSensitivityMatrix_Shape(self):
        # Pertubation measure
        h = 0.00001
        corr_factor = shape_perturbation_correction_factor(self.model_part.Nodes[1], self.model_part.Nodes[2])
        delta = corr_factor * h

        # Create perturbed elements
        perturbed_model_part = create_shape_perturbed_elements(self.model_part,delta, "TrussElement3D2N")

        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_model_part, delta)

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix()
        self.model_part.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True

        self.model_part.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = h
        self.adjoint_truss_element.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, PseudoLoadMatrix,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)

if __name__ == '__main__':
    KratosUnittest.main()
