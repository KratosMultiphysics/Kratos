from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

def get_displacement_vector(mp,disp):
    index=0
    for node in mp.Nodes:
        disp[index] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0)
        index = index + 1
        disp[index] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0)
        index = index + 1
        disp[index] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
        index = index + 1
 
def add_variables(mp):
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)

def apply_material_properties(mp,dim):
    #define properties
    mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,1000)
    mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,7850)
    mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.CROSS_AREA,10)
    mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,5)
    g = [0,0,0]
    mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)
    cl = StructuralMechanicsApplication.TrussConstitutiveLaw()
    mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

def copy_solution_step_data_of_node(node, mp_old, node_id, step=0):
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,step,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0))
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0))
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,
    mp_old.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0))

def shape_perturbation_correction_factor(node_1, node_2):
    dx = node_1.X - node_2.X
    dy = node_1.Y - node_2.Y
    dz = node_1.Z - node_2.Z
    l = math.sqrt(dx*dx + dy*dy + dz*dz)
    return l    
  
def zero_vector(size):
    v = KratosMultiphysics.Vector(size)
    for i in range(size):
        v[i] = 0.0
    return v

def create_property_perturbed_elements(mp, delta, new_element_name):
    dim = 3
    model_part = KratosMultiphysics.ModelPart("Property_Perturbed_Elements")
    model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,dim)
    add_variables(model_part)
    model_part.CreateNewNode(1, mp.Nodes[1].X, mp.Nodes[1].Y, mp.Nodes[1].Z)
    model_part.CreateNewNode(2, mp.Nodes[2].X, mp.Nodes[2].Y, mp.Nodes[2].Z)
    apply_material_properties(model_part,dim)

    A_initial = mp.GetProperties()[1][StructuralMechanicsApplication.CROSS_AREA]
    model_part.GetProperties()[1].SetValue(StructuralMechanicsApplication.CROSS_AREA, A_initial + delta )
    prop = model_part.GetProperties()[1]

    model_part.CreateNewElement(new_element_name, 1, [1, 2], prop)

    for i in range(2):
        copy_solution_step_data_of_node(model_part.Nodes[i+1], mp, i+1, 0)

    for element in model_part.Elements:
        element.Initialize()

    return model_part

def create_shape_perturbed_elements(mp, delta, new_element_name):
    dim=3
    model_part = KratosMultiphysics.ModelPart("Shape_Perturbed_Elements")
    add_variables(model_part)

    x1 = mp.Nodes[1].X
    y1 = mp.Nodes[1].Y
    z1 = mp.Nodes[1].Z
    x2 = mp.Nodes[2].X
    y2 = mp.Nodes[2].Y
    z2 = mp.Nodes[2].Z
    model_part.CreateNewNode(1, x1, y1, z1)
    model_part.CreateNewNode(2, x1+delta, y1, z1)
    model_part.CreateNewNode(3, x1, y1+delta, z1)
    model_part.CreateNewNode(4, x1, y1, z1+delta)
    model_part.CreateNewNode(5, x2, y2, z2)
    model_part.CreateNewNode(6, x2+delta, y2, z2)
    model_part.CreateNewNode(7, x2, y2+delta, z2)
    model_part.CreateNewNode(8, x2, y2, z2+delta)

    apply_material_properties(model_part,dim)
    prop = model_part.GetProperties()[1]

    model_part.CreateNewElement(new_element_name, 1, [2, 5], prop)
    model_part.CreateNewElement(new_element_name, 2, [3, 5], prop)
    model_part.CreateNewElement(new_element_name, 3, [4, 5], prop)
    model_part.CreateNewElement(new_element_name, 4, [1, 6], prop)
    model_part.CreateNewElement(new_element_name, 5, [1, 7], prop)
    model_part.CreateNewElement(new_element_name, 6, [1, 8], prop)

    index = 1
    for i in range(2):
        for _ in range(4):
            copy_solution_step_data_of_node(model_part.Nodes[index], mp, i + 1, 0)
            index += 1

    for element in model_part.Elements:
        element.Initialize()

    return model_part

def FD_calculate_sensitivity_matrix(primal_element, primal_mp, perturbed_mp, delta):
    # unperturbed residual
    LHS_dummy = KratosMultiphysics.Matrix(6,6)
    RHSUnperturbed = zero_vector(6)
    primal_element.CalculateLocalSystem(LHS_dummy, RHSUnperturbed, primal_mp.ProcessInfo)

    # pseudo-load by finite difference approximation
    num_derivatives = len(perturbed_mp.Elements)
    FDPseudoLoadMatrix = KratosMultiphysics.Matrix(num_derivatives,6)
    RHSPerturbed = zero_vector(6)

    row_index = 0
    for element in perturbed_mp.Elements:
        element.CalculateLocalSystem(LHS_dummy, RHSPerturbed, perturbed_mp.ProcessInfo)
        for j in range(6):
            FDPseudoLoadMatrix[row_index,j] = (RHSPerturbed[j] - RHSUnperturbed[j]) / delta
        row_index = row_index + 1

    return FDPseudoLoadMatrix

def assert_matrix_almost_equal(matrix1, matrix2, prec=7):
    KratosUnittest.TestCase().assertEqual(matrix1.Size1(), matrix2.Size1())
    KratosUnittest.TestCase().assertEqual(matrix1.Size2(), matrix2.Size2())
    for i in range(matrix1.Size1()):
        for j in range(matrix1.Size2()):
            KratosUnittest.TestCase().assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

class TestTrussLinearAdjointElement(KratosUnittest.TestCase):

    def setUp(self):
        # create test model part
        dim=3
        self.model_part = KratosMultiphysics.ModelPart("test")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,dim)
        add_variables(self.model_part)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 1.0, 0.3)
        apply_material_properties(self.model_part,dim)
        prop = self.model_part.GetProperties()[1]

        self.model_part.CreateNewElement("TrussLinearElement3D2N", 1, [1, 2], prop)
        StructuralMechanicsApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(
            self.model_part).Execute()
        self.adjoint_truss_element = self.model_part.GetElement(1)

        self.model_part.CreateNewElement("TrussLinearElement3D2N", 2, [1, 2], prop)
        self.truss_element = self.model_part.GetElement(2)

        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.014725)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.001200)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.0725715)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.019735)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.002400)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.377976)

        self.truss_element.Initialize()
        self.adjoint_truss_element.Initialize()


    def test_CalculateSensitivityMatrix_Property(self):
        # Pertubation measure
        h = 0.00001
        A_initial = self.model_part.GetProperties()[1][StructuralMechanicsApplication.CROSS_AREA]
        delta = A_initial * h

        # Create perturbed elements
        perturbed_mp = create_property_perturbed_elements(self.model_part,delta,"TrussLinearElement3D2N")
        
        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_mp, delta)  
        
        # Pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(1,6)
        self.adjoint_truss_element.SetValue(StructuralMechanicsApplication.PERTURBATION_SIZE, h)
        self.adjoint_truss_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.CROSS_AREA, PseudoLoadMatrix, self.model_part.ProcessInfo)
        assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)


    def test_CalculateSensitivityMatrix_Shape(self):
        # Pertubation measure
        h = 0.00001
        corr_factor = shape_perturbation_correction_factor(self.model_part.Nodes[1], self.model_part.Nodes[2])
        delta = corr_factor * h

        perturbed_mp = create_shape_perturbed_elements(self.model_part,delta,"TrussLinearElement3D2N")

        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_mp, delta)  

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(6,6)
        self.adjoint_truss_element.SetValue(StructuralMechanicsApplication.PERTURBATION_SIZE, h)
        self.adjoint_truss_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.SHAPE,PseudoLoadMatrix,self.model_part.ProcessInfo)
        assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)


class TestTrussAdjointElement(KratosUnittest.TestCase):

    def setUp(self):
        # create test model part
        dim=3
        self.model_part = KratosMultiphysics.ModelPart("test")
       
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,dim)
        add_variables(self.model_part)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 1.0, 0.3)
        apply_material_properties(self.model_part,dim)
        prop = self.model_part.GetProperties()[1]

        self.model_part.CreateNewElement("TrussElement3D2N", 1, [1, 2], prop)
        StructuralMechanicsApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(
            self.model_part).Execute()
        self.adjoint_truss_element = self.model_part.GetElement(1)

        self.model_part.CreateNewElement("TrussElement3D2N", 2, [1, 2], prop)
        self.truss_element = self.model_part.GetElement(2)
      
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.14725)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.01200)
        self.model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.725715)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 1.49735)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.2400)
        self.model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.377976)

        self.truss_element.Initialize()
        self.adjoint_truss_element.Initialize()


    def test_CalculateSensitivityMatrix_Property(self):
        # Pertubation measure
        h = 0.00001
        A_initial = self.model_part.GetProperties()[1][StructuralMechanicsApplication.CROSS_AREA]
        delta = A_initial * h

        # Create perturbed elements
        perturbed_mp = create_property_perturbed_elements(self.model_part,delta,"TrussElement3D2N")
        
        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_mp, delta)  
        
        # Pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(1,6)
        self.adjoint_truss_element.SetValue(StructuralMechanicsApplication.PERTURBATION_SIZE, h)
        self.adjoint_truss_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.CROSS_AREA, PseudoLoadMatrix, self.model_part.ProcessInfo)
        assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)


    def test_CalculateSensitivityMatrix_Shape(self):
        # Pertubation measure
        h = 0.00001
        corr_factor = shape_perturbation_correction_factor(self.model_part.Nodes[1], self.model_part.Nodes[2])
        delta = corr_factor * h

        # Create perturbed elements
        perturbed_mp = create_shape_perturbed_elements(self.model_part,delta,"TrussElement3D2N")

        # Derive RHS by finite differences
        FDPseudoLoadMatrix = FD_calculate_sensitivity_matrix(self.truss_element, self.model_part, perturbed_mp, delta)  

        # pseudo-load computation by adjoint element
        PseudoLoadMatrix = KratosMultiphysics.Matrix(6,6)
        self.adjoint_truss_element.SetValue(StructuralMechanicsApplication.PERTURBATION_SIZE, h)
        self.adjoint_truss_element.CalculateSensitivityMatrix(StructuralMechanicsApplication.SHAPE,PseudoLoadMatrix,self.model_part.ProcessInfo)
        assert_matrix_almost_equal(FDPseudoLoadMatrix, PseudoLoadMatrix, 5)

if __name__ == '__main__':
    KratosUnittest.main()
