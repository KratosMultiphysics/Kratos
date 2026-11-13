import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

# corner nodes of a distorted tetrahedron, mid-side nodes in the Kratos Tetrahedra3D10 order (0-1, 1-2, 2-0, 0-3, 1-3, 2-3),
# moved slightly off the edge mid-points so that the element is curved
CORNERS = [(0.0, 0.0, 0.0), (1.2, 0.1, 0.0), (0.2, 0.9, 0.1), (0.1, 0.2, 1.1)]
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
MID_OFFSETS = [(0.01, 0.0, 0.02), (0.0, -0.02, 0.0), (0.0, 0.0, 0.01), (-0.01, 0.0, 0.0), (0.0, 0.01, 0.0), (0.02, 0.0, 0.0)]


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)


def CreateNodes(model_part):
    for i, (x, y, z) in enumerate(CORNERS):
        model_part.CreateNewNode(i + 1, x, y, z)
    for i, ((a, b), (dx, dy, dz)) in enumerate(zip(EDGES, MID_OFFSETS)):
        p, q = CORNERS[a], CORNERS[b]
        model_part.CreateNewNode(i + 5, 0.5 * (p[0] + q[0]) + dx, 0.5 * (p[1] + q[1]) + dy, 0.5 * (p[2] + q[2]) + dz)
    for node in model_part.Nodes:                       # an arbitrary, non-rigid displacement field
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0,
                                  [0.01 * node.Id, -0.003 * node.Id * node.Id, 0.02 * node.X * node.Z + 0.001 * node.Id])


def SetProperties(properties, young_modulus):
    properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, young_modulus)
    properties.SetValue(KratosMultiphysics.POISSON_RATIO, 0.25)
    properties.SetValue(KratosMultiphysics.DENSITY, 2.4e-9)
    properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, StructuralMechanicsApplication.LinearElastic3DLaw())


class TestSmallDisplacementAdjointElement3D10N(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        AddVariables(self.model_part)
        CreateNodes(self.model_part)
        properties = self.model_part.GetProperties()[1]
        SetProperties(properties, 30000.0)
        node_ids = list(range(1, 11))
        self.adjoint_element = self.model_part.CreateNewElement("AdjointFiniteDifferencingSmallDisplacementElement3D10N", 1, node_ids, properties)
        self.primal_element = self.model_part.CreateNewElement("SmallDisplacementElement3D10N", 2, node_ids, properties)
        self.primal_element.Initialize(self.model_part.ProcessInfo)
        self.adjoint_element.Initialize(self.model_part.ProcessInfo)

    def test_CalculateLeftHandSide(self):
        primal_lhs, adjoint_lhs = KratosMultiphysics.Matrix(), KratosMultiphysics.Matrix()
        self.primal_element.CalculateLeftHandSide(primal_lhs, self.model_part.ProcessInfo)
        self.adjoint_element.CalculateLeftHandSide(adjoint_lhs, self.model_part.ProcessInfo)
        self.assertEqual(adjoint_lhs.Size1(), 30)
        self.assertMatrixAlmostEqual(primal_lhs, adjoint_lhs, 8)

    def test_CalculateSensitivityMatrix_YoungModulus(self):
        h = 1e-6
        young_modulus = self.model_part.GetProperties()[1][KratosMultiphysics.YOUNG_MODULUS]
        delta = young_modulus * h

        perturbed = self.model.CreateModelPart("perturbed")
        perturbed.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        AddVariables(perturbed)
        CreateNodes(perturbed)
        perturbed_properties = perturbed.GetProperties()[1]
        SetProperties(perturbed_properties, young_modulus + delta)
        perturbed_element = perturbed.CreateNewElement("SmallDisplacementElement3D10N", 1, list(range(1, 11)), perturbed_properties)
        perturbed_element.Initialize(perturbed.ProcessInfo)

        lhs, rhs, rhs_perturbed = KratosMultiphysics.Matrix(), KratosMultiphysics.Vector(), KratosMultiphysics.Vector()
        self.primal_element.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
        perturbed_element.CalculateLocalSystem(lhs, rhs_perturbed, perturbed.ProcessInfo)
        fd_pseudo_load = KratosMultiphysics.Matrix(1, rhs.Size())
        for j in range(rhs.Size()):
            fd_pseudo_load[0, j] = (rhs_perturbed[j] - rhs[j]) / delta

        pseudo_load = KratosMultiphysics.Matrix()
        self.model_part.ProcessInfo[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True
        self.model_part.ProcessInfo[StructuralMechanicsApplication.PERTURBATION_SIZE] = h
        self.adjoint_element.CalculateSensitivityMatrix(KratosMultiphysics.YOUNG_MODULUS, pseudo_load, self.model_part.ProcessInfo)
        self.assertEqual(pseudo_load.Size2(), 30)
        self.assertGreater(max(abs(fd_pseudo_load[0, j]) for j in range(30)), 1e-3)     # the check is not trivial
        self.assertMatrixAlmostEqual(fd_pseudo_load, pseudo_load, 5)


if __name__ == "__main__":
    KratosUnittest.main()
