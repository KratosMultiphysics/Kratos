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


APEX = (0.3, 0.3, -1.0)
# element 2 of the two-element mesh: corners 2,1,3 (shared face) and the apex 11, then its mid-side nodes
ELEMENT_2_NODES = [2, 1, 3, 11, 5, 7, 6, 13, 12, 14]
NUM_DOFS = 42
APEX_MID_OFFSETS = [(0.01, 0.0, 0.0), (0.0, -0.01, 0.0), (0.0, 0.0, 0.02)]


def CreateTwoElementNodes(model_part):
    """Nodes 1-10 as in CreateNodes plus the apex 11 and its mid-side nodes 12-14 (corner 1/2/3 to the apex)."""
    CreateNodes(model_part)
    model_part.CreateNewNode(11, *APEX)
    for i, (dx, dy, dz) in enumerate(APEX_MID_OFFSETS):
        p = CORNERS[i]
        model_part.CreateNewNode(12 + i, 0.5 * (p[0] + APEX[0]) + dx, 0.5 * (p[1] + APEX[1]) + dy, 0.5 * (p[2] + APEX[2]) + dz)
    for node in model_part.Nodes:
        if node.Id > 10:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0,
                                      [0.01 * node.Id, -0.003 * node.Id * node.Id, 0.02 * node.X * node.Z + 0.001 * node.Id])


def CreateTwoElementPrimal(model_part, properties):
    elements = [model_part.CreateNewElement("SmallDisplacementElement3D10N", 1, list(range(1, 11)), properties),
                model_part.CreateNewElement("SmallDisplacementElement3D10N", 2, ELEMENT_2_NODES, properties)]
    for element in elements:
        element.Initialize(model_part.ProcessInfo)
    return elements


class TestSmallDisplacementAdjointElement3D10NTwoElements(KratosUnittest.TestCase):
    """Two tetrahedra sharing the face (1,2,3): checks sensitivities that touch shared nodes."""
    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        AddVariables(self.model_part)
        CreateTwoElementNodes(self.model_part)
        self.properties = self.model_part.GetProperties()[1]
        SetProperties(self.properties, 30000.0)
        self.primal_elements = CreateTwoElementPrimal(self.model_part, self.properties)
        self.adjoint_elements = [
            self.model_part.CreateNewElement("AdjointFiniteDifferencingSmallDisplacementElement3D10N", 3, list(range(1, 11)), self.properties),
            self.model_part.CreateNewElement("AdjointFiniteDifferencingSmallDisplacementElement3D10N", 4, ELEMENT_2_NODES, self.properties)]
        for element in self.adjoint_elements:
            element.Initialize(self.model_part.ProcessInfo)

    def _GlobalRhs(self, elements):
        rhs_global = [0.0] * NUM_DOFS
        rhs = KratosMultiphysics.Vector()
        for element in elements:
            element.CalculateRightHandSide(rhs, self.model_part.ProcessInfo)
            for a, node in enumerate(element.GetGeometry()):
                for d in range(3):
                    rhs_global[3 * (node.Id - 1) + d] += rhs[3 * a + d]
        return rhs_global

    def test_ElementsHavePositiveOrientation(self):
        for element in self.primal_elements:
            p = [[n.X, n.Y, n.Z] for n in list(element.GetGeometry())[:4]]
            a, b, c = [[p[i][k] - p[0][k] for k in range(3)] for i in (1, 2, 3)]
            det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
            self.assertGreater(det / 6.0, 0.05)

    def test_CalculateLeftHandSide_Displacement(self):
        lhs_global = [[0.0] * NUM_DOFS for _ in range(NUM_DOFS)]
        lhs = KratosMultiphysics.Matrix()
        contributions = {}      # (global dof row, col) -> number of elements that contribute
        for element in self.adjoint_elements:
            element.CalculateLeftHandSide(lhs, self.model_part.ProcessInfo)
            self.assertEqual(lhs.Size1(), 30)
            ids = [node.Id for node in element.GetGeometry()]
            for a in range(30):
                for b in range(30):
                    row, col = 3 * (ids[a // 3] - 1) + a % 3, 3 * (ids[b // 3] - 1) + b % 3
                    lhs_global[row][col] += lhs[a, b]
                    contributions[(row, col)] = contributions.get((row, col), 0) + 1

        # shared node 1 (dofs 0-2): its diagonal block receives contributions from both elements
        self.assertEqual(contributions[(0, 0)], 2)
        self.assertEqual(contributions[(1, 1)], 2)
        self.assertEqual(contributions[(3 * (5 - 1), 3 * (5 - 1))], 2)         # shared mid-side node 5 (id 5)

        # central FD of the summed primal RHS w.r.t. every nodal displacement component (linear problem, so exact up to round-off)
        delta = 1e-3
        fd_lhs = [[0.0] * NUM_DOFS for _ in range(NUM_DOFS)]
        for node in self.model_part.Nodes:
            for d in range(3):
                u = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                u_plus, u_minus = list(u), list(u)
                u_plus[d] += delta
                u_minus[d] -= delta
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, u_plus)
                rhs_plus = self._GlobalRhs(self.primal_elements)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, u_minus)
                rhs_minus = self._GlobalRhs(self.primal_elements)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, u)
                col = 3 * (node.Id - 1) + d
                for row in range(NUM_DOFS):
                    fd_lhs[row][col] = -(rhs_plus[row] - rhs_minus[row]) / (2.0 * delta)

        max_k = max(abs(v) for row in lhs_global for v in row)
        self.assertGreater(max_k, 1.0)
        max_error = max(abs(lhs_global[i][j] - fd_lhs[i][j]) for i in range(NUM_DOFS) for j in range(NUM_DOFS))
        self.assertLess(max_error, 1e-11 * max_k)

    def _CreatePerturbedModelPart(self, name, node_id, direction, alpha):
        perturbed = self.model.CreateModelPart(name)
        perturbed.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        AddVariables(perturbed)
        CreateTwoElementNodes(perturbed)
        properties = perturbed.GetProperties()[1]
        SetProperties(properties, 30000.0)
        node = perturbed.GetNode(node_id)
        # the element is evaluated in the initial configuration, so move X0 as well (as the adjoint FD utility does)
        initial = [node.X0, node.Y0, node.Z0]
        initial[direction] += alpha
        node.X0, node.Y0, node.Z0 = initial
        coords = [node.X, node.Y, node.Z]
        coords[direction] += alpha
        node.X, node.Y, node.Z = coords
        return perturbed, CreateTwoElementPrimal(perturbed, properties)

    def test_CalculateSensitivityMatrix_Shape(self):
        h = 1e-6
        process_info = self.model_part.ProcessInfo
        process_info[StructuralMechanicsApplication.ADAPT_PERTURBATION_SIZE] = True
        process_info[StructuralMechanicsApplication.PERTURBATION_SIZE] = h

        rhs0, rhs_plus, rhs_minus = KratosMultiphysics.Vector(), KratosMultiphysics.Vector(), KratosMultiphysics.Vector()
        nodes_before = [(n.X, n.Y, n.Z, n.X0, n.Y0, n.Z0) for n in self.model_part.Nodes]
        rhs_before = []
        for element in self.primal_elements:
            element.CalculateRightHandSide(rhs0, process_info)
            rhs_before.append([rhs0[j] for j in range(30)])

        counter = 0
        for e, adjoint_element in enumerate(self.adjoint_elements):
            # same step as GetPerturbationSize: PERTURBATION_SIZE * element domain size
            alpha = h * self.primal_elements[e].GetGeometry().DomainSize()
            alpha_ref = 100.0 * alpha       # larger central reference, so the check measures the C++ forward-FD truncation, not only plumbing
            pseudo_load = KratosMultiphysics.Matrix()
            adjoint_element.CalculateSensitivityMatrix(KratosMultiphysics.SHAPE_SENSITIVITY, pseudo_load, process_info)
            self.assertEqual(pseudo_load.Size1(), 30)
            self.assertEqual(pseudo_load.Size2(), 30)
            max_value = 0.0
            max_error = 0.0
            for a, node in enumerate(adjoint_element.GetGeometry()):
                for d in range(3):
                    counter += 1
                    plus, elements_plus = self._CreatePerturbedModelPart("perturbed_p_" + str(counter), node.Id, d, alpha_ref)
                    minus, elements_minus = self._CreatePerturbedModelPart("perturbed_m_" + str(counter), node.Id, d, -alpha_ref)
                    elements_plus[e].CalculateRightHandSide(rhs_plus, plus.ProcessInfo)
                    elements_minus[e].CalculateRightHandSide(rhs_minus, minus.ProcessInfo)
                    for j in range(30):
                        fd = (rhs_plus[j] - rhs_minus[j]) / (2.0 * alpha_ref)
                        max_value = max(max_value, abs(fd))
                        max_error = max(max_error, abs(fd - pseudo_load[3 * a + d, j]))
            self.assertGreater(max_value, 1.0)      # the check is not trivial
            self.assertLess(max_error, 1e-6 * max_value)

        # the +=/-= perturbation of the shared nodes must be undone (up to round-off), in X and in X0 (the element reads X0)
        for node, before in zip(self.model_part.Nodes, nodes_before):
            for value, reference in zip((node.X, node.Y, node.Z, node.X0, node.Y0, node.Z0), before):
                self.assertAlmostEqual(value, reference, 12)
        for element, expected in zip(self.primal_elements, rhs_before):
            element.CalculateRightHandSide(rhs0, process_info)
            max_rhs = max(abs(v) for v in expected)
            for j in range(30):
                self.assertLess(abs(rhs0[j] - expected[j]), 1e-9 * max_rhs)


if __name__ == "__main__":
    KratosUnittest.main()
