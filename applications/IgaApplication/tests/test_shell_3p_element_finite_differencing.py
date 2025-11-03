#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os


def _report_matrix_diff(reference_matrix, comparison_matrix, diff_threshold=1e-4,
                        reference_label="lhs", comparison_label="adjoint_LHS"):
    """Print matrix entries whose absolute difference exceeds diff_threshold."""
    if reference_matrix.Size1() != comparison_matrix.Size1() or reference_matrix.Size2() != comparison_matrix.Size2():
        raise RuntimeError(f"{reference_label} and {comparison_label} have different shapes")

    large_differences = []
    for row in range(reference_matrix.Size1()):
        for col in range(reference_matrix.Size2()):
            diff = abs(reference_matrix[row, col] - comparison_matrix[row, col])
            if diff > diff_threshold:
                large_differences.append((row, col, reference_matrix[row, col], comparison_matrix[row, col], diff))

    if large_differences:
        print(f"Entries with |{reference_label} - {comparison_label}| > {diff_threshold}:")
        for row, col, ref_val, cmp_val, diff in large_differences:
            print(f"  row {row}, col {col}: {reference_label}={ref_val}, {comparison_label}={cmp_val}, diff={diff}")
    else:
        print(f"{reference_label.upper()} vs {comparison_label.upper()} - No entries differ by more than {diff_threshold}")


def _report_vector_diff(reference_vector, comparison_vector, diff_threshold=1e-4,
                        reference_label="primal_rhs", comparison_label="adjoint_rhs"):
    """Print vector entries whose absolute difference exceeds diff_threshold."""
    if reference_vector.Size() != comparison_vector.Size():
        raise RuntimeError(f"{reference_label} and {comparison_label} have different sizes")

    large_differences = []
    for idx in range(reference_vector.Size()):
        diff = abs(reference_vector[idx] - comparison_vector[idx])
        if diff > diff_threshold:
            large_differences.append((idx, reference_vector[idx], comparison_vector[idx], diff))

    if large_differences:
        print(f"Entries with |{reference_label} - {comparison_label}| > {diff_threshold}:")
        for idx, ref_val, cmp_val, diff in large_differences:
            print(f"  index {idx}: {reference_label}={ref_val}, {comparison_label}={cmp_val}, diff={diff}")
    else:
        print(f"{reference_label.upper()} vs {comparison_label.upper()} - No entries differ by more than {diff_threshold}")


def solve_cantilever(create_geometry):
    model = KM.Model()
    model_part = model.CreateModelPart('Model')

    model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KM.REACTION)
    model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

    # create property for shell elements

    shell_properties = model_part.GetProperties()[1]
    shell_properties.SetValue(KM.THICKNESS, 0.1)
    shell_properties.SetValue(KM.YOUNG_MODULUS, 100000)
    shell_properties.SetValue(KM.POISSON_RATIO, 0)
    shell_properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElasticPlaneStress2DLaw())

    # create a nurbs surface
    surface = create_geometry(model_part)

    # create quadrature_point_geometries
    quadrature_point_geometries = KM.GeometriesVector()

    surface.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)

    element_id = 1
    for i in range(0, len(quadrature_point_geometries)):
        model_part.CreateNewElement('ActiveShell3pElement', element_id, quadrature_point_geometries[i], shell_properties)
        model_part.CreateNewElement('ActiveAdjointFiniteDifferenceBaseElement', element_id + 100, quadrature_point_geometries[i], shell_properties)
        model_part.CreateNewElement('Shell3pElement', element_id + 1000, quadrature_point_geometries[i], shell_properties)
        element_id += 1


    params = KM.Parameters("""{
        "iga_model_part_name"         : "Model",
        "active_shell_model_part_name": "ACTIVE_MP",
        "applied_actuation_list" : [
"alpha",
"beta",
"gamma",
"kappa_1",
"kappa_2",
"kappa_12"
                           ],
        "applied_actuation_value": [
0.50,
0.00,
0.00,
0.00,
0.00,
0.00
        ],
        "unfixed_actuation_list": [
"fix",
"fix",
"fix",
"fix",
"fix",
"fix"
                           ]
    }""")
    IGA.ActiveShellElementDofAssignmentProcess(model, params).ExecuteInitialize()

    # divisor = 50
    # for node in model_part.Nodes:
    #     node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Array3([node.Id / divisor, (node.Id + 1) / divisor, (node.Id + 2) / divisor]))

    # add dofs
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

    
    print(model_part)
    delta = 1e-6
    for element in model_part.Elements:
        lhs = KM.Matrix()
        primal_rhs = KM.Vector()

        element.Initialize(model_part.ProcessInfo)
        element.InitializeNonLinearIteration(model_part.ProcessInfo)

        if element.Id < 100:
            element.CalculateLocalSystem(lhs, primal_rhs, model_part.ProcessInfo)

            print("primal_rhs\n:", primal_rhs)
            print("-" * 40)
            print("lhs\n:", lhs)
            print("-" * 40)


            # adjoint calculation
            adjoint_element: KM.Element = model_part.GetElement(element.Id + 100)
            adjoint_LHS = KM.Matrix()
            adjoint_element.Initialize(model_part.ProcessInfo)
            adjoint_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            adjoint_element.CalculateLeftHandSide(adjoint_LHS, model_part.ProcessInfo)

            adjoint_rhs = KM.Vector()
            adjoint_element.CalculateRightHandSide(adjoint_rhs, model_part.ProcessInfo)
            # adjoint_element.CalculateSensitivityMatrix(adjoint_LHS, model_part.ProcessInfo)
            

            #normal Shell3pElement
            Shell3p_element: KM.Element = model_part.GetElement(element.Id + 1000)
            Shell3p_LHS = KM.Matrix()
            Shell3p_element.Initialize(model_part.ProcessInfo)
            Shell3p_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            Shell3p_element.CalculateLeftHandSide(Shell3p_LHS, model_part.ProcessInfo)

            print("size Shell3p_LHS:", Shell3p_LHS.Size1(), Shell3p_LHS.Size2())

            _report_matrix_diff(lhs, adjoint_LHS, diff_threshold=1e-4,
                                reference_label="lhs", comparison_label="adjoint_LHS")

            _report_vector_diff(primal_rhs, adjoint_rhs, diff_threshold=1e-4,
                                reference_label="primal_rhs", comparison_label="adjoint_rhs")


            perturbed_rhs = KM.Vector()
            shell3p_reference_rhs = KM.Vector()
            Shell3p_element.CalculateLocalSystem(Shell3p_LHS, shell3p_reference_rhs, model_part.ProcessInfo)
            shell3p_rhs = KM.Vector()

            for node_idx, node in enumerate(element.GetGeometry()):
                node.X += delta

                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(lhs, perturbed_rhs, model_part.ProcessInfo)
                sensitivities = (perturbed_rhs - primal_rhs) / delta

                Shell3p_element.InitializeNonLinearIteration(model_part.ProcessInfo)
                Shell3p_element.CalculateLocalSystem(Shell3p_LHS, shell3p_rhs, model_part.ProcessInfo)
                shell3p_sensitivities = (shell3p_rhs - shell3p_reference_rhs) / delta

                node.X -= delta

                print(f"perturbed node {node_idx + 1}")
                print(perturbed_rhs)

                max_dofs = min(27, sensitivities.Size(), shell3p_sensitivities.Size())
                for dof_idx in range(max_dofs):
                    print(sensitivities[dof_idx], shell3p_sensitivities[dof_idx], adjoint_LHS[dof_idx * 3, dof_idx])




    return surface

def create_geometry(model_part):
    node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0)
    node2 = model_part.CreateNewNode(2, 0.5, 0.0, 0)
    node3 = model_part.CreateNewNode(3, 1.0, 0.0, 0)

    node4 = model_part.CreateNewNode(4, 0.0, 0.5, 0)
    node5 = model_part.CreateNewNode(5, 0.5, 0.5, 0)
    node6 = model_part.CreateNewNode(6, 1.0, 0.5, 0)

    node7 = model_part.CreateNewNode(7, 0.0, 1.0, 0)
    node8 = model_part.CreateNewNode(8, 0.5, 1.0, 0)
    node9 = model_part.CreateNewNode(9, 1.0, 1.0, 0)

    nodes = KM.NodesVector()
    nodes.append(node1)
    nodes.append(node2)
    nodes.append(node3)
    nodes.append(node4)
    nodes.append(node5)
    nodes.append(node6)
    nodes.append(node7)
    nodes.append(node8)
    nodes.append(node9)

    knots_u = KM.Vector(4)
    knots_u[0] = 0.0
    knots_u[1] = 0.0
    knots_u[2] = 1.0
    knots_u[3] = 1.0

    knots_v = KM.Vector(4)
    knots_v[0] = 0.0
    knots_v[1] = 0.0
    knots_v[2] = 1.0
    knots_v[3] = 1.0

    surface = KM.NurbsSurfaceGeometry3D(
        nodes, 2, 2, knots_u, knots_v)

    return surface

######################################################################################################

if __name__ == "__main__":
    surface = solve_cantilever(create_geometry)

    for node in surface:
        KratosUnittest.TestCase().assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)
        KratosUnittest.TestCase().assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Y), 0)

    # KratosUnittest.TestCase().assertAlmostEqual(surface[0].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[3].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[6].Z, 0.0)

    # KratosUnittest.TestCase().assertAlmostEqual(surface[1].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[4].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[7].Z, 0.0)

    # KratosUnittest.TestCase().assertAlmostEqual(surface[2].Z, -0.223839729168301)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[5].Z, -0.223958628927084)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[8].Z, -0.223839729168301)
