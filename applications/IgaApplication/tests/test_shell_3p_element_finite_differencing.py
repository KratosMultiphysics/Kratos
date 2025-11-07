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
import math


def _report_matrix_diff(reference_matrix, comparison_matrix, diff_threshold=1e-4,
                        reference_label="primal_LHS", comparison_label="adjoint_LHS"):
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


def _report_matrix_diff_first_entries(reference_matrix, comparison_matrix, block_size=27,
                                      diff_threshold=1e-4, reference_label="primal_LHS",
                                      comparison_label="shell3p_LHS"):
    """Print differences for the top-left block of two matrices."""
    if reference_matrix.Size1() < block_size or reference_matrix.Size2() < block_size:
        raise RuntimeError(f"{reference_label} is smaller than the requested block size {block_size}")
    if comparison_matrix.Size1() < block_size or comparison_matrix.Size2() < block_size:
        raise RuntimeError(f"{comparison_label} is smaller than the requested block size {block_size}")

    large_differences = []
    for row in range(block_size):
        for col in range(block_size):
            diff = abs(reference_matrix[row, col] - comparison_matrix[row, col])
            if diff > diff_threshold:
                large_differences.append((row, col, reference_matrix[row, col], comparison_matrix[row, col], diff))

    if large_differences:
        print(f"Entries with |{reference_label}[0:{block_size}] - {comparison_label}[0:{block_size}]| > {diff_threshold}:")
        for row, col, ref_val, cmp_val, diff in large_differences:
            print(f"  row {row}, col {col}: {reference_label}={ref_val}, {comparison_label}={cmp_val}, diff={diff}")
    else:
        print(f"{reference_label.upper()} vs {comparison_label.upper()} (first {block_size}) - No entries differ by more than {diff_threshold}")


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


def _report_vector_diff_first_entries(reference_vector, comparison_vector, vector_size=27,
                                      diff_threshold=1e-4, reference_label="primal_rhs",
                                      comparison_label="shell3p_rhs"):
    """Print differences for the leading entries of two vectors."""
    if reference_vector.Size() < vector_size:
        raise RuntimeError(f"{reference_label} has fewer than {vector_size} entries")
    if comparison_vector.Size() < vector_size:
        raise RuntimeError(f"{comparison_label} has fewer than {vector_size} entries")

    large_differences = []
    for idx in range(vector_size):
        diff = abs(reference_vector[idx] - comparison_vector[idx])
        if diff > diff_threshold:
            large_differences.append((idx, reference_vector[idx], comparison_vector[idx], diff))

    if large_differences:
        print(f"Entries with |{reference_label}[0:{vector_size}] - {comparison_label}[0:{vector_size}]| > {diff_threshold}:")
        for idx, ref_val, cmp_val, diff in large_differences:
            print(f"  index {idx}: {reference_label}={ref_val}, {comparison_label}={cmp_val}, diff={diff}")
    else:
        print(f"{reference_label.upper()} vs {comparison_label.upper()} (first {vector_size}) - No entries differ by more than {diff_threshold}")


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
    print(len(quadrature_point_geometries), " quadrature point geometries created.")

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
0.00,
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
    delta = 1e-8
    for element in model_part.Elements:
        primal_LHS = KM.Matrix()
        primal_rhs = KM.Vector()

        element.Initialize(model_part.ProcessInfo)
        element.InitializeNonLinearIteration(model_part.ProcessInfo)

        if element.Id < 100:
            element.CalculateLocalSystem(primal_LHS, primal_rhs, model_part.ProcessInfo)

            # print("primal_rhs\n:", primal_rhs)
            # print("-" * 40)
            # print("primal_LHS\n:", primal_LHS)
            # print("-" * 40)


            # adjoint calculation
            adjoint_element: KM.Element = model_part.GetElement(element.Id + 100)
            adjoint_LHS = KM.Matrix()
            adjoint_element.Initialize(model_part.ProcessInfo)
            adjoint_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            adjoint_element.CalculateLeftHandSide(adjoint_LHS, model_part.ProcessInfo)

            adjoint_rhs = KM.Vector()
            adjoint_element.CalculateRightHandSide(adjoint_rhs, model_part.ProcessInfo)
            
            # _report_matrix_diff(primal_LHS, adjoint_LHS, diff_threshold=1e-4,
            #                     reference_label="primal_LHS", comparison_label="adjoint_LHS")

            # _report_vector_diff(primal_rhs, adjoint_rhs, diff_threshold=1e-4,
            #                     reference_label="primal_rhs", comparison_label="adjoint_rhs")

            #normal Shell3pElement
            Shell3p_element: KM.Element = model_part.GetElement(element.Id + 1000)
            Shell3p_LHS = KM.Matrix()
            shell3p_reference_rhs = KM.Vector()
            Shell3p_element.Initialize(model_part.ProcessInfo)
            Shell3p_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            Shell3p_element.CalculateLocalSystem(Shell3p_LHS, shell3p_reference_rhs, model_part.ProcessInfo)


            # print("size Shell3p_LHS:", Shell3p_LHS.Size1(), Shell3p_LHS.Size2())
            # _report_matrix_diff_first_entries(primal_LHS, Shell3p_LHS, block_size=27, diff_threshold=1e-4,
            #                                   reference_label="primal_LHS", comparison_label="shell3p_LHS")

            # _report_vector_diff_first_entries(primal_rhs, shell3p_reference_rhs, vector_size=27,
            #                                   diff_threshold=1e-4,
            #                                   reference_label="primal_rhs",
            #                                   comparison_label="shell3p_reference_rhs")


            # Finite differencing loop over node perturbations
            perturbed_primal_rhs = KM.Vector()
            shell3p_perturbed_rhs = KM.Vector()

            target_node_id = 9

            if not model_part.HasNode(target_node_id):
                raise RuntimeError(f"Model part does not contain node {target_node_id}")

            target_node = model_part.GetNode(target_node_id)
            
            #  node.X  #######################################
            target_idx = 24  # column index associated with node 9 displacement X

            target_node.X += delta

            element.InitializeNonLinearIteration(model_part.ProcessInfo)
            element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)

            Shell3p_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            Shell3p_element.CalculateLocalSystem(Shell3p_LHS, shell3p_perturbed_rhs, model_part.ProcessInfo)

            fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
            shell3p_fd_sensitivities = (shell3p_perturbed_rhs - shell3p_reference_rhs) / delta

            target_node.X -= delta

            #print(f"perturbed node {target_node_id} Coordinate .X")
            #print("sensitivities comparison: Fdiff_RHS_primal; Fdiff_RHS_Shell3p; LHS_adjoint -------------------")
            for eq_idx in range(33):
                print(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx])
                    
                # print(-adjoint_LHS[eq_idx, 24])
                KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
                #print(fd_sensitivities[eq_idx], shell3p_fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
            # print("-" * 40)

            #  node.Y  #######################################
            target_idx = 25  # column index associated with node 9 displacement Y

            target_node.Y += delta

            element.InitializeNonLinearIteration(model_part.ProcessInfo)
            element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)

            Shell3p_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            Shell3p_element.CalculateLocalSystem(Shell3p_LHS, shell3p_perturbed_rhs, model_part.ProcessInfo)

            fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
            shell3p_fd_sensitivities = (shell3p_perturbed_rhs - shell3p_reference_rhs) / delta

            target_node.Y -= delta

            # print(f"perturbed node {target_node_id} Coordinate .Y")
            # print("sensitivities comparison: Fdiff_RHS_primal; Fdiff_RHS_Shell3p; LHS_adjoint -------------------")
            for eq_idx in range(26):
                KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
            #     print(fd_sensitivities[eq_idx], shell3p_fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
            # print("-" * 40)


            #  node.Z  #######################################
            target_idx = 26  # column index associated with node 9 displacement Z

            target_node.Z += delta

            element.InitializeNonLinearIteration(model_part.ProcessInfo)
            element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)

            Shell3p_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            Shell3p_element.CalculateLocalSystem(Shell3p_LHS, shell3p_perturbed_rhs, model_part.ProcessInfo)

            fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
            shell3p_fd_sensitivities = (shell3p_perturbed_rhs - shell3p_reference_rhs) / delta

            target_node.Z -= delta

            # print(f"perturbed node {target_node_id} Coordinate .Z")
            # print("sensitivities comparison: Fdiff_RHS_primal; Fdiff_RHS_Shell3p; LHS_adjoint -------------------")
            for eq_idx in range(26):
                KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=4)
            #     print(fd_sensitivities[eq_idx], shell3p_fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
            # print("-" * 40)

            #  ALPHA  #######################################
            target_idx = 27  # column index associated with Global Point - actuation variable Alpha

            # There is only one node per background geometry -  the loop is redundant here but kept for generality
            print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            for node in model["ACTIVE_MP"].Nodes:
                print("perturbing node: ", node,"  - node id:", node.Id)
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
                fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA) - delta)

                # print("...............adjoint_LHS:",adjoint_LHS)

                print("perturbed ActiveGlobalNode - Variable ALPHA")
                print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
                for eq_idx in range(33):
                    # print(fd_sensitivities[eq_idx], - adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
            #     print("-" * 40)

            #  BETA  #######################################
            target_idx = 28  # column index associated with Global Point - actuation variable Alpha

            # There is only one node per background geometry -  the loop is redundant here but kept for generality
            print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            for node in model["ACTIVE_MP"].Nodes:
                print("perturbing node: ", node,"  - node id:", node.Id)
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_BETA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_BETA) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
                fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_BETA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_BETA) - delta)

                print("perturbed ActiveGlobalNode - Variable BETA")
                print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
                for eq_idx in range(33):
                    # print(fd_sensitivities[eq_idx], - adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
            #     print("-" * 40)

            #  GAMMA  #######################################
            target_idx = 29  # column index associated with Global Point - actuation variable Alpha

            # There is only one node per background geometry -  the loop is redundant here but kept for generality
            print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            for node in model["ACTIVE_MP"].Nodes:
                print("perturbing node: ", node,"  - node id:", node.Id)
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_GAMMA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_GAMMA) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
                fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_GAMMA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_GAMMA) - delta)

                print("perturbed ActiveGlobalNode - Variable GAMMA")
                print("fd_sensitivities:",fd_sensitivities)
                print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
                for eq_idx in range(33):
                    print(fd_sensitivities[eq_idx], - adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
            #     print("-" * 40)

            # #  KAPPA_1  #######################################
            # target_idx = 30  # column index associated with Global Point - actuation variable Kappa_1

            # # There is only one node per background geometry -  the loop is redundant here but kept for generality
            # print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            # for node in model["ACTIVE_MP"].Nodes:
            #     print("perturbing node: ", node,"  - node id:", node.Id)
            #     node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1) + delta)
            #     element.InitializeNonLinearIteration(model_part.ProcessInfo)
            #     element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
            #     fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
            #     node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1) - delta)

            #     print("perturbed ActiveGlobalNode - Variable KAPPA_1")
            #     print("fd_sensitivities:",fd_sensitivities)
            #     print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
            #     for eq_idx in range(33):
            #         print(fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
            #         KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
            #     # print("-" * 40)

            # #  KAPPA_2  #######################################
            # target_idx = 31  # column index associated with Global Point - actuation variable Kappa_2

            # # There is only one node per background geometry -  the loop is redundant here but kept for generality
            # print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            # for node in model["ACTIVE_MP"].Nodes:
            #     print("perturbing node: ", node,"  - node id:", node.Id)
            #     node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2) + delta)
            #     element.InitializeNonLinearIteration(model_part.ProcessInfo)
            #     element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
            #     fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
            #     node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2) - delta)

            #     print("perturbed ActiveGlobalNode - Variable KAPPA_2")
            #     print("fd_sensitivities:",fd_sensitivities)
            #     print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
            #     for eq_idx in range(33):
            #         print(fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
            #         KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
            #     # print("-" * 40)

            #  KAPPA_12  #######################################
            target_idx = 32  # column index associated with Global Point - actuation variable Kappa_12

            # There is only one node per background geometry -  the loop is redundant here but kept for generality
            print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            for node in model["ACTIVE_MP"].Nodes:
                print("perturbing node: ", node,"  - node id:", node.Id)
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_12, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_12) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
                fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_12, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_12) - delta)

                print("perturbed ActiveGlobalNode - Variable KAPPA_12")
                print("fd_sensitivities:",fd_sensitivities)
                print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
                for eq_idx in range(33):
                    print(fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)
                # print("-" * 40)

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

    #checks whether the displacement values are zero
    for node in surface:
        KratosUnittest.TestCase().assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)
        KratosUnittest.TestCase().assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Y), 0)
        KratosUnittest.TestCase().assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Z), 0)

    # KratosUnittest.TestCase().assertAlmostEqual(surface[0].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[3].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[6].Z, 0.0)

    # KratosUnittest.TestCase().assertAlmostEqual(surface[1].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[4].Z, 0.0)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[7].Z, 0.0)

    # KratosUnittest.TestCase().assertAlmostEqual(surface[2].Z, -0.223839729168301)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[5].Z, -0.223958628927084)
    # KratosUnittest.TestCase().assertAlmostEqual(surface[8].Z, -0.223839729168301)
