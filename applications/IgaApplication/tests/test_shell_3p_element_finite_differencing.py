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
        element_id += 1

    # assign active shell element inputs
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

            # adjoint calculation
            adjoint_element: KM.Element = model_part.GetElement(element.Id + 100)
            adjoint_LHS = KM.Matrix()
            adjoint_element.Initialize(model_part.ProcessInfo)
            adjoint_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            adjoint_element.CalculateLeftHandSide(adjoint_LHS, model_part.ProcessInfo)

            adjoint_rhs = KM.Vector()
            adjoint_element.CalculateRightHandSide(adjoint_rhs, model_part.ProcessInfo)


            # Finite differencing loop over node perturbations
            perturbed_primal_rhs = KM.Vector()

            target_node_id = 9

            if not model_part.HasNode(target_node_id):
                raise RuntimeError(f"Model part does not contain node {target_node_id}")

            target_node = model_part.GetNode(target_node_id)
            
            #  node.X  #######################################
            target_idx = 24  # column index associated with node 9 displacement X

            target_node.X += delta

            element.InitializeNonLinearIteration(model_part.ProcessInfo)
            element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)

            fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta

            target_node.X -= delta

            for eq_idx in range(33):
                print(adjoint_LHS[eq_idx, target_idx], -fd_sensitivities[eq_idx])
                KratosUnittest.TestCase().assertAlmostEqual(adjoint_LHS[eq_idx, target_idx], -fd_sensitivities[eq_idx],places=3)


            #  node.Y  #######################################
            target_idx = 25  # column index associated with node 9 displacement Y

            target_node.Y += delta

            element.InitializeNonLinearIteration(model_part.ProcessInfo)
            element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)

            fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta

            target_node.Y -= delta

            for eq_idx in range(26):
                KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)


            #  node.Z  #######################################
            target_idx = 26  # column index associated with node 9 displacement Z

            target_node.Z += delta

            element.InitializeNonLinearIteration(model_part.ProcessInfo)
            element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)

            fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta

            target_node.Z -= delta

            for eq_idx in range(26):
                KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=4)

            #  ALPHA  #######################################
            print("++++++++++++++Starting actuation variable sensitivity FD checks...")
            target_idx = 27  # column index associated with Global Point - actuation variable Alpha
            delta = 1e-8
            # There is only one node (ACTIVE_MP) in this background geometry -  the loop is redundant here but kept for generality
            print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            for node in model["ACTIVE_MP"].Nodes:
                print("perturbing node: ", node,"  - node id:", node.Id)
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
                fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA) - delta)

                print("perturbed ActiveGlobalNode - Variable ALPHA")
                print("LHS_adjoint versus Fdiff_RHS_primal :  -------------------")
                for eq_idx in range(33):
                    print(adjoint_LHS[eq_idx, target_idx],fd_sensitivities[eq_idx])
                    
                for eq_idx in range(33):
                    KratosUnittest.TestCase().assertAlmostEqual(adjoint_LHS[eq_idx, target_idx], fd_sensitivities[eq_idx], places=3)

                # Note: A recurring mismatch is observed here, e.g.
                #   1253.858024691359 != 2218.364249984006
                # for eq_idx = 27.
                #
                # Root cause is not yet identified and requires further investigation.
                # See Leonhard Rieder, "Implementation and inverse analysis of active shells for soft robotics",
                # Master’s thesis report, pp. 71–73.
                #
                # The report includes an example LHS entry calculation (pp. 33–35) indicating that 1253.85
                # corresponds to the contribution of the first quadrature point to the element’s local stiffness
                # matrix (LHS). This suggests that quadrature-point contributions may not be accumulated
                # correctly in the ActiveShell element.


            #raise Exception("Alpha check complete: stop here.")

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
                    #print(fd_sensitivities[eq_idx], - adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)

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
                    #print(fd_sensitivities[eq_idx], - adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)

            #  KAPPA_1  #######################################
            target_idx = 30  # column index associated with Global Point - actuation variable Kappa_1

            # There is only one node per background geometry -  the loop is redundant here but kept for generality
            print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            for node in model["ACTIVE_MP"].Nodes:
                print("perturbing node: ", node,"  - node id:", node.Id)
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
                fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_1) - delta)

                print("perturbed ActiveGlobalNode - Variable KAPPA_1")
                print("fd_sensitivities:",fd_sensitivities)
                print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
                for eq_idx in range(33):
                    #print(fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)

            #  KAPPA_2  #######################################
            target_idx = 31  # column index associated with Global Point - actuation variable Kappa_2

            # There is only one node per background geometry -  the loop is redundant here but kept for generality
            print("number of nodes in ACTIVE_MP:", model["ACTIVE_MP"].NumberOfNodes())
            for node in model["ACTIVE_MP"].Nodes:
                print("perturbing node: ", node,"  - node id:", node.Id)
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(primal_LHS, perturbed_primal_rhs, model_part.ProcessInfo)
                fd_sensitivities = (perturbed_primal_rhs - primal_rhs) / delta
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_KAPPA_2) - delta)

                print("perturbed ActiveGlobalNode - Variable KAPPA_2")
                print("fd_sensitivities:",fd_sensitivities)
                print("Fdiff_RHS_primal versus LHS_adjoint:  -------------------")
                for eq_idx in range(33):
                    # print(fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)

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
                    # print(fd_sensitivities[eq_idx], adjoint_LHS[eq_idx, target_idx])
                    KratosUnittest.TestCase().assertAlmostEqual(fd_sensitivities[eq_idx], -adjoint_LHS[eq_idx, target_idx], places=3)

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
