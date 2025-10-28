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


def solve_cantilever(create_geometry):
    model = KM.Model()
    model_part = model.CreateModelPart('Model')

    model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KM.REACTION)
    model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

    # create property for shell elements

    shell_properties = model_part.GetProperties()[1]
    shell_properties.SetValue(KM.THICKNESS, 0.1)
    shell_properties.SetValue(KM.YOUNG_MODULUS, 200000000)
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
0.000000005,
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

    delta = 1e-6
    print(model_part)

    for element in model_part.Elements:
        lhs = KM.Matrix()
        ref_rhs = KM.Vector()

        element.Initialize(model_part.ProcessInfo)
        element.InitializeNonLinearIteration(model_part.ProcessInfo)

        if element.Id < 100:
            element.CalculateLocalSystem(lhs, ref_rhs, model_part.ProcessInfo)

            print("ref_rhs")
            print(ref_rhs)



            # adjoint calculation
            adjoint_element: KM.Element = model_part.GetElement(element.Id + 100)
            adjoint_sensitivities = KM.Matrix()
            adjoint_element.Initialize(model_part.ProcessInfo)
            adjoint_element.InitializeNonLinearIteration(model_part.ProcessInfo)
            adjoint_element.CalculateLeftHandSide(adjoint_sensitivities, model_part.ProcessInfo)

            perturbed_rhs = KM.Vector()
            # for i, node in enumerate(element.GetGeometry()):

                # # print(node)

                # node.X += delta

                # # node.SetSolutionStepValue(KM.DISPLACEMENT_X, node.GetSolutionStepValue(KM.DISPLACEMENT_X) + delta)

                # # print(node.GetSolutionStepValue(KM.DISPLACEMENT_X))

                # element.InitializeNonLinearIteration(model_part.ProcessInfo)
                # element.CalculateLocalSystem(lhs, perturbed_rhs, model_part.ProcessInfo)
                # sensitivities = (perturbed_rhs - ref_rhs) / delta
                # # node.SetSolutionStepValue(KM.DISPLACEMENT_X, node.GetSolutionStepValue(KM.DISPLACEMENT_X) - delta)

                # node.X -= delta

                # print("perturbed")
                # print(perturbed_rhs)

                #         # print(sensitivities)

                # for i in range(27):
                #     print(sensitivities[i], adjoint_sensitivities[i * 3, i])

            for node in model["ACTIVE_MP"].Nodes:
                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA) + delta)
                element.InitializeNonLinearIteration(model_part.ProcessInfo)
                element.CalculateLocalSystem(lhs, perturbed_rhs, model_part.ProcessInfo)
                sensitivities = (perturbed_rhs - ref_rhs) / delta

                # rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_ALPHA)); 27
                # rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_BETA)); 28
                # rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_GAMMA)); 29
                # rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_KAPPA_1)); 30
                # rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_KAPPA_2)); 31
                # rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_KAPPA_12)); 32

                for i in range(33):
                    print(sensitivities[i], adjoint_sensitivities[27, i])

                node.SetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA, node.GetSolutionStepValue(IGA.ACTIVE_SHELL_ALPHA) - delta)



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
    knots_u[2] = 5.0
    knots_u[3] = 5.0

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
