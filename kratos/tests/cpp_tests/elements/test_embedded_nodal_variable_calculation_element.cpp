//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kratos default authors
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "elements/embedded_nodal_variable_calculation_element_simplex.h"
#include "elements/levelset_convection_element_simplex_algebraic_stabilization.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
namespace Testing
{
    KRATOS_TEST_CASE_IN_SUITE(EmbeddedNodalVariableCalculationElementVectorialRHS2D, KratosCoreFastSuite)
    {
        // The penalty term of the vector-valued specialization must couple node i
        // with the values of node j: rhs_i = f_i - sum_j K_ij g_j. The buggy version
        // used g_i in place of g_j inside the j-loop, which is only invisible when
        // both nodes carry the same g: hence the deliberately distinct nodal values.
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");
        model_part.SetBufferSize(1);

        model_part.AddNodalSolutionStepVariable(NODAL_MAUX);
        model_part.AddNodalSolutionStepVariable(NODAL_VAUX);
        model_part.GetProcessInfo().SetValue(GRADIENT_PENALTY_COEFFICIENT, 0.3);

        using ElementType = EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>;
        auto p_node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = model_part.CreateNewNode(2, 2.0, 0.0, 0.0);
        ElementType::NodesArrayType nodes;
        nodes.push_back(p_node_1);
        nodes.push_back(p_node_2);

        ElementType element(1, nodes);
        auto& r_geom = element.GetGeometry();

        // Element-level interpolation data (stored on the first node, non-historical)
        const double dist = 0.25;
        element.SetValue(DISTANCE, dist);
        element.SetValue(NODAL_VAUX, array_1d<double,3>{1.5, -2.0, 0.5});

        // Distinct per-node historical values: the only configuration where a wrong j-index shows up
        noalias(r_geom[0].FastGetSolutionStepValue(NODAL_VAUX)) = array_1d<double,3>{1.0, 0.0, 0.0};
        noalias(r_geom[1].FastGetSolutionStepValue(NODAL_VAUX)) = array_1d<double,3>{0.0, 3.0, 0.0};

        Vector rhs = ZeroVector(6);
        element.CalculateRightHandSide(rhs, model_part.GetProcessInfo());

        // Independent reference, replicating the intended formula directly from
        // the values set above (element-level values live in the geometry's own
        // container, so they are not re-read from the nodes here).
        const array_1d<double,3> rData{1.5, -2.0, 0.5};
        const array_1d<double, 2> rN{{1.0 - dist, dist}};
        const double penalty = 0.3;
        const std::array<double, 2> aux_penalty{{penalty, -penalty}};

        array_1d<double,3> g0, g1;
        noalias(g0) = r_geom[0].FastGetSolutionStepValue(NODAL_VAUX);
        noalias(g1) = r_geom[1].FastGetSolutionStepValue(NODAL_VAUX);

        for (unsigned int i = 0; i < 2; ++i) {
            for (unsigned int k = 0; k < 3; ++k) {
                double expected = rN[i] * rData(k);
                for (unsigned int j = 0; j < 2; ++j) {
                    const array_1d<double,3>& g_j = (j == 0) ? g0 : g1;
                    expected -= (rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j]) * g_j[k];
                }
                KRATOS_EXPECT_NEAR(expected, rhs(i * 3 + k), 1e-12);
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(LevelSetConvectionAlgebraicStabilizationLHS2D, KratosCoreFastSuite)
    {
        // S_vector_LHS accumulates a per-Gauss-point contribution (its companion
        // RHS stabilization S_vector integrates with += over the same loop). The
        // buggy implementation used a plain '=', leaving only the last Gauss
        // point's contribution in the LHS of the theta-scheme.
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");
        model_part.SetBufferSize(3);

        model_part.AddNodalSolutionStepVariable(VELOCITY);
        model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
        model_part.AddNodalSolutionStepVariable(DISTANCE);

        model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.1);
        model_part.GetProcessInfo().SetValue(TIME_INTEGRATION_THETA, 0.7);

        ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
        p_conv_diff_settings->SetUnknownVariable(DISTANCE);
        p_conv_diff_settings->SetConvectionVariable(VELOCITY);
        p_conv_diff_settings->SetMeshVelocityVariable(MESH_VELOCITY);
        // The element reads grad_phi_mean from the NON-historical values of the
        // gradient variable (r_node.GetValue(r_grad_var)); DISTANCE_GRADIENT is
        // used here so that the test is self-contained.
        p_conv_diff_settings->SetGradientVariable(DISTANCE_GRADIENT);
        model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);

        auto p_elem_prop = model_part.CreateNewProperties(0);

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
        auto p_element = model_part.CreateNewElement(
            "LevelSetConvectionElementSimplexAlgebraicStabilization2D3N", 1, elem_nodes, p_elem_prop);

        auto& r_geom = p_element->GetGeometry();
        Matrix vel(3,2);
        vel(0,0) = 0.1;  vel(0,1) = 0.2;
        vel(1,0) = 0.3;  vel(1,1) = 0.1;
        vel(2,0) = -0.2; vel(2,1) = 0.4;
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int k = 0; k < 2; ++k) {
                r_geom[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel(i,k);
                r_geom[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.5 * vel(i,k);
            }
            r_geom[i].FastGetSolutionStepValue(DISTANCE)    = -0.5 + 0.75 * i;
            r_geom[i].FastGetSolutionStepValue(DISTANCE, 1) = -0.3 + 0.4  * i;
            // Non-historical gradient variable values read by the element
            r_geom[i].SetValue(DISTANCE_GRADIENT, array_1d<double,3>{0.2 - 0.1 * i, 0.5 + 0.3 * i, 0.0});
        }

        // Activate the high-order (stabilization) branch
        p_element->SetValue(LIMITER_COEFFICIENT, 1.0);

        const ProcessInfo& r_process_info = model_part.GetProcessInfo();
        p_element->Initialize(r_process_info);

        Matrix lhs = ZeroMatrix(3,3);
        Vector rhs = ZeroVector(3);
        p_element->CalculateLocalSystem(lhs, rhs, r_process_info);

        // ------------------------------------------------------------------
        // Independent replication of the (correct) element algebra
        // ------------------------------------------------------------------
        BoundedMatrix<double, 3, 2> DN_DX;
        array_1d<double, 3> N;
        double Volume = 0.0;
        GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, Volume);

        array_1d<double, 3> phi, phi_old;
        array_1d<array_1d<double, 3>, 3> v, vold;
        array_1d<array_1d<double, 3>, 3> grad_phi_nodes;
        array_1d<array_1d<double, 3>, 3> x_node;
        array_1d<double, 3> x_mean_tmp = ZeroVector(3);
        double phi_mean = 0.0;
        for (unsigned int i = 0; i < 3; ++i) {
            phi[i]      = r_geom[i].FastGetSolutionStepValue(DISTANCE);
            phi_old[i]  = r_geom[i].FastGetSolutionStepValue(DISTANCE, 1);
            v[i]        = r_geom[i].FastGetSolutionStepValue(VELOCITY)    - r_geom[i].FastGetSolutionStepValue(MESH_VELOCITY);
            vold[i]     = r_geom[i].FastGetSolutionStepValue(VELOCITY, 1) - r_geom[i].FastGetSolutionStepValue(MESH_VELOCITY, 1);
            grad_phi_nodes[i] = r_geom[i].GetValue(DISTANCE_GRADIENT); // non-historical, as in the element
            x_node[i]   = r_geom[i].Coordinates();
            x_mean_tmp += r_geom[i].Coordinates();
            phi_mean += r_geom[i].FastGetSolutionStepValue(DISTANCE);
        }
        const double aux_weight = 1.0 / 3.0;
        phi_mean *= aux_weight;
        array_1d<double, 2> x_mean = ZeroVector(2);
        for (unsigned int k = 0; k < 2; ++k) x_mean[k] = aux_weight * x_mean_tmp[k];

        // Gauss shape functions (same data as GetShapeFunctionsOnGauss for 2D3N)
        BoundedMatrix<double, 3, 3> Ncontainer;
        const double one_sixt = 1.0/6.0, two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt;  Ncontainer(0,1) = one_sixt;  Ncontainer(0,2) = two_third;
        Ncontainer(1,0) = one_sixt;  Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt;  Ncontainer(2,2) = one_sixt;

        const double theta = 0.7;
        const double dt_inv = 10.0;
        const double limiter = 1.0;

        BoundedMatrix<double, 3, 3> K_matrix = ZeroMatrix(3,3);
        BoundedMatrix<double, 3, 3> Mc_matrix = ZeroMatrix(3,3);
        BoundedMatrix<double, 3, 3> S_vector_LHS = ZeroMatrix(3,3);
        Vector S_vector = ZeroVector(3);

        for (unsigned int g = 0; g < 3; ++g) {
            array_1d<double, 3> Ng = row(Ncontainer, g);

            array_1d<double, 2> vel_gauss = ZeroVector(2);
            array_1d<double, 2> x_gauss = ZeroVector(2);
            double phi_gauss = 0.0;
            for (unsigned int i = 0; i < 3; ++i) {
                for (unsigned int k = 0; k < 2; ++k) {
                    vel_gauss[k] += Ng[i] * ((1.0 - theta) * vold[i][k] + theta * v[i][k]);
                    x_gauss[k]   += Ng[i] * x_node[i][k];
                }
                phi_gauss += Ng[i] * phi[i];
            }

            array_1d<double, 3> v_dot_grad_N = prod(DN_DX, vel_gauss);
            noalias(Mc_matrix) += outer_prod(Ng, Ng);
            noalias(K_matrix)  += outer_prod(Ng, v_dot_grad_N);

            array_1d<double, 3> grad_phi_mean_tmp = ZeroVector(3);
            for (unsigned int i = 0; i < 3; ++i)
                grad_phi_mean_tmp += grad_phi_nodes[i];

            array_1d<double, 2> grad_phi_mean = ZeroVector(2);
            for (unsigned int k = 0; k < 2; ++k) grad_phi_mean[k] = aux_weight * grad_phi_mean_tmp[k];

            array_1d<double, 3> dx_dot_grad_N = prod(DN_DX, x_gauss - x_mean);

            array_1d<double, 2> x_gauss_minus_mean = x_gauss - x_mean;
            array_1d<double, 2> scaled_grad = ZeroVector(2);
            for (unsigned int k = 0; k < 2; ++k) scaled_grad[k] = (1.0 - theta) * grad_phi_mean[k];

            for (unsigned int i = 0; i < 3; ++i) {
                S_vector[i] += ((1.0 - theta) * (phi_gauss - phi_mean)
                                - inner_prod(scaled_grad, x_gauss_minus_mean)) * Ng[i];
                for (unsigned int j = 0; j < 3; ++j) {
                    S_vector_LHS(i, j) += (aux_weight - Ng[j] + dx_dot_grad_N[j]) * Ng[i];
                }
            }
        }

        BoundedMatrix<double, 3, 3> S_matrix = aux_weight * (IdentityMatrix(3,3) - Mc_matrix);

        double nu_e = 0.0;
        for (unsigned int i = 0; i < 3; ++i)
            for (unsigned int j = 0; j < 3; ++j)
                if (i != j)
                    nu_e = std::max(nu_e, std::max(0.0, K_matrix(i,j)) / std::abs(S_matrix(i,j)));

        BoundedMatrix<double, 3, 3> M_matrix, L_matrix;
        noalias(M_matrix) = dt_inv * ((1.0 - limiter) * IdentityMatrix(3,3) + limiter * Mc_matrix);
        noalias(L_matrix) = K_matrix + (1.0 - limiter) * nu_e * S_matrix;

        Matrix lhs_reference = ZeroMatrix(3,3);
        Vector rhs_reference = ZeroVector(3);
        noalias(lhs_reference) = M_matrix + theta * (L_matrix + S_vector_LHS);
        noalias(rhs_reference) = prod(M_matrix - (1.0 - theta) * L_matrix, phi_old) - limiter * nu_e * S_vector;
        noalias(rhs_reference) -= prod(lhs_reference, phi);

        const double gauss_pt_weight = Volume * aux_weight;
        lhs_reference *= gauss_pt_weight;
        rhs_reference *= gauss_pt_weight;

        KRATOS_EXPECT_MATRIX_NEAR(lhs, lhs_reference, 1e-10);
        KRATOS_EXPECT_VECTOR_NEAR(rhs, rhs_reference, 1e-10);
    }

} // namespace Testing.
} // namespace Kratos.
