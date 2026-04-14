// test_condition.cpp

#include "test_condition.h"

namespace Kratos
{

    // =======================
    // EquationIdVector 2D2N
    // =======================

    template <>
    void TestCondition<2, 2>::EquationIdVector(
        EquationIdVectorType &rResult,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_WATCH("this");

        const auto &r_geom = GetGeometry();

        const unsigned int num_nodes = 2;
        const unsigned int block_size = 4;

        rResult.resize(num_nodes * block_size, false);

        unsigned int k = 0;

        for (unsigned int i = 0; i < num_nodes; ++i)
        {
            const unsigned int den_pos = r_geom[i].GetDofPosition(DENSITY);
            const unsigned int mom_pos = r_geom[i].GetDofPosition(MOMENTUM);
            const unsigned int enr_pos = r_geom[i].GetDofPosition(TOTAL_ENERGY);

            rResult[k++] = r_geom[i].GetDof(DENSITY, den_pos).EquationId();
            rResult[k++] = r_geom[i].GetDof(MOMENTUM_X, mom_pos).EquationId();
            rResult[k++] = r_geom[i].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
            rResult[k++] = r_geom[i].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
        }
    }

    // =======================
    // EquationIdVector 3D3N
    // =======================

    template <>
    void TestCondition<3, 3>::EquationIdVector(
        EquationIdVectorType &rResult,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        const auto &r_geom = GetGeometry();

        const unsigned int num_nodes = 3;
        const unsigned int block_size = 5;

        rResult.resize(num_nodes * block_size, false);

        unsigned int k = 0;

        for (unsigned int i = 0; i < num_nodes; ++i)
        {
            const unsigned int den_pos = r_geom[i].GetDofPosition(DENSITY);
            const unsigned int mom_pos = r_geom[i].GetDofPosition(MOMENTUM);
            const unsigned int enr_pos = r_geom[i].GetDofPosition(TOTAL_ENERGY);

            rResult[k++] = r_geom[i].GetDof(DENSITY, den_pos).EquationId();
            rResult[k++] = r_geom[i].GetDof(MOMENTUM_X, mom_pos).EquationId();
            rResult[k++] = r_geom[i].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
            rResult[k++] = r_geom[i].GetDof(MOMENTUM_Z, mom_pos + 2).EquationId();
            rResult[k++] = r_geom[i].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
        }
    }

    // =======================
    // GetDofList 2D2N
    // =======================

    template <>
    void TestCondition<2, 2>::GetDofList(
        DofsVectorType &rDofs,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        const auto &r_geom = GetGeometry();

        const unsigned int num_nodes = 2;
        const unsigned int block_size = 4;

        rDofs.resize(num_nodes * block_size);

        unsigned int k = 0;

        for (unsigned int i = 0; i < num_nodes; ++i)
        {
            const unsigned int den_pos = r_geom[i].GetDofPosition(DENSITY);
            const unsigned int mom_pos = r_geom[i].GetDofPosition(MOMENTUM);
            const unsigned int enr_pos = r_geom[i].GetDofPosition(TOTAL_ENERGY);

            rDofs[k++] = r_geom[i].pGetDof(DENSITY, den_pos);
            rDofs[k++] = r_geom[i].pGetDof(MOMENTUM_X, mom_pos);
            rDofs[k++] = r_geom[i].pGetDof(MOMENTUM_Y, mom_pos + 1);
            rDofs[k++] = r_geom[i].pGetDof(TOTAL_ENERGY, enr_pos);
        }
    }

    // =======================
    // GetDofList 3D3N
    // =======================

    template <>
    void TestCondition<3, 3>::GetDofList(
        DofsVectorType &rDofs,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        const auto &r_geom = GetGeometry();

        const unsigned int num_nodes = 3;
        const unsigned int block_size = 5;

        rDofs.resize(num_nodes * block_size);

        unsigned int k = 0;

        for (unsigned int i = 0; i < num_nodes; ++i)
        {
            const unsigned int den_pos = r_geom[i].GetDofPosition(DENSITY);
            const unsigned int mom_pos = r_geom[i].GetDofPosition(MOMENTUM);
            const unsigned int enr_pos = r_geom[i].GetDofPosition(TOTAL_ENERGY);

            rDofs[k++] = r_geom[i].pGetDof(DENSITY, den_pos);
            rDofs[k++] = r_geom[i].pGetDof(MOMENTUM_X, mom_pos);
            rDofs[k++] = r_geom[i].pGetDof(MOMENTUM_Y, mom_pos + 1);
            rDofs[k++] = r_geom[i].pGetDof(MOMENTUM_Z, mom_pos + 2);
            rDofs[k++] = r_geom[i].pGetDof(TOTAL_ENERGY, enr_pos);
        }
    }

    // =======================
    // CalculateRightHandSide
    // =======================

    template <unsigned int TDim, unsigned int TNumNodes>
    void TestCondition<TDim, TNumNodes>::CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
    {
        GeometryType &r_geom = GetGeometry();

        const unsigned int dim = TDim;
        const unsigned int num_nodes = TNumNodes;
        const unsigned int block_size = 4;

        // ---- Initialize RHS
        rRightHandSideVector.resize(num_nodes * block_size, false);
        noalias(rRightHandSideVector) = ZeroVector(num_nodes * block_size);

        // Integration
        const auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto &r_integration_points = r_geom.IntegrationPoints(integration_method);
        const Matrix &Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

        // Viscosity
        const double mu = GetProperties()[DYNAMIC_VISCOSITY];

        // ---- Normal
        array_1d<double, 2> normal;
        const double dx = r_geom[1].X() - r_geom[0].X();
        const double dy = r_geom[1].Y() - r_geom[0].Y();

        normal[0] = -dy;
        normal[1] = dx;

        const double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
        if (length > 1e-12)
            normal /= length;

        // Parent element
        if (!r_geom[0].Has(NEIGHBOUR_ELEMENTS))
        {
            KRATOS_WARNING("TestCondition") << "" << std::endl;
            return;
        }

        const auto &r_neigh_elems = r_geom[0].GetValue(NEIGHBOUR_ELEMENTS);

        // Robust search for the element that contains BOTH nodes of the edge
        const Element *p_elem = nullptr;

        for (unsigned int k = 0; k < r_neigh_elems.size(); ++k)
        {
            const Element &r_candidate = r_neigh_elems[k];
            const auto &geom = r_candidate.GetGeometry();

            bool has_node0 = false;
            bool has_node1 = false;

            for (unsigned int i = 0; i < geom.size(); ++i)
            {
                if (geom[i].Id() == r_geom[0].Id())
                    has_node0 = true;
                if (geom[i].Id() == r_geom[1].Id())
                    has_node1 = true;
            }

            if (has_node0 && has_node1)
            {
                p_elem = &r_candidate;
                break;
            }
        }

        if (!p_elem)
        {
            KRATOS_WARNING("TestCondition") << "No parent element found for edge ("
                                            << r_geom[0].Id() << ", "
                                            << r_geom[1].Id() << ")" << std::endl;
            return;
        }

        const Element &r_elem = *p_elem;

        const auto elem_method = GeometryData::IntegrationMethod::GI_GAUSS_2;

        // if (this->Id() == 100004)
        // {
        //     KRATOS_WATCH(r_elem.GetGeometry()[0].FastGetSolutionStepValue(DENSITY));
        //     KRATOS_WATCH(r_elem.GetGeometry()[0].FastGetSolutionStepValue(MOMENTUM));
        //     KRATOS_WATCH(r_elem.GetGeometry()[0].FastGetSolutionStepValue(TOTAL_ENERGY));

        //     KRATOS_WATCH(r_elem.GetGeometry()[1].FastGetSolutionStepValue(DENSITY));
        //     KRATOS_WATCH(r_elem.GetGeometry()[1].FastGetSolutionStepValue(MOMENTUM));
        //     KRATOS_WATCH(r_elem.GetGeometry()[1].FastGetSolutionStepValue(TOTAL_ENERGY));

        //     KRATOS_WATCH(r_elem.GetGeometry()[2].FastGetSolutionStepValue(DENSITY));
        //     KRATOS_WATCH(r_elem.GetGeometry()[2].FastGetSolutionStepValue(MOMENTUM));
        //     KRATOS_WATCH(r_elem.GetGeometry()[2].FastGetSolutionStepValue(TOTAL_ENERGY));
        // }

        // MANUAL DN_DX RECONSTRUCTION

        const auto &r_geom_elem = r_elem.GetGeometry();

        // Node coordinates of the triangle
        const array_1d<double, 3> &node_0_coord = r_geom_elem[0].Coordinates();
        const array_1d<double, 3> &node_1_coord = r_geom_elem[1].Coordinates();
        const array_1d<double, 3> &node_2_coord = r_geom_elem[2].Coordinates();

        // Local derivatives of linear triangle
        const double dN0_dxi = -1.0;
        const double dN0_deta = -1.0;
        const double dN1_dxi = 1.0;
        const double dN1_deta = 0.0;
        const double dN2_dxi = 0.0;
        const double dN2_deta = 1.0;

        // Build Jacobian manually
        Matrix J(2, 2);
        J.clear();
        J(0, 0) = dN0_dxi * node_0_coord[0] + dN1_dxi * node_1_coord[0] + dN2_dxi * node_2_coord[0];
        J(0, 1) = dN0_deta * node_0_coord[0] + dN1_deta * node_1_coord[0] + dN2_deta * node_2_coord[0];
        J(1, 0) = dN0_dxi * node_0_coord[1] + dN1_dxi * node_1_coord[1] + dN2_dxi * node_2_coord[1];
        J(1, 1) = dN0_deta * node_0_coord[1] + dN1_deta * node_1_coord[1] + dN2_deta * node_2_coord[1];

        // Invert Jacobian
        Matrix inv_J(2, 2);
        double detJ_manual = 0.0;
        MathUtils<double>::InvertMatrix(J, inv_J, detJ_manual);

        // Compute DN_DX manually
        Matrix DN_DX_manual(3, 2);
        DN_DX_manual.clear();

        // node 0
        Vector dN_De(2);
        Vector dN_DX(2);

        dN_De[0] = dN0_dxi;
        dN_De[1] = dN0_deta;

        noalias(dN_DX) = prod(inv_J, dN_De);

        DN_DX_manual(0, 0) = dN_DX[0];
        DN_DX_manual(0, 1) = dN_DX[1];

        // node 1
        dN_De[0] = dN1_dxi;
        dN_De[1] = dN1_deta;

        noalias(dN_DX) = prod(inv_J, dN_De);

        DN_DX_manual(1, 0) = dN_DX[0];
        DN_DX_manual(1, 1) = dN_DX[1];

        // node 2
        dN_De[0] = dN2_dxi;
        dN_De[1] = dN2_deta;

        noalias(dN_DX) = prod(inv_J, dN_De);

        DN_DX_manual(2, 0) = dN_DX[0];
        DN_DX_manual(2, 1) = dN_DX[1];

        // Velocity gradient (for linear triangle)

        Matrix grad_u(2, 2);
        grad_u.clear();

        const auto &v0 = r_elem.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const auto &v1 = r_elem.GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
        const auto &v2 = r_elem.GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

        // du/dx
        grad_u(0, 0) =
            v0[0] * DN_DX_manual(0, 0) + v1[0] * DN_DX_manual(1, 0) + v2[0] * DN_DX_manual(2, 0);
        // du/dy
        grad_u(0, 1) =
            v0[0] * DN_DX_manual(0, 1) + v1[0] * DN_DX_manual(1, 1) + v2[0] * DN_DX_manual(2, 1);
        // dv/dx
        grad_u(1, 0) =
            v0[1] * DN_DX_manual(0, 0) + v1[1] * DN_DX_manual(1, 0) + v2[1] * DN_DX_manual(2, 0);
        // dv/dy
        grad_u(1, 1) =
            v0[1] * DN_DX_manual(0, 1) + v1[1] * DN_DX_manual(1, 1) + v2[1] * DN_DX_manual(2, 1);

        // Stress tensor τ (explicit 2D form)

        Matrix tau(2, 2);
        tau.clear();

        // divergence
        const double div_u = grad_u(0, 0) + grad_u(1, 1);

        // tau_xx
        tau(0, 0) = 2.0 * mu * grad_u(0, 0) - (2.0 / 3.0) * mu * div_u;

        // tau_xy
        tau(0, 1) = mu * (grad_u(0, 1) + grad_u(1, 0));

        // tau_yx
        tau(1, 0) = tau(0, 1);

        // tau_yy
        tau(1, 1) = 2.0 * mu * grad_u(1, 1) - (2.0 / 3.0) * mu * div_u;

        // Penalty
        const double gamma = 50.0;

        // Boundary Jacobian
        const double dx_line = r_geom[1].X() - r_geom[0].X();
        const double dy_line = r_geom[1].Y() - r_geom[0].Y();
        const double edge_length = std::sqrt(dx_line * dx_line + dy_line * dy_line);
        const double detJ_line = edge_length * 0.5;

        // if (this->Id() == 100004)
        // {
        //     KRATOS_WATCH(mu);
        //     KRATOS_WATCH(length);
        //     KRATOS_WATCH(detJ_manual);
        //     KRATOS_WATCH(normal);
        //     KRATOS_WATCH(v0);
        //     KRATOS_WATCH(v1);
        //     KRATOS_WATCH(v2);
        // }
        // ---- Gauss loop
        for (unsigned int g = 0; g < r_integration_points.size(); ++g)
        {
            const double weight = r_integration_points[g].Weight();
            const double coeff = weight * detJ_line;

            // τ · n
            array_1d<double, 2> traction;
            noalias(traction) = prod(tau, normal);

            // Build nodal momentum matrix M (num_nodes x dim)
            Matrix M(num_nodes, dim);

            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                const auto &mom = r_geom[i].FastGetSolutionStepValue(MOMENTUM);
                M(i, 0) = mom[0];
                M(i, 1) = mom[1];
            }

            // Shape function row (1 x num_nodes)
            Vector N_row(num_nodes);
            for (unsigned int i = 0; i < num_nodes; ++i)
                N_row[i] = Ncontainer(g, i);

            //  Compute m_gp = N * M -> (1x2)
            Vector m_gp_vec = prod(N_row, M);

            //  Convert to array_1d
            array_1d<double, 2> m_gp;
            m_gp[0] = m_gp_vec[0];
            m_gp[1] = m_gp_vec[1];

            //  m · n
            const double m_dot_n = inner_prod(m_gp, normal);

            //  Build Nw^T (8x2)
            Matrix NwT(num_nodes * block_size, 2, 0.0);

            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                const double Ni = Ncontainer(g, i);
                const unsigned int base = i * block_size;

                NwT(base + 1, 0) = Ni;
                NwT(base + 2, 1) = Ni;
            }

            // Convert vectors to Kratos Vector
            Vector traction_vec(2);
            traction_vec[0] = traction[0];
            traction_vec[1] = traction[1];

            Vector penalty_vec(2);
            penalty_vec[0] = gamma * m_dot_n * normal[0];
            penalty_vec[1] = gamma * m_dot_n * normal[1];

            Vector R_total = -prod(NwT, traction_vec) + prod(NwT, penalty_vec);

            noalias(rRightHandSideVector) -= R_total * coeff;

            // if (this->Id() == 100004)
            //{
            //     KRATOS_WATCH(m_gp);
            //     KRATOS_WATCH(m_dot_n);
            //     KRATOS_WATCH(traction);
            //     KRATOS_WATCH(coeff);
            // }
            //
            // if (this->Id() == 100004)
            //{
            //     KRATOS_WATCH(R_total.size());
            //     KRATOS_WATCH(rRightHandSideVector);
            // }
        }
    }

    // =======================
    // Explicit instantiation
    // =======================

    template class TestCondition<2, 2>;
    template class TestCondition<3, 3>;

} // namespace Kratos