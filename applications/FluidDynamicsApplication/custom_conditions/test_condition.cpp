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
        const ProcessInfo &rProcessInfo) const
    {
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
        const ProcessInfo &rProcessInfo) const
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
        const ProcessInfo &rProcessInfo) const
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
        const ProcessInfo &rProcessInfo) const
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
        const ProcessInfo &rProcessInfo)
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
            return;

        const auto &r_geom_elem = p_elem->GetGeometry();

        // ---- DN_DX
        Matrix DN_DX(3, 2);

        const double x0 = r_geom_elem[0].X(), y0 = r_geom_elem[0].Y();
        const double x1 = r_geom_elem[1].X(), y1 = r_geom_elem[1].Y();
        const double x2 = r_geom_elem[2].X(), y2 = r_geom_elem[2].Y();

        const double detJ = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

        DN_DX(0, 0) = (y1 - y2) / detJ;
        DN_DX(0, 1) = (x2 - x1) / detJ;

        DN_DX(1, 0) = (y2 - y0) / detJ;
        DN_DX(1, 1) = (x0 - x2) / detJ;

        DN_DX(2, 0) = (y0 - y1) / detJ;
        DN_DX(2, 1) = (x1 - x0) / detJ;

        // ---- grad_u
        Matrix grad_u(2, 2, 0.0);

        for (unsigned int i = 0; i < 3; ++i)
        {
            const auto &v = r_geom_elem[i].FastGetSolutionStepValue(VELOCITY);

            grad_u(0, 0) += v[0] * DN_DX(i, 0);
            grad_u(0, 1) += v[0] * DN_DX(i, 1);
            grad_u(1, 0) += v[1] * DN_DX(i, 0);
            grad_u(1, 1) += v[1] * DN_DX(i, 1);
        }

        // ---- stress
        Matrix tau(2, 2);

        const double div_u = grad_u(0, 0) + grad_u(1, 1);

        tau(0, 0) = 2.0 * mu * grad_u(0, 0) - (2.0 / 3.0) * mu * div_u;
        tau(0, 1) = mu * (grad_u(0, 1) + grad_u(1, 0));
        tau(1, 0) = tau(0, 1);
        tau(1, 1) = 2.0 * mu * grad_u(1, 1) - (2.0 / 3.0) * mu * div_u;

        // Penalty
        const double gamma = 160.0;

        const double edge_length = std::sqrt(dx * dx + dy * dy);
        const double detJ_line = edge_length * 0.5;

        // ---- Gauss loop
        for (unsigned int g = 0; g < r_integration_points.size(); ++g)
        {
            const double weight = r_integration_points[g].Weight();
            const double coeff = weight * detJ_line;

            array_1d<double, 2> traction;
            noalias(traction) = prod(tau, normal);

            // M
            Matrix M(num_nodes, dim);
            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                const auto &mom = r_geom[i].FastGetSolutionStepValue(MOMENTUM);
                M(i, 0) = mom[0];
                M(i, 1) = mom[1];
            }

            Vector N_row(num_nodes);
            for (unsigned int i = 0; i < num_nodes; ++i)
                N_row[i] = Ncontainer(g, i);

            Vector m_gp = prod(N_row, M);
            const double m_dot_n = m_gp[0] * normal[0] + m_gp[1] * normal[1];

            Matrix NwT(num_nodes * block_size, 2, 0.0);

            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                const double Ni = Ncontainer(g, i);
                const unsigned int base = i * block_size;

                NwT(base + 1, 0) = Ni;
                NwT(base + 2, 1) = Ni;
            }

            Vector traction_vec(2);
            traction_vec[0] = traction[0];
            traction_vec[1] = traction[1];

            Vector penalty_vec(2);
            penalty_vec[0] = gamma * m_dot_n * normal[0];
            penalty_vec[1] = gamma * m_dot_n * normal[1];

            Vector R_total = prod(NwT, traction_vec) + prod(NwT, penalty_vec);

            noalias(rRightHandSideVector) += R_total * coeff;

            if (this->Id() == 100004)
            {
                KRATOS_WATCH(R_total.size());
                KRATOS_WATCH(rRightHandSideVector);
            }
        }
    }

    // =======================
    // Explicit instantiation
    // =======================

    template class TestCondition<2, 2>;
    template class TestCondition<3, 3>;

} // namespace Kratos