//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <array>
#include <limits>
#include <vector>

// Project includes
#include "custom_conditions/stabilization_penalty_condition_3p.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
    void StabilizationPenaltyCondition3P::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const double stabilization_translation_factor =
            GetProperties()[STABILIZATION_TRANSLATION_FACTOR];

        auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        KRATOS_ERROR_IF(number_of_nodes < 3 || number_of_nodes % 2 == 0)
            << "StabilizationPenaltyCondition3P requires node ordering "
            << "[A, R1_1, R2_1, ..., R1_nd, R2_nd], "
            << "therefore the number of nodes must be odd and >= 3. Got "
            << number_of_nodes << std::endl;

        const SizeType number_of_directions = (number_of_nodes - 1) / 2;
        const double inv_nd = 1.0 / static_cast<double>(number_of_directions);

        const SizeType mat_size = 3 * number_of_nodes;

        if (CalculateStiffnessMatrixFlag) {
            if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size) {
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
        }

        if (CalculateResidualVectorFlag) {
            if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            noalias(rRightHandSideVector) = ZeroVector(mat_size);
        }

        // node 0 = A
        // nodes 1,2 = R1_1, R2_1
        // nodes 3,4 = R1_2, R2_2
        // ...
        std::vector<double> c_t(number_of_nodes, 0.0);
        c_t[0] = 1.0;

        array_1d<double, 3> XA;
        XA[0] = r_geometry[0].X0();
        XA[1] = r_geometry[0].Y0();
        XA[2] = r_geometry[0].Z0();

        for (SizeType k = 0; k < number_of_directions; ++k) {
            const SizeType iR1 = 1 + 2 * k;
            const SizeType iR2 = 2 + 2 * k;

            array_1d<double, 3> XR1;
            array_1d<double, 3> XR2;

            XR1[0] = r_geometry[iR1].X0();
            XR1[1] = r_geometry[iR1].Y0();
            XR1[2] = r_geometry[iR1].Z0();

            XR2[0] = r_geometry[iR2].X0();
            XR2[1] = r_geometry[iR2].Y0();
            XR2[2] = r_geometry[iR2].Z0();

            array_1d<double, 3> dA1;
            array_1d<double, 3> d12;

            for (IndexType d = 0; d < 3; ++d) {
                dA1[d] = XA[d] - XR1[d];
                d12[d] = XR1[d] - XR2[d];
            }

            const double L1  = norm_2(dA1);
            const double L12 = norm_2(d12);

            KRATOS_ERROR_IF(L12 <= std::numeric_limits<double>::epsilon())
                << "Distance L12 must be > 0.0 for direction " << k
                << " in StabilizationPenaltyCondition3P" << std::endl;

            const double r = L1 / L12;

            // g_t = u_A - (1/n_d) * sum_k [ (1+r_k) u_R1,k - r_k u_R2,k ]
            c_t[iR1] += -(1.0 + r) * inv_nd;
            c_t[iR2] +=  (r)       * inv_nd;
        }

        array_1d<double, 3> g_t = ZeroVector(3);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const array_1d<double, 3>& r_disp =
                r_geometry[i].FastGetCurrentSolutionStepValue(DISPLACEMENT, 0);

            for (IndexType d = 0; d < 3; ++d) {
                g_t[d] += c_t[i] * r_disp[d];
            }
        }

        const double gt_norm = norm_2(g_t);
        KRATOS_WATCH(number_of_directions)
        KRATOS_WATCH(gt_norm)

        // K = beta_t * (dg/dd)^T * (dg/dd)
        if (CalculateStiffnessMatrixFlag) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType iu = 3 * i; // ux, uy, uz

                for (IndexType j = 0; j < number_of_nodes; ++j) {
                    const IndexType ju = 3 * j;

                    const double k_t = stabilization_translation_factor * c_t[i] * c_t[j];

                    for (IndexType d = 0; d < 3; ++d) {
                        rLeftHandSideMatrix(iu + d, ju + d) += k_t;
                    }
                }
            }
        }

        // RHS = - beta_t * (dg/dd)^T * g_t
        if (CalculateResidualVectorFlag) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType iu = 3 * i;

                for (IndexType d = 0; d < 3; ++d) {
                    rRightHandSideVector[iu + d] -=
                        stabilization_translation_factor * c_t[i] * g_t[d];
                }
            }

            const double rhs_norm = norm_2(rRightHandSideVector);
            KRATOS_WATCH(rhs_norm)
        }

        KRATOS_CATCH("")
    }

    void StabilizationPenaltyCondition3P::AddExplicitContribution(
        const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();

        #pragma omp critical
        {
            if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
                for (SizeType i = 0; i < number_of_nodes; ++i) {
                    const SizeType index = 3 * i;

                    array_1d<double, 3>& r_force_residual =
                        GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                    for (SizeType j = 0; j < dimension; ++j) {
                        AtomicAdd(r_force_residual[j], rRHS[index + j]);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    int StabilizationPenaltyCondition3P::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(STABILIZATION_TRANSLATION_FACTOR))
            << "Missing STABILIZATION_TRANSLATION_FACTOR in property of StabilizationPenaltyCondition3P" << std::endl;

        const SizeType number_of_nodes = GetGeometry().size();

        KRATOS_ERROR_IF(number_of_nodes < 3 || number_of_nodes % 2 == 0)
            << "StabilizationPenaltyCondition3P requires node ordering "
            << "[A, R1_1, R2_1, ..., R1_nd, R2_nd], "
            << "therefore the number of nodes must be odd and >= 3. Got "
            << number_of_nodes << std::endl;

        return 0;
    }

    void StabilizationPenaltyCondition3P::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes) {
            rResult.resize(3 * number_of_nodes, false);
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];

            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    void StabilizationPenaltyCondition3P::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    }

} // namespace Kratos