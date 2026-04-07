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

        const double initial_distance_to_first_stable_point =
            GetProperties()[INITIAL_DISTANCE_TO_FIRST_STABLE_POINT];

        const double initial_distance_first_to_second_stable_point =
            GetProperties()[INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT];

        auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        KRATOS_ERROR_IF(number_of_nodes != 3)
            << "StabilizationPenaltyCondition3P requires exactly 3 nodes [A, R1, R2]." << std::endl;

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

        KRATOS_ERROR_IF(initial_distance_first_to_second_stable_point <= std::numeric_limits<double>::epsilon())
            << "INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT must be > 0.0" << std::endl;

        const double r = initial_distance_to_first_stable_point /
                         initial_distance_first_to_second_stable_point;

        // node order must be [A, R1, R2]
        // g_t = u_A - (1+r) u_R1 + r u_R2
        const std::array<double, 3> c_t{{1.0, -(1.0 + r), r}};

        array_1d<double, 3> g_t = ZeroVector(3);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const array_1d<double, 3>& r_disp =
                r_geometry[i].FastGetCurrentSolutionStepValue(DISPLACEMENT, 0);

            for (IndexType d = 0; d < 3; ++d) {
                g_t[d] += c_t[i] * r_disp[d];
            }
        }

        const double gt_norm = norm_2(g_t);
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

            KRATOS_WATCH(rLeftHandSideMatrix(0,0))
            KRATOS_WATCH(rLeftHandSideMatrix(0,3))
            KRATOS_WATCH(rLeftHandSideMatrix(0,6))
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

        KRATOS_ERROR_IF_NOT(GetProperties().Has(INITIAL_DISTANCE_TO_FIRST_STABLE_POINT))
            << "Missing INITIAL_DISTANCE_TO_FIRST_STABLE_POINT in property of StabilizationPenaltyCondition3P" << std::endl;

        KRATOS_ERROR_IF_NOT(GetProperties().Has(INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT))
            << "Missing INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT in property of StabilizationPenaltyCondition3P" << std::endl;

        KRATOS_ERROR_IF(GetGeometry().size() != 3)
            << "StabilizationPenaltyCondition3P requires exactly 3 nodes [A, R1, R2]." << std::endl;

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