//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//                   Tobias Teschemacher
//

// System includes
#include <array>
#include <limits>

// External includes

// Project includes
#include "custom_conditions/stabilization_penalty_condition.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
    void StabilizationPenaltyCondition::CalculateAll(
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
        const double stabilization_rotation_factor =
            GetProperties()[STABILIZATION_ROTATION_FACTOR];
        const double initial_distance_to_first_stable_point =
            GetProperties()[INITIAL_DISTANCE_TO_FIRST_STABLE_POINT];
        const double initial_distance_first_to_second_stable_point =
            GetProperties()[INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT];

        auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        KRATOS_ERROR_IF(number_of_nodes != 3)
            << "StabilizationPenaltyCondition requires exactly 3 nodes [A, R1, R2]." << std::endl;

        const SizeType mat_size = 6 * number_of_nodes;

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
        // g_r = th_A - th_R1
        const std::array<double, 3> c_t{{1.0, -(1.0 + r), r}};
        const std::array<double, 3> c_r{{1.0, -1.0, 0.0}};

        array_1d<double, 3> g_t = ZeroVector(3);
        array_1d<double, 3> g_r = ZeroVector(3);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const array_1d<double, 3>& r_disp =
                r_geometry[i].FastGetCurrentSolutionStepValue(DISPLACEMENT, 0);

            array_1d<double, 3> r_rot;
            r_rot[0] = r_geometry[i].FastGetCurrentSolutionStepValue(ROTATION_X, 0);
            r_rot[1] = r_geometry[i].FastGetCurrentSolutionStepValue(ROTATION_Y, 0);
            r_rot[2] = r_geometry[i].FastGetCurrentSolutionStepValue(ROTATION_Z, 0);

            for (IndexType d = 0; d < 3; ++d) {
                g_t[d] += c_t[i] * r_disp[d];
                g_r[d] += c_r[i] * r_rot[d];
            }
        }

       const double gt_norm = norm_2(g_t);
       const double gr_norm = norm_2(g_r);
       KRATOS_WATCH(gt_norm)
       KRATOS_WATCH(gr_norm)

        // K = beta * (dg/dd)^T * (dg/dd)
        if (CalculateStiffnessMatrixFlag) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType iu = 6 * i;       // ux, uy, uz
                const IndexType ir = 6 * i + 3;   // rx, ry, rz

                for (IndexType j = 0; j < number_of_nodes; ++j) {
                    const IndexType ju = 6 * j;
                    const IndexType jr = 6 * j + 3;

                    const double k_t = stabilization_translation_factor * c_t[i] * c_t[j];
                    const double k_r = stabilization_rotation_factor * c_r[i] * c_r[j];

                    for (IndexType d = 0; d < 3; ++d) {
                        rLeftHandSideMatrix(iu + d, ju + d) += k_t;
                        rLeftHandSideMatrix(ir + d, jr + d) += k_r;
                    }
                }
            }

            KRATOS_WATCH(rLeftHandSideMatrix(0,0))
            KRATOS_WATCH(rLeftHandSideMatrix(0,6))
            KRATOS_WATCH(rLeftHandSideMatrix(0,12))
            KRATOS_WATCH(rLeftHandSideMatrix(3,3))
            KRATOS_WATCH(rLeftHandSideMatrix(3,9))
            KRATOS_WATCH(rLeftHandSideMatrix(3,15))
        }

        // RHS = - beta * (dg/dd)^T * g
        if (CalculateResidualVectorFlag) {
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType iu = 6 * i;
                const IndexType ir = 6 * i + 3;

                for (IndexType d = 0; d < 3; ++d) {
                    rRightHandSideVector[iu + d] -=
                        stabilization_translation_factor * c_t[i] * g_t[d];
                    rRightHandSideVector[ir + d] -=
                        stabilization_rotation_factor * c_r[i] * g_r[d];
                }
            }

            const double rhs_norm = norm_2(rRightHandSideVector);
            KRATOS_WATCH(rhs_norm)
        }

        KRATOS_CATCH("")
    }

    void StabilizationPenaltyCondition::AddExplicitContribution(
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
                    const SizeType index = 6 * i;

                    array_1d<double, 3>& r_force_residual =
                        GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                    for (SizeType j = 0; j < dimension; ++j) {
                        AtomicAdd(r_force_residual[j], rRHS[index + j]);
                    }
                }
            }
            else if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == MOMENT_RESIDUAL) {
                for (IndexType i = 0; i < number_of_nodes; ++i) {
                    const IndexType index = 6 * i;

                    array_1d<double, 3>& r_moment_residual =
                        GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

                    for (IndexType j = 0; j < dimension; ++j) {
                        AtomicAdd(r_moment_residual[j], rRHS[index + j + 3]);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    int StabilizationPenaltyCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(STABILIZATION_TRANSLATION_FACTOR))
            << "Missing STABILIZATION_TRANSLATION_FACTOR in property of StabilizationPenaltyCondition" << std::endl;

        KRATOS_ERROR_IF_NOT(GetProperties().Has(STABILIZATION_ROTATION_FACTOR))
            << "Missing STABILIZATION_ROTATION_FACTOR in property of StabilizationPenaltyCondition" << std::endl;

        KRATOS_ERROR_IF_NOT(GetProperties().Has(INITIAL_DISTANCE_TO_FIRST_STABLE_POINT))
            << "Missing INITIAL_DISTANCE_TO_FIRST_STABLE_POINT in property of StabilizationPenaltyCondition" << std::endl;

        KRATOS_ERROR_IF_NOT(GetProperties().Has(INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT))
            << "Missing INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT in property of StabilizationPenaltyCondition" << std::endl;

        KRATOS_ERROR_IF(GetGeometry().size() != 3)
            << "StabilizationPenaltyCondition requires exactly 3 nodes [A, R1, R2]." << std::endl;

        return 0;
    }

    void StabilizationPenaltyCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 6 * number_of_nodes) {
            rResult.resize(6 * number_of_nodes, false);
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 6;
            const auto& r_node = r_geometry[i];

            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
        }
    }

    void StabilizationPenaltyCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_X));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Y));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Z));
        }
    }

} // namespace Kratos






