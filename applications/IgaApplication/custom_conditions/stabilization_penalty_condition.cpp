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

        const double stabilization_translation_factor = GetProperties()[STABILIZATION_TRANSLATION_FACTOR];
        const double stabilization_rotation_factor = GetProperties()[STABILIZATION_ROTATION_FACTOR];
        const double initial_distance_to_first_stable_point = GetProperties()[INITIAL_DISTANCE_TO_FIRST_STABLE_POINT];
        const double initial_distance_first_to_second_stable_point = GetProperties()[INITIAL_DISTANCE_FIRST_TO_SECOND_STABLE_POINT];

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = r_geometry.WorkingSpaceDimension() * number_of_nodes;

        // Assembly
        if (CalculateStiffnessMatrixFlag) {
            // noalias(rLeftHandSideMatrix) += 0.0; //TO DO
        }
        if (CalculateResidualVectorFlag) {
            // noalias(rRightHandSideVector) -= 0.0; //TO DO
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
        const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

        #pragma omp critical
        {
            if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
                for(SizeType i=0; i< number_of_nodes; ++i) {
                    SizeType index = 6 * i;

                    array_1d<double, 3 >& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                    for(SizeType j = 0; j < dimension; ++j) {
                        AtomicAdd(r_force_residual[j], rRHS[index + j]);
                    }
                }
            }
            else if(rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == MOMENT_RESIDUAL) {
             
                for (IndexType i = 0; i < number_of_nodes; ++i) {
                    IndexType index = 6 * i;
                    array_1d<double, 3>& r_moment_residual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

                    for (IndexType j = 0; j < dimension; ++j) {
                        AtomicAdd(r_moment_residual[j], rRHS[index + j + 3]);
                    }
                }
            }
        }

        KRATOS_CATCH( "" )
    }


    int StabilizationPenaltyCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of StabilizationPenaltyCondition" << std::endl;
        return 0;
    }

    void StabilizationPenaltyCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 6 * number_of_nodes)
            rResult.resize(6 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
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
    };

} // Namespace Kratos
