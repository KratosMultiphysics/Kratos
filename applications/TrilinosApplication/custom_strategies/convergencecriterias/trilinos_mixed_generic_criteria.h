//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela, Riccardo Rossi, Carlos Roig and Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"

// Application includes

namespace Kratos
{
///@addtogroup TrilinosApplication
///@{

///@name Kratos Classes
///@{

/// Convergence criteria for mixed vector-scalar problems.
/**
 This class implements a convergence control based on a nodal vector variable and
 a nodal scalar variable. The error is evaluated separately for each of them, and
 relative and absolute tolerances for both must be specified.
 */
template< class TSparseSpace, class TDenseSpace >
class TrilinosMixedGenericCriteria : public MixedGenericCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(TrilinosMixedGenericCriteria);

    using BaseType = MixedGenericCriteria< TSparseSpace, TDenseSpace >;

    using TDataType = typename BaseType::TDataType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using ConvergenceVariableListType = typename BaseType::ConvergenceVariableListType;

    using KeyType = typename BaseType::KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    /**
     * @brief Construct a new Trilinos Mixed Generic Criteria object
     * Construct the mixed generic convergence criteria from a convergence variables list.
     * The convergence variable list contains for each variable the variable itself as well as the corresponding relative and absolute tolerances.
     * @param rConvergenceVariablesList List containing tuples with the convergence variables to be checked. The tuples are set as <Variable, relative tolerance, absolute tolerance>
     */
    TrilinosMixedGenericCriteria(const ConvergenceVariableListType& rConvergenceVariablesList)
        : MixedGenericCriteria<TSparseSpace, TDenseSpace>(rConvergenceVariablesList)
    {}

    /// Destructor.
    ~TrilinosMixedGenericCriteria() override
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Compute relative and absoute error.
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b) override
    {
        // Serial base implementation convergence check call
        return BaseType::PostCriteria(rModelPart, rDofSet, A, Dx, b);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void GetNormValues(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rDx,
        std::vector<int>& rDofsCount,
        std::vector<TDataType>& rSolutionNormsVector,
        std::vector<TDataType>& rIncreaseNormsVector) const override
    {
        const int n_dofs = rDofSet.size();
        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int rank = r_data_comm.Rank();

        std::vector<int> equation_ids;
        equation_ids.reserve(n_dofs);

        std::vector<int> var_local_keys;
        var_local_keys.reserve(n_dofs);

        std::vector<TDataType> dof_values;
        dof_values.reserve(n_dofs);

        // Gather equation ids and corresponding information for free dofs in this rank
        for (int i = 0; i < n_dofs; i++) {
            auto it_dof = rDofSet.begin() + i;
            if (it_dof->IsFree() && it_dof->GetSolutionStepValue(PARTITION_INDEX) == rank) {
                int var_local_key;
                bool key_found = BaseType::FindVarLocalKey(it_dof,var_local_key);
                if (!key_found) {
                    // the dof does not belong to the list of variables
                    // we are checking for convergence, so we skip it
                    continue;
                }

                equation_ids.push_back(static_cast<int>(it_dof->EquationId()));
                var_local_keys.push_back(var_local_key);
                dof_values.push_back(it_dof->GetSolutionStepValue(0));
            }
        }

        std::vector<double> dof_increments(equation_ids.size(), 0.0);
        if (!equation_ids.empty()) {
            TSparseSpace::GatherValues(rDx, equation_ids, dof_increments.data());
        }

        const std::size_t n_local_dofs = equation_ids.size();
        for (std::size_t i = 0; i < n_local_dofs; ++i) {
            const int var_local_key = var_local_keys[i];
            const TDataType& r_dof_value = dof_values[i];
            const TDataType dof_dx = static_cast<TDataType>(dof_increments[i]);

            rSolutionNormsVector[var_local_key] += r_dof_value * r_dof_value;
            rIncreaseNormsVector[var_local_key] += dof_dx * dof_dx;
            rDofsCount[var_local_key]++;
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
};
///@} // Kratos classes

///@} // Application group
}