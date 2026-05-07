//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela, Riccardo Rossi, Carlos Roig and Ruben Zorrilla
//

#pragma once

// System includes

// External includes

/* Epetra includes */
#include <Epetra_Vector.h>
#include <Epetra_Import.h>

/* Tpetra includes (only when Tpetra support is compiled in) */
#ifdef HAVE_TPETRA
#include <Tpetra_MultiVector.hpp>
#endif

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

/**
 * @class TrilinosMixedGenericCriteria
 * @ingroup TrilinosApplication
 * @brief Convergence criteria for mixed vector-scalar problems.
 * @details This class implements a convergence control based on a nodal vector variable and a nodal scalar variable. The error is evaluated separately for each of them, and relative and absolute tolerances for both must be specified.
 * @author Jordi Cotela, Riccardo Rossi, Carlos Roig and Ruben Zorrilla
 * @tparam TSparseSpace The sparse space considered (e.g. for the system matrix and the solution vector)
 * @tparam TDenseSpace The dense space considered (e.g. for the local element matrices and vectors)
 */
template< class TSparseSpace, class TDenseSpace >
class TrilinosMixedGenericCriteria : public MixedGenericCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of TrilinosMixedGenericCriteria
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosMixedGenericCriteria);

    /// The definition of the base class
    using BaseType = MixedGenericCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the current class
    using TDataType = typename BaseType::TDataType;

    /// DoF array type definition
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// System matrix type definition
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    /// Definition of the solution vector type
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    /// Convergence variable list type definition
    using ConvergenceVariableListType = typename BaseType::ConvergenceVariableListType;

    /// Key type definition (used for the map of the norms)
    using KeyType = typename BaseType::KeyType;

    /// Definition of the linear algebra library
    static constexpr TrilinosLinearAlgebraLibrary LinearAlgebraLibrary = TSparseSpace::LinearAlgebraLibrary();

    ///@}
    ///@name Life Cycle
    ///@{

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

    /**
     * @brief  Compute relative and absoute error. 
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
        // If required, initialize the vector import
        if(mpDofImport == nullptr) {
            mpDofImport = TSparseSpace::CreateImport(rDofSet, Dx);
        }

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

    typename TSparseSpace::ImportPointerType mpDofImport = nullptr; /// Importer for the Dx vector to be used in the convergence check

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Compute the norms of the solution and the increase for the convergence check
     * @details This method is called by the base class PostCriteria after importing the Dx vector to the local vector. It loops over the dofs, and for those that are free and belong to the current rank, it adds the contribution of the corresponding variable to the solution and increase norms.
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rDx Vector of results (variations on nodal variables) already imported to the local vector
     * @param rDofsCount Vector containing the number of dofs for each variable to be used in the averaging of the norms
     * @param rSolutionNormsVector Vector containing the solution norms for each variable. The contribution of each dof to the corresponding variable is added to the corresponding entry of this vector
     * @param rIncreaseNormsVector Vector containing the increase norms for each variable. The contribution of each dof to the corresponding variable is added to the corresponding entry of this vector
     */
    void GetNormValues(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rDx,
        std::vector<int>& rDofsCount,
        std::vector<TDataType>& rSolutionNormsVector,
        std::vector<TDataType>& rIncreaseNormsVector) const override
    {
        const std::size_t n_dofs = rDofSet.size();
        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int rank = r_data_comm.Rank();

        // Import Dx into a local vector and compute norms.
        // Use compile-time branching so this works for both EPETRA and TPETRA spaces.
        if constexpr (LinearAlgebraLibrary == TrilinosLinearAlgebraLibrary::EPETRA) {
            // Epetra implementation
            Epetra_Vector local_dx(mpDofImport->TargetMap());
            int i_err = local_dx.Import(rDx, *mpDofImport, Insert);
            KRATOS_ERROR_IF_NOT(i_err == 0) << "Local Dx import failed!" << std::endl;

            // Local thread variables
            int dof_id;
            TDataType dof_dx;

            // Loop over Dofs
            for (std::size_t i = 0; i < n_dofs; i++) {
                auto it_dof = rDofSet.begin() + i;
                if (it_dof->IsFree() && it_dof->GetSolutionStepValue(PARTITION_INDEX) == rank) {
                    dof_id = it_dof->EquationId();
                    const TDataType& r_dof_value = it_dof->GetSolutionStepValue(0);
                    dof_dx = local_dx[mpDofImport->TargetMap().LID(dof_id)];

                    int var_local_key;
                    bool key_found = BaseType::FindVarLocalKey(it_dof,var_local_key);
                    if (!key_found) {
                        continue;
                    }

                    rSolutionNormsVector[var_local_key] += r_dof_value * r_dof_value;
                    rIncreaseNormsVector[var_local_key] += dof_dx * dof_dx;
                    rDofsCount[var_local_key]++;
                }
            }
        } else {
            KRATOS_ERROR << "Only EPETRA is supported for now" << std::endl;
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
}; /// TrilinosMixedGenericCriteria
///@} // Kratos classes

///@} // Application group
}