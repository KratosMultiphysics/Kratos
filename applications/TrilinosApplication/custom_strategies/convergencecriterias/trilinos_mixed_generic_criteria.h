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

#ifndef KRATOS_TRILINOS_MIXED_GENERIC_CRITERIA_H
#define	KRATOS_TRILINOS_MIXED_GENERIC_CRITERIA_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
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

    typedef MixedGenericCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::ConvergenceVariableListType ConvergenceVariableListType;

    typedef typename BaseType::KeyType KeyType;

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
        // If required, initialize the Epetra vector import
        if(!mEpetraImportIsInitialized) {
            InitializeEpetraImport(rDofSet, Dx);
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

    bool mEpetraImportIsInitialized = false;
    std::unique_ptr<Epetra_Import> mpDofImport = nullptr;

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
        int n_dofs = rDofSet.size();
        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int rank = r_data_comm.Rank();

        // Do the local Dx vector import
        Epetra_Vector local_dx(mpDofImport->TargetMap());
        int i_err = local_dx.Import(rDx, *mpDofImport, Insert);
        KRATOS_ERROR_IF_NOT(i_err == 0) << "Local Dx import failed!" << std::endl;

        // Local thread variables
        int dof_id;
        TDataType dof_dx;

        // Loop over Dofs
        for (int i = 0; i < n_dofs; i++) {
            auto it_dof = rDofSet.begin() + i;
            if (it_dof->IsFree() && it_dof->GetSolutionStepValue(PARTITION_INDEX) == rank) {
                dof_id = it_dof->EquationId();
                const TDataType& r_dof_value = it_dof->GetSolutionStepValue(0);
                dof_dx = local_dx[mpDofImport->TargetMap().LID(dof_id)];

                int var_local_key;
                bool key_found = BaseType::FindVarLocalKey(it_dof,var_local_key);
                if (!key_found) {
                    // the dof does not belong to the list of variables
                    // we are checking for convergence, so we skip it
                    continue;
                }

                rSolutionNormsVector[var_local_key] += r_dof_value * r_dof_value;
                rIncreaseNormsVector[var_local_key] += dof_dx * dof_dx;
                rDofsCount[var_local_key]++;
            }
        }
    }

    /**
     * @brief Initialize the DofUpdater in preparation for a subsequent UpdateDofs call.
     *  The DofUpdater needs to be initialized only if the dofset changes.
     *  If the problem does not require creating/destroying nodes or changing the
     *  mesh graph, it is in general enough to intialize this tool once at the
     *  begining of the problem.
     *  If the dofset only changes under certain conditions (for example because
     *  the domain is remeshed every N iterations), it is enough to call the
     *  Clear method to let this class know that its auxiliary data has to be re-generated
     *  and Initialize will be called as part of the next UpdateDofs call.
     * @param rDofSet rDofSet The list of degrees of freedom.
     * @param rDx rDx The update vector.
     */
    void InitializeEpetraImport(
        const DofsArrayType &rDofSet,
        const TSystemVectorType &rDx)
    {
        int number_of_dofs = rDofSet.size();
        int system_size = TSparseSpace::Size(rDx);
        std::vector<int> index_array(number_of_dofs);

        // Filling the array with the global ids
        unsigned int counter = 0;
        for (typename DofsArrayType::const_iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof) {
            const int id = i_dof->EquationId();
            if (id < system_size) {
                index_array[counter++] = id;
            }
        }

        std::sort(index_array.begin(), index_array.end());
        std::vector<int>::iterator new_end = std::unique(index_array.begin(), index_array.end());
        index_array.resize(new_end - index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        rDx.Comm().SumAll(&tot_update_dofs, &check_size, 1);
        KRATOS_ERROR_IF(check_size < system_size)
            << "DOF count is not correct. There are less dofs then expected.\n"
            << "Expected number of active dofs: " << system_size << ", DOFs found: " << check_size << std::endl;

        // Defining a map as needed
        Epetra_Map dof_update_map(-1, index_array.size(), &(*(index_array.begin())), 0, rDx.Comm());

        // Defining the import instance
        std::unique_ptr<Epetra_Import> p_dof_import(new Epetra_Import(dof_update_map, rDx.Map()));
        mpDofImport.swap(p_dof_import);

        mEpetraImportIsInitialized = true;
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

#endif // KRATOS_TRILINOS_MIXED_GENERIC_CRITERIA_H