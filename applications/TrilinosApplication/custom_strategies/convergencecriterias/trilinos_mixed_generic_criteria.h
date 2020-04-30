//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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
class TrilinosMixedGenericCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(TrilinosMixedGenericCriteria);

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename DofUpdater<TSparseSpace >::UniquePointer DofUpdaterPointerType;

    typedef std::vector<std::tuple<VariableData, TDataType, TDataType>> ConvergenceVariableListType;

    typedef std::size_t KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    /**
     * @brief Construct a new Trilinos Mixed Generic Criteria object
     *
     * @param rConvergenceVariablesList List containing tuples with the convergence variables to be checked. The tuples are set as <Variable, relative tolerance, absolute tolerance>
     */
    TrilinosMixedGenericCriteria(const ConvergenceVariableListType& rConvergenceVariablesList)
        : ConvergenceCriteria<TSparseSpace, TDenseSpace>()
        , mVariableSize([&] (const ConvergenceVariableListType& rList) -> int {return rList.size();} (rConvergenceVariablesList))
        , mVariableDataVector([&] (const ConvergenceVariableListType& rList) -> std::vector<std::shared_ptr<VariableData>> {
            int i = 0;
            std::vector<std::shared_ptr<VariableData>> aux_vect(mVariableSize);
            for (const auto &r_tup : rList) {
                aux_vect[i++] = Kratos::make_shared<VariableData>(std::get<0>(r_tup));
            }
            return aux_vect;
        } (rConvergenceVariablesList))
        , mRatioToleranceVector([&] (const ConvergenceVariableListType& rList) -> std::vector<TDataType> {
            int i = 0;
            std::vector<TDataType> aux_vect(mVariableSize);
            for (const auto &r_tup : rList) {
                aux_vect[i++] = std::get<1>(r_tup);
            }
            return aux_vect;
        } (rConvergenceVariablesList))
        , mAbsToleranceVector([&] (const ConvergenceVariableListType& rList) -> std::vector<TDataType> {
            int i = 0;
            std::vector<TDataType> aux_vect(mVariableSize);
            for (const auto &r_tup : rList) {
                aux_vect[i++] = std::get<2>(r_tup);
            }
            return aux_vect;
        } (rConvergenceVariablesList))
        , mLocalKeyMap([&] (const ConvergenceVariableListType& rList) -> std::unordered_map<KeyType, KeyType> {
            KeyType local_key = 0;
            std::unordered_map<KeyType, KeyType> aux_map;
            for (const auto &r_tup : rList) {
                const auto *p_var_data = &(std::get<0>(r_tup));
                aux_map[p_var_data->Key()] = local_key++;
            }
            return aux_map;
        } (rConvergenceVariablesList))
    {}

    /// Destructor.
    ~TrilinosMixedGenericCriteria() override {}

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

        // Check if we are solving for something
        if (TSparseSpace::Size(Dx) != 0) {

            // Initialize
            std::vector<int> dofs_count(mVariableSize, 0);
            std::vector<TDataType> solution_norms_vector(mVariableSize, 0.0);
            std::vector<TDataType> increase_norms_vector(mVariableSize, 0.0);

            int n_dofs = rDofSet.size();
            const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
            const int rank = r_data_comm.Rank();

            // Do the local Dx vector import
            Epetra_Vector local_dx( mpDofImport->TargetMap() );
            int i_err = local_dx.Import(Dx, *mpDofImport, Insert);
            KRATOS_ERROR_IF_NOT(i_err == 0) << "Local Dx import failed!" << std::endl;

            // Loop over Dofs
#pragma omp parallel
            {
                // Local thread variables
                int dof_id;
                TDataType dof_dx;
                TDataType dof_value;

                // Local reduction variables
                std::vector<TDataType> var_solution_norm_reduction(mVariableSize);
                std::vector<TDataType> var_correction_norm_reduction(mVariableSize);
                std::vector<int> dofs_counter_reduction(mVariableSize);
                for (int i = 0; i < mVariableSize; i++) {
                    var_solution_norm_reduction[i] = 0.0;
                    var_correction_norm_reduction[i] = 0.0;
                    dofs_counter_reduction[i] = 0;
                }

#pragma omp for
                for (int i = 0; i < n_dofs; i++) {
                    auto it_dof = rDofSet.begin() + i;
                    if (it_dof->IsFree() && it_dof->GetSolutionStepValue(PARTITION_INDEX) == rank) {
                        dof_id = it_dof->EquationId();
                        dof_value = it_dof->GetSolutionStepValue(0);
                        dof_dx = local_dx[mpDofImport->TargetMap().LID(dof_id)];

                        const auto &r_current_variable = it_dof->GetVariable();
                        int var_local_key = mLocalKeyMap[r_current_variable.IsComponent() ? r_current_variable.GetSourceVariable().Key() : r_current_variable.Key()];

                        var_solution_norm_reduction[var_local_key] += dof_value * dof_value;
                        var_correction_norm_reduction[var_local_key] += dof_dx * dof_dx;
                        dofs_counter_reduction[var_local_key]++;
                    }
                }

#pragma omp critical
                {
                    for (int i = 0; i < mVariableSize; i++) {
                        solution_norms_vector[i] += var_solution_norm_reduction[i];
                        increase_norms_vector[i] += var_correction_norm_reduction[i];
                        dofs_count[i] += dofs_counter_reduction[i];
                    }
                }
            }

            auto global_solution_norms_vector = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(solution_norms_vector);
            auto global_increase_norms_vector = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(increase_norms_vector);
            auto global_dofs_count = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(dofs_count);

            const double zero_tol = 1.0e-12;
            for(int i = 0; i < mVariableSize; i++) {
                if (global_solution_norms_vector[i] < zero_tol) {
                    global_solution_norms_vector[i] = 1.0;
                }
            }

            std::vector<TDataType> var_ratio(mVariableSize, 0.0);
            std::vector<TDataType> var_abs(mVariableSize, 0.0);
            for(int i = 0; i < mVariableSize; i++) {
                var_ratio[i] = std::sqrt(global_increase_norms_vector[i]) / std::sqrt(global_solution_norms_vector[i]);
                var_abs[i] = std::sqrt(global_increase_norms_vector[i]) / static_cast<TDataType>(global_dofs_count[i]);
            }

            // Output convergence status
            if (this->GetEchoLevel() > 0) {
                std::ostringstream stringbuf;
                stringbuf << "CONVERGENCE CHECK:\n";
                for(int i = 0; i < mVariableSize; i++) {
                    const auto r_var_data = mVariableDataVector[i];
                    const int key_map = mLocalKeyMap[r_var_data->Key()];
                    stringbuf << " " << r_var_data->Name() << " : ratio = " << var_ratio[key_map] << "; exp.ratio = " << mRatioToleranceVector[key_map] << " abs = " << var_abs[key_map] << " exp.abs = " << mAbsToleranceVector[key_map] << "\n";
                }
                KRATOS_INFO("") << stringbuf.str();
            }

            // Check convergence
            int is_converged = 0;
            for (int i = 0; i < mVariableSize; i++) {
                const auto r_var_data = mVariableDataVector[i];
                const int key_map = mLocalKeyMap[r_var_data->Key()];
                is_converged += var_ratio[key_map] <= mRatioToleranceVector[key_map] || var_abs[key_map] <= mAbsToleranceVector[key_map];
            }

            if (is_converged) {
                KRATOS_INFO_IF("", this->GetEchoLevel() > 0) << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
                return true;
            } else {
                return false;
            }
        } else {
            // Case in which all the DOFs are constrained!
            return true;
        }
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

    const int mVariableSize;
    const std::vector<std::shared_ptr<VariableData>> mVariableDataVector;
    const std::vector<TDataType> mRatioToleranceVector;
    const std::vector<TDataType> mAbsToleranceVector;
    std::unordered_map<KeyType, KeyType> mLocalKeyMap;

    bool mEpetraImportIsInitialized = false;
    std::unique_ptr<Epetra_Import> mpDofImport = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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