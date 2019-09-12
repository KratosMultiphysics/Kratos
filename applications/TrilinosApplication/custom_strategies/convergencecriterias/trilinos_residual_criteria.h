//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_TRILINOS_RESIDUAL_CRITERIA_H_INCLUDED)
#define KRATOS_TRILINOS_RESIDUAL_CRITERIA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"

namespace Kratos
{
///@addtogroup TrilinosApplication
///@{

///@name Kratos Classes
///@{

/// MPI version of the ResidualCriteria.
/** Implements a convergence criteria based on the norm of the (free rows of) the RHS vector.
 *  @see ResidualCriteria
 */
template <class TSparseSpace, class TDenseSpace>
class TrilinosResidualCriteria : public ResidualCriteria<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosResidualCriteria
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosResidualCriteria);

    typedef ResidualCriteria<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    explicit TrilinosResidualCriteria(TDataType NewRatioTolerance, TDataType AlwaysConvergedNorm)
        : ResidualCriteria<TSparseSpace, TDenseSpace>(NewRatioTolerance, AlwaysConvergedNorm)
    {
    }

    /// Copy constructor
    explicit TrilinosResidualCriteria(const TrilinosResidualCriteria& rOther)
        : ResidualCriteria<TSparseSpace, TDenseSpace>(rOther)
    {
    }

    /// Destructor.
    ~TrilinosResidualCriteria() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    TrilinosResidualCriteria& operator=(TrilinosResidualCriteria const& rOther) = delete;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void CalculateResidualNorm(ModelPart& rModelPart,
                               TDataType& rResidualSolutionNorm,
                               typename BaseType::SizeType& rDofNum,
                               typename BaseType::DofsArrayType& rDofSet,
                               const typename BaseType::TSystemVectorType& rB) override
    {
        // Initialize
        if (!this->mImportIsInitialized)
            this->InitializeMapping(rDofSet, rB);

        TDataType residual_solution_norm = TDataType();
        long int local_dof_num = 0;

        int system_size = TSparseSpace::Size(rB);

        // defining a temporary vector to gather all of the values needed
        Epetra_Vector local_dx(mpDofImport->TargetMap());

        // importing in the new temp vector the values
        int ierr = local_dx.Import(rB, *mpDofImport, Insert);
        KRATOS_ERROR_IF(ierr != 0)
            << "Epetra failure found while trying to import Dx." << std::endl;

        int num_dof = rDofSet.size();

// Loop over Dofs
#pragma omp parallel for reduction(+ : residual_solution_norm, local_dof_num)
        for (int i = 0; i < num_dof; i++)
        {
            auto it_dof = rDofSet.begin() + i;

            int dof_id;
            TDataType residual_dof_value;

            if (it_dof->IsFree())
            {
                dof_id = static_cast<int>(it_dof->EquationId());
                if (dof_id < system_size)
                {
                    residual_dof_value = local_dx[mpDofImport->TargetMap().LID(dof_id)];
                    residual_solution_norm += residual_dof_value * residual_dof_value;
                    local_dof_num++;
                }
            }
        }

        // Combine local contributions
        // Note that I'm not merging the two calls because one adds doubles and the other ints (JC)
        rB.Comm().SumAll(&residual_solution_norm, &rResidualSolutionNorm, 1);
        // SizeType is long unsigned int in linux, but EpetraComm does not support unsigned types
        long int global_dof_num = 0;
        rB.Comm().SumAll(&local_dof_num, &global_dof_num, 1);
        rDofNum = static_cast<typename BaseType::SizeType>(global_dof_num);

        rResidualSolutionNorm = std::sqrt(rResidualSolutionNorm);
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// This lets the class control if Initialize() was properly called.
    bool mImportIsInitialized = false;

    /// Auxiliary trilinos data structure to import out-of-process data in the update vector.
    std::unique_ptr<Epetra_Import> mpDofImport = nullptr;

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeMapping(typename BaseType::DofsArrayType& rDofSet,
                           const typename BaseType::TSystemVectorType& rB)
    {
        int system_size = TSparseSpace::Size(rB);
        int number_of_dofs = rDofSet.size();
        std::vector<int> index_array(number_of_dofs);

        // filling the array with the global ids
        unsigned int counter = 0;
        for (typename BaseType::DofsArrayType::const_iterator i_dof = rDofSet.begin();
             i_dof != rDofSet.end(); ++i_dof)
        {
            int id = i_dof->EquationId();
            if (id < system_size)
            {
                index_array[counter++] = id;
            }
        }

        std::sort(index_array.begin(), index_array.end());
        std::vector<int>::iterator new_end =
            std::unique(index_array.begin(), index_array.end());
        index_array.resize(new_end - index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        rB.Comm().SumAll(&tot_update_dofs, &check_size, 1);
        if ((check_size < system_size) && (rB.Comm().MyPID() == 0))
        {
            std::stringstream msg;
            msg << "Dof count is not correct. There are less dofs then "
                   "expected."
                << std::endl;
            msg << "Expected number of active dofs: " << system_size
                << ", dofs found: " << check_size << std::endl;
            KRATOS_ERROR << msg.str();
        }

        // defining a map as needed
        Epetra_Map dof_update_map(-1, index_array.size(),
                                  &(*(index_array.begin())), 0, rB.Comm());

        // defining the import instance
        std::unique_ptr<Epetra_Import> p_dof_import(
            new Epetra_Import(dof_update_map, rB.Map()));
        mpDofImport.swap(p_dof_import);

        mImportIsInitialized = true;
    }

    ///@}

}; // Class TrilinosResidualCriteria

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TRILINOS_RESIDUAL_CRITERIA_H_INCLUDED  defined
