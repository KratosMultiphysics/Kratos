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
#define  KRATOS_TRILINOS_RESIDUAL_CRITERIA_H_INCLUDED

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
template< class TSparseSpace, class TDenseSpace >
class TrilinosResidualCriteria : public ResidualCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosResidualCriteria
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosResidualCriteria);

    typedef ResidualCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef typename BaseType::TDataType TDataType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    explicit TrilinosResidualCriteria(TDataType NewRatioTolerance,TDataType AlwaysConvergedNorm):
        ResidualCriteria<TSparseSpace,TDenseSpace>(NewRatioTolerance, AlwaysConvergedNorm)
    {}

    /// Copy constructor
    explicit TrilinosResidualCriteria(const TrilinosResidualCriteria& rOther):
        ResidualCriteria<TSparseSpace,TDenseSpace>(rOther)
    {}

    /// Destructor.
    ~TrilinosResidualCriteria() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    TrilinosResidualCriteria& operator=(TrilinosResidualCriteria const& rOther) = delete;

    ///@}

protected:

    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the norm of the residual
     * @details It checks if the dof is fixed
     * @param rResidualSolutionNorm The norm of the residual
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param b RHS vector (residual + reactions)
     */
    void CalculateResidualNorm(
        TDataType& rResidualSolutionNorm,
        typename BaseType::SizeType& rDofNum,
        typename BaseType::DofsArrayType& rDofSet,
        const typename BaseType::TSystemVectorType& rB) override
    {
        // Initialize
        TDataType residual_solution_norm = TDataType();
        long int local_dof_num = 0;

        const double rank = rB.Comm().MyPID(); // To compare with PARTITION_INDEX, which is a double variable

        // Loop over Dofs
        #pragma omp parallel for reduction(+:residual_solution_norm,local_dof_num)
        for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
            auto it_dof = rDofSet.begin() + i;

            typename BaseType::IndexType dof_id;
            TDataType residual_dof_value;

            if (it_dof->IsFree() && (it_dof->GetSolutionStepValue(PARTITION_INDEX) == rank)) {
                dof_id = it_dof->EquationId();
                residual_dof_value = TSparseSpace::GetValue(rB,dof_id);
                residual_solution_norm += residual_dof_value * residual_dof_value;
                local_dof_num++;
            }
        }

        // Combine local contributions
        // Note that I'm not merging the two calls because one adds doubles and the other ints (JC)
        rB.Comm().SumAll(&residual_solution_norm,&rResidualSolutionNorm,1);
        // SizeType is long unsigned int in linux, but EpetraComm does not support unsigned types
        long int global_dof_num = 0;
        rB.Comm().SumAll(&local_dof_num,&global_dof_num,1);
        rDofNum = static_cast<typename BaseType::SizeType>(global_dof_num);

        rResidualSolutionNorm = std::sqrt(rResidualSolutionNorm);
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}

}; // Class TrilinosResidualCriteria

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_RESIDUAL_CRITERIA_H_INCLUDED  defined
