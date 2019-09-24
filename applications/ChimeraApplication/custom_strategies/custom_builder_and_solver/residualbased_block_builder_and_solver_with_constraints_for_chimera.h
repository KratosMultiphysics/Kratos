//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
//
#if !defined(KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS_FOR_CHIMERA)
#define KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS__RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS_FOR_CHIMERA

/* System includes */
#include <unordered_set>
#include <unordered_map>

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/master_slave_constraint.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera
 * @ingroup KratosChimeraApplication
 * @brief Current class provides an implementation for applying the chimera constraints that is enforcing continuity.
 * @details This implementation enforces continuity necessary for chimera in the following way :
 *
 *
 *
 *   L =        [I  0  0 ]
 *              [0  I  0 ]
 *
 *  K_mod = L'KT
 *  F_mod = L'(F-K*g)
 * 
 *  Where T has the same definition as that of the classical master slave constraints
 *
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver>
class KRATOS_API(CHIMERA_APPLICATION) ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera
    : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base class
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    // The size_t types
    typedef typename BaseType::IndexType IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;


    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera);

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();
        mL.resize(0,0,false);
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    TSystemMatrixType mL; /// This is L matrix described above (at class definition)
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name ProtectedM_LN2
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ConstructMasterSlaveConstraintsStructure(ModelPart &rModelPart) override
    {
        BaseType::ConstructMasterSlaveConstraintsStructure(rModelPart);
        if (rModelPart.MasterSlaveConstraints().size() > 0)
        {
            mL = BaseType::mT;
        }
    }

    void BuildMasterSlaveConstraints(ModelPart &rModelPart) override
    {

        KRATOS_TRY

        BaseType::BuildMasterSlaveConstraints(rModelPart);

        // Setting the master dofs into the T and C system
        for (auto eq_id : BaseType::mMasterIds)
        {
            mL(eq_id, eq_id) = 1.0;
        }

        // Setting inactive slave dofs in the T and C system
        for (auto eq_id : BaseType::mInactiveSlaveDofs)
        {
            mL(eq_id, eq_id) = 1.0;
        }

        KRATOS_CATCH("")
    }

    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0)
        {
            BuildMasterSlaveConstraints(rModelPart);
            // We compute the transposed matrix of the global relation matrix
            // TSystemMatrixType L_transpose_matrix(mL.size2(), mL.size1());
            // auto L_transpose_matrix =  boost::numeric::ublas::trans(mL);
            // SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(L_transpose_matrix, mL, 1.0);

            TSystemVectorType b_modified(rb.size());
            boost::numeric::ublas::axpy_prod(boost::numeric::ublas::trans(mL), rb, b_modified, true);
            //TSparseSpace::Mult(L_transpose_matrix, rb, b_modified);
            TSparseSpace::Copy(b_modified, rb);
            b_modified.resize(0, false); //free memory

            TSystemMatrixType auxiliar_A_matrix(BaseType::mT.size2(), rA.size2());
            boost::numeric::ublas::axpy_prod(boost::numeric::ublas::trans(mL), rA, auxiliar_A_matrix, true);
            //SparseMatrixMultiplicationUtility::MatrixMultiplication(L_transpose_matrix, rA, auxiliar_A_matrix); //auxiliar = T_transpose * rA
            // L_transpose_matrix.resize(0, 0, false);                                                             //free memory

            boost::numeric::ublas::axpy_prod(auxiliar_A_matrix, BaseType::mT, rA, true);
            //SparseMatrixMultiplicationUtility::MatrixMultiplication(auxiliar_A_matrix, BaseType::mT, rA); //A = auxilar * T   NOTE: here we are overwriting the old A matrix!
            auxiliar_A_matrix.resize(0, 0, false);                                                        //free memory

            double max_diag = 0.0;
            for (IndexType i = 0; i < rA.size1(); ++i)
            {
                max_diag = std::max(std::abs(rA(i, i)), max_diag);
            }
// Apply diagonal values on slaves  BaseType::mDofSet.size()
#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(BaseType::mSlaveIds.size()); ++i)
            {
                const IndexType slave_equation_id = BaseType::mSlaveIds[i];
                if (BaseType::mInactiveSlaveDofs.find(slave_equation_id) == BaseType::mInactiveSlaveDofs.end())
                {
                    rA(slave_equation_id, slave_equation_id) = max_diag;
                    rb[slave_equation_id] = 0.0;
                }
            }
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
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

}; /* Class ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
