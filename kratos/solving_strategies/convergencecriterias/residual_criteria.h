//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    

#if !defined(KRATOS_NEW_RESIDUAL_CRITERIA )
#define  KRATOS_NEW_RESIDUAL_CRITERIA

// System includes 

// External includes 

// Project includes 
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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
 * @class ResidualCriteria
 * @ingroup KratosCore
 * @brief This is a convergence criteria that employes the residual as criteria
 * @details The reactions from the RHS are not computed in the residual
 * @author Riccardo Rossi
*/
template<class TSparseSpace,
         class TDenseSpace
         >
class ResidualCriteria
    : public  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions 
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ResidualCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@} 
    ///@name Life Cycle
    
    ///@{

    //* Constructor.
    
    ResidualCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRatioTolerance       = NewRatioTolerance;
        mAlwaysConvergedNorm  = AlwaysConvergedNorm;
        mInitialResidualIsSet = false;
    }

    //* Copy constructor.
    
    ResidualCriteria( ResidualCriteria const& rOther )
      :BaseType(rOther) 
      ,mInitialResidualIsSet(rOther.mInitialResidualIsSet)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mInitialResidualNorm(rOther.mInitialResidualNorm)
      ,mCurrentResidualNorm(rOther.mCurrentResidualNorm)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
    {
    }

    //* Destructor.
    
    ~ResidualCriteria() override {}


    ///@} 
    ///@name Operators
    
    ///@{

    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        const std::size_t size_b = b.size();
        if (size_b != 0) { //if we are solving for something

            std::size_t size_residual;
            if (mInitialResidualIsSet == false) {
//                 mInitialResidualNorm = TSparseSpace::TwoNorm(b); // NOTE: This doesn't take into account the reaction dofs
                GetResidualNorm(mInitialResidualNorm, size_residual, rDofSet, b);
                mInitialResidualIsSet = true;

                //KRATOS_INFO("RESIDUAL CRITERION") << "Initial Residual: " << mInitialResidualNorm <<std::endl;
            }

            TDataType ratio;
//             mCurrentResidualNorm = TSparseSpace::TwoNorm(b); // NOTE: This doesn't take into account the reaction dofs
            GetResidualNorm(mInitialResidualNorm, size_residual, rDofSet, b);

            if(mInitialResidualNorm < std::numeric_limits<TDataType>::epsilon()) {
                ratio = 0.0;
            } else {
                ratio = mCurrentResidualNorm/mInitialResidualNorm;
            }

            //KRATOS_INFO("RESIDUAL CRITERION") << "Current Residual: " <<  mCurrentResidualNorm << " ratio: "<< ratio << std::endl;
            
            const TDataType float_size_residual = static_cast<TDataType>(size_residual);
//             const TDataType float_size_residual = static_cast<TDataType>(size_b);
            const TDataType absolute_norm = (mCurrentResidualNorm/float_size_residual);

            KRATOS_INFO_IF("RESIDUAL CRITERION", this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0) << " :: Ratio = "<< ratio  << ";  Norm = " << absolute_norm << std::endl;

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

            if (ratio <= mRatioTolerance || absolute_norm < mAlwaysConvergedNorm) {
                KRATOS_INFO_IF("RESIDUAL CRITERION", this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0) << "Convergence is achieved" << std::endl;
                return true;
            } else {
                return false;
            }
        } else {
            return true;
        }
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the problem. (unused)
     */
    void Initialize(
        ModelPart& rModelPart
        ) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        mInitialResidualIsSet = false;
    }

    /**
     * This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rDofSet, A, Dx, b);
    }

    ///@} 
    ///@name Operations 
    ///@{


    ///@} 
    ///@name Access 
    ///@{


    ///@} 
    ///@name Inquiry 
    ///@{


    ///@} 
    ///@name Friends 
    ///@{


    ///@} 

protected:
    ///@name Protected static Member Variables 
    ///@{


    ///@} 
    ///@name Protected member Variables 
    ///@{


    ///@} 
    ///@name Protected Operators
    ///@{


    ///@} 
    ///@name Protected Operations
    ///@{


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


    bool mInitialResidualIsSet; /// This "flag" is set in order to set that the initial residual is already computed

    TDataType mRatioTolerance; /// The ratio threshold for the norm of the residual

    TDataType mInitialResidualNorm; /// The reference norm of the residual

    TDataType mCurrentResidualNorm; /// The current norm of the residual

    TDataType mAlwaysConvergedNorm; /// The absolute value threshold for the norm of the residual

    TDataType mReferenceDispNorm; /// The ratio threshold for the norm of the residual


    ///@} 
    ///@name Private Operators
    ///@{

    ///@} 
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the norm of the residual
     * @details It checks if the dof is fixed
     * @param rResidualSolutionNorm The norm of the residual
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param b RHS vector (residual + reactions)
     */
    void GetResidualNorm(
        TDataType& rResidualSolutionNorm,
        IndexType& rDofNum,
        DofsArrayType& rDofSet,
        const TSystemVectorType& b
        )
    {
        // Initialize
        rResidualSolutionNorm = 0.0;
        rDofNum = 0;

        // Loop over Dofs
        #pragma omp parallel for reduction(+:rResidualSolutionNorm,rDofNum)
        for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
            auto it_dof = rDofSet.begin() + i;

            IndexType dof_id;
            TDataType residual_dof_value;

            if (it_dof->IsFree()) {
                dof_id = it_dof->EquationId();
                residual_dof_value = b[dof_id];

                const auto curr_var = it_dof->GetVariable();
                rResidualSolutionNorm += residual_dof_value * residual_dof_value;
                rDofNum++;
            }
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

}; // Class ClassName 

///@} 

///@name Type Definitions 
///@{


///@} 

}  // namespace Kratos.

#endif // KRATOS_NEW_DISPLACEMENT_CRITERIA  defined 

