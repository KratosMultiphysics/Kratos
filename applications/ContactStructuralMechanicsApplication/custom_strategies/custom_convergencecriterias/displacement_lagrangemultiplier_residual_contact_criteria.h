// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_CONTACT_CRITERIA_H)
#define	KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_CONTACT_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_utilities/color_utilities.h"

namespace Kratos
{
///@addtogroup ContacStructuralMechanicsApplication
///@{

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
    
///@name Kratos Classes
///@{

/**
 * @brief Convergence criteria for contact problems
 * This class implements a convergence control based on nodal displacement and
 * lagrange multiplier values. The error is evaluated separately for each of them, and
 * relative and absolute tolerances for both must be specified.
 */
template<   class TSparseSpace,
            class TDenseSpace >
class DisplacementLagrangeMultiplierResidualContactCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementLagrangeMultiplierResidualContactCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef OpenMPUtils::PartitionVector PartitionVector;

    typedef std::size_t KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param DispRatioTolerance Relative tolerance for displacement residual error
     * @param DispAbsTolerance Absolute tolerance for displacement residual error
     * @param LMRatioTolerance Relative tolerance for lagrange multiplier residual  error
     * @param LMAbsTolerance Absolute tolerance for lagrange multiplier residual error
     */
    
    DisplacementLagrangeMultiplierResidualContactCriteria(  
                    TDataType DispRatioTolerance,
                    TDataType DispAbsTolerance,
                    TDataType LMRatioTolerance,
                    TDataType LMAbsTolerance)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mDispRatioTolerance = DispRatioTolerance;
        mDispAbsTolerance = DispAbsTolerance;

        mLMRatioTolerance = LMRatioTolerance;
        mLMAbsTolerance = LMAbsTolerance;
        
        mInitialResidualIsSet = false;
    }

    //* Copy constructor.
    DisplacementLagrangeMultiplierResidualContactCriteria( DisplacementLagrangeMultiplierResidualContactCriteria const& rOther )
      :BaseType(rOther) 
      ,mInitialResidualIsSet(rOther.mInitialResidualIsSet)
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mDispInitialResidualNorm(rOther.mDispInitialResidualNorm)
      ,mDispCurrentResidualNorm(rOther.mDispCurrentResidualNorm)
      ,mLMRatioTolerance(rOther.mLMRatioTolerance)
      ,mLMAbsTolerance(rOther.mLMAbsTolerance)
      ,mLMInitialResidualNorm(rOther.mLMInitialResidualNorm)
      ,mLMCurrentResidualNorm(rOther.mLMCurrentResidualNorm)
    {
    }
    
    /// Destructor.
    virtual ~DisplacementLagrangeMultiplierResidualContactCriteria() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
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
        const TSystemVectorType& b 
        ) override
    {
        if (SparseSpaceType::Size(b) != 0) //if we are solving for something
        {
            // Initialize
            TDataType DispResidualSolutionNorm = 0.0;
            TDataType LMResidualSolutionNorm = 0.0;
            unsigned int DispDofNum(0),LMDofNum(0);

            // Set a partition for OpenMP
            int NumDofs = rDofSet.size();
            PartitionVector DofPartition;
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(NumDofs,NumThreads,DofPartition);

            // Loop over Dofs
            #pragma omp parallel reduction(+:DispResidualSolutionNorm,LMResidualSolutionNorm,DispDofNum,LMDofNum)
            {
                int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
                typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k+1];

                std::size_t DofId;
                TDataType ResidualDofValue;

                for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
                {
                    if (itDof->IsFree())
                    {
                        DofId = itDof->EquationId();
                        ResidualDofValue = b[DofId];

                        KeyType CurrVar = itDof->GetVariable().Key();
                        if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z))
                        {
                            DispResidualSolutionNorm += ResidualDofValue * ResidualDofValue;
                            ++DispDofNum;
                        }
                        else
                        {
                            LMResidualSolutionNorm += ResidualDofValue * ResidualDofValue;
                            ++LMDofNum;
                        }
                    }
                }
            }

            mDispCurrentResidualNorm = DispResidualSolutionNorm;
            mLMCurrentResidualNorm = LMResidualSolutionNorm;
            
            TDataType ResidualDispRatio; 
            TDataType ResidualLMRatio;
            
            // We initialize the solution
            if (mInitialResidualIsSet == false)
            {
                mDispInitialResidualNorm = DispResidualSolutionNorm;
                mLMInitialResidualNorm = LMResidualSolutionNorm;
                ResidualDispRatio = 1.0;
                ResidualLMRatio = 1.0;
                mInitialResidualIsSet = true;
            }
            
            // We calculate the ratio of the displacements
            if(mDispInitialResidualNorm == 0.00)
            {
                ResidualDispRatio = 0.00;
            }
            else
            {
                ResidualDispRatio = mDispCurrentResidualNorm/mDispInitialResidualNorm;
            }
            
            // We calculate the ratio of the LM
            if(mLMInitialResidualNorm == 0.00)
            {
                ResidualLMRatio = 0.00;
            }
            else
            {
                ResidualLMRatio = mLMCurrentResidualNorm/mLMInitialResidualNorm;
            }

            // We calculate the absolute norms
            TDataType ResidualDispAbs = mDispCurrentResidualNorm/DispDofNum;
            TDataType ResidualLMAbs = mLMCurrentResidualNorm/LMDofNum;

            // We print the results
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout.precision(4);
                std::cout << BOLD("RESIDUAL CONVERGENCE CHECK:") << std::endl << std::scientific;
                std::cout << BOLD("\tDISPLACEMENT: RATIO = ") << ResidualDispRatio << BOLD(" EXP.RATIO = ") << mDispRatioTolerance << BOLD(" ABS = ") << ResidualDispAbs  << BOLD(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                std::cout << BOLD("\tLAGRANGE MUL: RATIO = ") << ResidualLMRatio  << BOLD(" EXP.RATIO = ") << mLMRatioTolerance << BOLD(" ABS = ") << ResidualLMAbs << BOLD(" EXP.ABS = ") << mLMAbsTolerance << std::endl;
            }

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = (ResidualDispRatio > ResidualLMRatio) ? ResidualDispRatio : ResidualLMRatio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = (ResidualLMAbs > mLMAbsTolerance) ? ResidualLMAbs : mLMAbsTolerance;
            
            if ((ResidualDispRatio <= mDispRatioTolerance || ResidualDispAbs <= mDispAbsTolerance) &&
                    (ResidualLMRatio <= mLMRatioTolerance || ResidualLMAbs <= mLMAbsTolerance) )
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    std::cout << BOLD("\tResidual") << " convergence is " << BOLD(FGRN("achieved")) << std::endl;
                }
                return true;
            }
            else
            {
                return false;
            }
        }
        else // In this case all the displacements are imposed!
        {
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << BOLD("\tResidual") << " convergence is " << BOLD(FRED(" not achieved")) << std::endl;
            }
            return true;
        }
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the contact problem. (unused)
     */
    
    void Initialize( ModelPart& rModelPart) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
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
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
        
    void FinalizeSolutionStep(  
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b 
        ) override
    {
        
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

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
    
    bool mInitialResidualIsSet;
    
    TDataType mDispRatioTolerance;
    TDataType mDispAbsTolerance;
    TDataType mDispInitialResidualNorm;
    TDataType mDispCurrentResidualNorm;
    
    TDataType mLMRatioTolerance;
    TDataType mLMAbsTolerance;
    TDataType mLMInitialResidualNorm;
    TDataType mLMCurrentResidualNorm;
    
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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};

///@} // Kratos classes

///@} // Application group
}

#endif	/* KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_CONTACT_CRITERIA_H */

