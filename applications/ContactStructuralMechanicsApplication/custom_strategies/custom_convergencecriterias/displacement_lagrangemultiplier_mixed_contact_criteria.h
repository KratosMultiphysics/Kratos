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

#if !defined(KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_MIXED_CONTACT_CRITERIA_H)
#define	KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_MIXED_CONTACT_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "custom_utilities/bprinter_utility.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#if !defined(_WIN32)
	#include "custom_utilities/color_utilities.h"
#endif

namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
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
class DisplacementLagrangeMultiplierMixedContactCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementLagrangeMultiplierMixedContactCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >  BaseType;

    typedef TSparseSpace                               SparseSpaceType;

    typedef typename BaseType::TDataType                     TDataType;

    typedef typename BaseType::DofsArrayType             DofsArrayType;

    typedef typename BaseType::TSystemMatrixType     TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType     TSystemVectorType;

    typedef OpenMPUtils::PartitionVector               PartitionVector;

    typedef std::size_t                                        KeyType;
    
    typedef boost::shared_ptr<BprinterUtility> TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param DispRatioTolerance Relative tolerance for displacement residual error
     * @param DispAbsTolerance Absolute tolerance for displacement residual error
     * @param LMRatioTolerance Relative tolerance for lagrange multiplier residual  error
     * @param LMAbsTolerance Absolute tolerance for lagrange multiplier residual error
     * @param EnsureContact: To check if the contact is lost
     */
    
    DisplacementLagrangeMultiplierMixedContactCriteria(  
        TDataType DispRatioTolerance,
        TDataType DispAbsTolerance,
        TDataType LMRatioTolerance,
        TDataType LMAbsTolerance,
        bool EnsureContact = false,
        TablePrinterPointerType pTable = nullptr
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mEnsureContact(EnsureContact),
          mpTable(pTable),
          mTableIsInitialized(false)
    {
        mDispRatioTolerance = DispRatioTolerance;
        mDispAbsTolerance = DispAbsTolerance;

        mLMRatioTolerance = LMRatioTolerance;
        mLMAbsTolerance = LMAbsTolerance;
        
        mInitialResidualIsSet = false;
    }

    //* Copy constructor.
    DisplacementLagrangeMultiplierMixedContactCriteria( DisplacementLagrangeMultiplierMixedContactCriteria const& rOther )
      :BaseType(rOther) 
      ,mInitialResidualIsSet(rOther.mInitialResidualIsSet)
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mDispInitialResidualNorm(rOther.mDispInitialResidualNorm)
      ,mDispCurrentResidualNorm(rOther.mDispCurrentResidualNorm)
      ,mLMRatioTolerance(rOther.mLMRatioTolerance)
      ,mLMAbsTolerance(rOther.mLMAbsTolerance)
    {
    }
    
    /// Destructor.
    ~DisplacementLagrangeMultiplierMixedContactCriteria() override = default;

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
            TDataType LMSolutionNorm = 0.0;
            TDataType LMIncreaseNorm = 0.0;
            unsigned int DispDofNum(0),LMDofNum(0);

            // Set a partition for OpenMP
            int NumDofs = rDofSet.size();
            PartitionVector DofPartition;
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(NumDofs,NumThreads,DofPartition);

            // Loop over Dofs
            #pragma omp parallel reduction(+:DispResidualSolutionNorm,LMSolutionNorm,LMIncreaseNorm,DispDofNum,LMDofNum)
            {
                int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
                typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k+1];

                std::size_t DofId;
                TDataType ResidualDofValue;
                TDataType DofValue;
                TDataType DofIncr;

                for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
                {
                    if (itDof->IsFree())
                    {
                        DofId = itDof->EquationId();
                        
                        KeyType CurrVar = itDof->GetVariable().Key();
                        if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z))
                        {
                            ResidualDofValue = b[DofId];
                            DispResidualSolutionNorm += ResidualDofValue * ResidualDofValue;
                            ++DispDofNum;
                        }
                        else
                        {
                            DofValue = itDof->GetSolutionStepValue(0);
                            DofIncr = Dx[DofId];
                            LMSolutionNorm += DofValue * DofValue;
                            LMIncreaseNorm += DofIncr * DofIncr;
                            ++LMDofNum;
                        }
                    }
                }
            }

            if(LMIncreaseNorm == 0.0) LMIncreaseNorm = 1.0;
            if(LMSolutionNorm == 0.0)
            {
                LMSolutionNorm = 1.0;
                if (mEnsureContact == true)
                {
                    KRATOS_ERROR << "WARNING::CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;
                }
            }
            
            mDispCurrentResidualNorm = DispResidualSolutionNorm;
            TDataType LMRatio = std::sqrt(LMIncreaseNorm/LMSolutionNorm);
            TDataType LMAbs = std::sqrt(LMIncreaseNorm)/ static_cast<TDataType>(LMDofNum);
            
            TDataType ResidualDispRatio; 
            
            // We initialize the solution
            if (mInitialResidualIsSet == false)
            {
                if (DispResidualSolutionNorm == 0.0)
                {
                    mDispInitialResidualNorm = 1.0;
                }
                else
                {
                    mDispInitialResidualNorm = DispResidualSolutionNorm;
                }
                
                ResidualDispRatio = 1.0;
                mInitialResidualIsSet = true;
            }
            
            // We calculate the ratio of the displacements
            ResidualDispRatio = mDispCurrentResidualNorm/mDispInitialResidualNorm;

            // We calculate the absolute norms
            TDataType ResidualDispAbs = mDispCurrentResidualNorm/DispDofNum;

            // We print the results
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                if (mpTable != nullptr)
                {
                    std::cout.precision(4);
                    auto& Table = mpTable->GetTable();
                    Table << ResidualDispRatio << mDispRatioTolerance << ResidualDispAbs << mDispAbsTolerance << LMRatio << mLMRatioTolerance << LMAbs << mLMAbsTolerance;
                }
                else
                {
                    std::cout.precision(4);
                    #if !defined(_WIN32)
                        std::cout << BOLDFONT("MIXED CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        std::cout << BOLDFONT("\tDISPLACEMENT: RATIO = ") << ResidualDispRatio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << ResidualDispAbs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        std::cout << BOLDFONT("\tLAGRANGE MUL: RATIO = ") << LMRatio << BOLDFONT(" EXP.RATIO = ") << mLMRatioTolerance << BOLDFONT(" ABS = ") << LMAbs << BOLDFONT(" EXP.ABS = ") << mLMAbsTolerance << std::endl;
                    #else
                        std::cout << "MIXED CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        std::cout << "\tDISPLACEMENT: RATIO = " << ResidualDispRatio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << ResidualDispAbs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        std::cout << "\tLAGRANGE MUL: RATIO = " << LMRatio << " EXP.RATIO = " << mLMRatioTolerance << " ABS = " << LMAbs << " EXP.ABS = " << mLMAbsTolerance << std::endl;
                    #endif
                }
            }

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = (ResidualDispRatio > LMRatio) ? ResidualDispRatio : LMRatio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = (LMAbs > mLMAbsTolerance) ? LMAbs : mLMAbsTolerance;
            
            if ((ResidualDispRatio <= mDispRatioTolerance || ResidualDispAbs <= mDispAbsTolerance) &&
                    (LMRatio <= mLMRatioTolerance || LMAbs <= mLMAbsTolerance) )
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    if (mpTable != nullptr)
                    {
                        auto& Table = mpTable->GetTable();
                        #if !defined(_WIN32)
                            Table << BOLDFONT(FGRN("       Achieved"));
                        #else
                            Table << "Achieved";
                        #endif
                    }
                    else
                    {
                        #if !defined(_WIN32)
                            std::cout << BOLDFONT("\tConvergence") << " is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        #else
                            std::cout << "\tConvergence is achieved" << std::endl;
                        #endif
                    }
                }
                return true;
            }
            else
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    if (mpTable != nullptr)
                    {
                        auto& Table = mpTable->GetTable();
                        #if !defined(_WIN32)
                            Table << BOLDFONT(FRED("   Not achieved"));
                        #else
                            Table << "Not achieved";
                        #endif
                    }
                    else
                    {
                        #if !defined(_WIN32)
                            std::cout << BOLDFONT("\tConvergence") << " is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        #else
                            std::cout << "\tConvergence is not achieved" << std::endl;
                        #endif
                    }
                }
                return false;
            }
        }
        else // In this case all the displacements are imposed!
        {
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
        
        if (mpTable != nullptr && mTableIsInitialized == false)
        {
            auto& Table = mpTable->GetTable();
            Table.AddColumn("DP RATIO", 10);
            Table.AddColumn("EXP. RAT", 10);
            Table.AddColumn("ABS", 10);
            Table.AddColumn("EXP. ABS", 10);
            Table.AddColumn("LM RATIO", 10);
            Table.AddColumn("EXP. RAT", 10);
            Table.AddColumn("ABS", 10);
            Table.AddColumn("EXP. ABS", 10);
            Table.AddColumn("CONVERGENCE", 15);
            mTableIsInitialized = true;
        }
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
    
    const bool mEnsureContact;
    
    TablePrinterPointerType mpTable; // Pointer to the fancy table 
    bool mTableIsInitialized;        // If the table is already initialized
    
    TDataType mDispRatioTolerance;
    TDataType mDispAbsTolerance;
    TDataType mDispInitialResidualNorm;
    TDataType mDispCurrentResidualNorm;
    
    TDataType mLMRatioTolerance;
    TDataType mLMAbsTolerance;
    
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

#endif	/* KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_MIXED_CONTACT_CRITERIA_H */

