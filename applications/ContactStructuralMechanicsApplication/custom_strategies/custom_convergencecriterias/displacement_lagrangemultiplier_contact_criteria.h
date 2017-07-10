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

#if !defined(KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_CONTACT_CRITERIA_H)
#define	KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_CONTACT_CRITERIA_H

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
class DisplacementLagrangeMultiplierContactCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementLagrangeMultiplierContactCriteria );

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
     * @param DispRatioTolerance: Relative tolerance for displacement error
     * @param DispAbsTolerance: Absolute tolerance for displacement error
     * @param LMRatioTolerance: Relative tolerance for lagrange multiplier error
     * @param LMAbsTolerance: Absolute tolerance for lagrange multiplier error
     * @param EnsureContact: To check if the contact is lost
     */
    
    DisplacementLagrangeMultiplierContactCriteria(  
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
    }

    //* Copy constructor.
    
    DisplacementLagrangeMultiplierContactCriteria( DisplacementLagrangeMultiplierContactCriteria const& rOther )
      :BaseType(rOther) 
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mLMRatioTolerance(rOther.mLMRatioTolerance)
      ,mLMAbsTolerance(rOther.mLMAbsTolerance)
    {
    }
    
    /// Destructor.
    ~DisplacementLagrangeMultiplierContactCriteria() override {}

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
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            // Initialize
            TDataType DispSolutionNorm = 0.0;
            TDataType LMSolutionNorm = 0.0;
            TDataType DispIncreaseNorm = 0.0;
            TDataType LMIncreaseNorm = 0.0;
            unsigned int DispDofNum(0),LMDofNum(0);

            // Set a partition for OpenMP
            int NumDofs = rDofSet.size();
            PartitionVector DofPartition;
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(NumDofs,NumThreads,DofPartition);

            // Loop over Dofs
            #pragma omp parallel reduction(+:DispSolutionNorm,LMSolutionNorm,DispIncreaseNorm,LMIncreaseNorm,DispDofNum,LMDofNum)
            {
                int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
                typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k+1];

                std::size_t DofId;
                TDataType DofValue;
                TDataType DofIncr;

                for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
                {
                    if (itDof->IsFree())
                    {
                        DofId = itDof->EquationId();
                        DofValue = itDof->GetSolutionStepValue(0);
                        DofIncr = Dx[DofId];

                        KeyType CurrVar = itDof->GetVariable().Key();
                        if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z))
                        {
                            DispSolutionNorm += DofValue * DofValue;
                            DispIncreaseNorm += DofIncr * DofIncr;
                            ++DispDofNum;
                        }
                        else
                        {
                            LMSolutionNorm += DofValue * DofValue;
                            LMIncreaseNorm += DofIncr * DofIncr;
                            ++LMDofNum;
                        }
                    }
                }
            }

            if(DispSolutionNorm == 0.0)
            {
                DispSolutionNorm = 1.0;
            }
            if(LMSolutionNorm == 0.0)
            {
                LMSolutionNorm = 1.0;
                if (mEnsureContact == true)
                {
                    KRATOS_ERROR << "WARNING::CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;
                }
            }
            
            TDataType DispRatio = std::sqrt(DispIncreaseNorm/DispSolutionNorm);
            TDataType LMRatio = std::sqrt(LMIncreaseNorm/LMSolutionNorm);

            TDataType DispAbs = std::sqrt(DispIncreaseNorm)/ static_cast<TDataType>(DispDofNum);
            TDataType LMAbs = std::sqrt(LMIncreaseNorm)/ static_cast<TDataType>(LMDofNum);

            // We print the results
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                if (mpTable != nullptr)
                {
                    std::cout.precision(4);
                    auto& Table = mpTable->GetTable();
                    Table  << DispRatio  << mDispRatioTolerance  << DispAbs  << mDispAbsTolerance  << LMRatio  << mLMRatioTolerance  << LMAbs  << mLMAbsTolerance;
                }
                else
                {
                    std::cout.precision(4);
                    #if !defined(_WIN32)
                        std::cout << BOLDFONT("DoF ONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                        std::cout << BOLDFONT("\tDISPLACEMENT: RATIO = ") << DispRatio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << DispAbs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        std::cout << BOLDFONT(" LAGRANGE MUL:\tRATIO = ") << LMRatio << BOLDFONT(" EXP.RATIO = ") << mLMRatioTolerance << BOLDFONT(" ABS = ") << LMAbs << BOLDFONT(" EXP.ABS = ") << mLMAbsTolerance << std::endl;
                    #else
                        std::cout << "DoF ONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                        std::cout << "\tDISPLACEMENT: RATIO = " << DispRatio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << DispAbs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        std::cout << " LAGRANGE MUL:\tRATIO = " << LMRatio << " EXP.RATIO = " << mLMRatioTolerance << " ABS = " << LMAbs << " EXP.ABS = " << mLMAbsTolerance << std::endl;
                    #endif
                }
            }

            if ((DispRatio <= mDispRatioTolerance || DispAbs <= mDispAbsTolerance) &&
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
                            std::cout << BOLDFONT("\tDoF") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        #else
                            std::cout << "\tDoF convergence is achieved" << std::endl;
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
                            std::cout << BOLDFONT("\tDoF") << " convergence is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        #else
                            std::cout << "\tDoF convergence is not achieved" << std::endl;
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
    
    void Initialize( ModelPart& rModelPart ) override
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
    
    const bool mEnsureContact;
    
    TablePrinterPointerType mpTable; // Pointer to the fancy table 
    bool mTableIsInitialized;        // If the table is already initialized
    
    TDataType mDispRatioTolerance;
    TDataType mDispAbsTolerance;

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

#endif	/* KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_CONTACT_CRITERIA_H */

