// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_MIXED_CONTACT_CRITERIA_H)
#define	KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_MIXED_CONTACT_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "utilities/table_stream_utility.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/color_utilities.h"

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
class DisplacementLagrangeMultiplierMixedContactCriteria 
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementLagrangeMultiplierMixedContactCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >     BaseType;

    typedef TSparseSpace                                  SparseSpaceType;

    typedef typename BaseType::TDataType                        TDataType;

    typedef typename BaseType::DofsArrayType                DofsArrayType;

    typedef typename BaseType::TSystemMatrixType        TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType        TSystemVectorType;

    typedef OpenMPUtils::PartitionVector                  PartitionVector;

    typedef std::size_t                                           KeyType;
    
    typedef TableStreamUtility::Pointer           TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param DispRatioTolerance Relative tolerance for displacement residual error
     * @param DispAbsTolerance Absolute tolerance for displacement residual error
     * @param LMRatioTolerance Relative tolerance for lagrange multiplier residual  error
     * @param LMAbsTolerance Absolute tolerance for lagrange multiplier residual error
     * @param EnsureContact To check if the contact is lost
     * @param pTable The pointer to the output table
     * @param PrintingOutput If the output is going to be printed in a txt file
     */
    
    DisplacementLagrangeMultiplierMixedContactCriteria(  
        TDataType DispRatioTolerance,
        TDataType DispAbsTolerance,
        TDataType LMRatioTolerance,
        TDataType LMAbsTolerance,
        bool EnsureContact = false,
        TablePrinterPointerType pTable = nullptr,
        const bool PrintingOutput = false 
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mEnsureContact(EnsureContact),
          mpTable(pTable),
          mPrintingOutput(PrintingOutput),
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
      ,mpTable(rOther.mpTable)
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
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
            TDataType disp_residual_solution_norm = 0.0;
            TDataType lm_solution_norm = 0.0;
            TDataType lm_increase_norm = 0.0;
            unsigned int disp_dof_num(0),lm_dof_num(0);

            // Set a partition for OpenMP
            const int num_dofs = rDofSet.size();
            PartitionVector dof_partition;
            const int num_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(num_dofs,num_threads,dof_partition);

            // Loop over Dofs
            #pragma omp parallel reduction(+:disp_residual_solution_norm,lm_solution_norm,lm_increase_norm,disp_dof_num,lm_dof_num)
            {
                const int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator dof_begin = rDofSet.begin() + dof_partition[k];
                typename DofsArrayType::iterator dof_end = rDofSet.begin() + dof_partition[k+1];

                std::size_t dof_id;
                TDataType residual_dof_value;
                TDataType dof_value;
                TDataType dof_incr;

                for (typename DofsArrayType::iterator it_dof = dof_begin; it_dof != dof_end; ++it_dof)
                {
                    if (it_dof->IsFree())
                    {
                        dof_id = it_dof->EquationId();
                        
                        const KeyType curr_var = it_dof->GetVariable().Key();
                        if ((curr_var == DISPLACEMENT_X) || (curr_var == DISPLACEMENT_Y) || (curr_var == DISPLACEMENT_Z))
                        {
                            residual_dof_value = b[dof_id];
                            disp_residual_solution_norm += residual_dof_value * residual_dof_value;
                            ++disp_dof_num;
                        }
                        else
                        {
                            dof_value = it_dof->GetSolutionStepValue(0);
                            dof_incr = Dx[dof_id];
                            lm_solution_norm += dof_value * dof_value;
                            lm_increase_norm += dof_incr * dof_incr;
                            ++lm_dof_num;
                        }
                    }
                }
            }

            if(lm_increase_norm == 0.0) lm_increase_norm = 1.0;
            KRATOS_ERROR_IF(mEnsureContact == true && lm_solution_norm == 0.0) << "WARNING::CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;
            
            mDispCurrentResidualNorm = disp_residual_solution_norm;
            TDataType lm_ratio = std::sqrt(lm_increase_norm/lm_solution_norm);
            TDataType lm_abs = std::sqrt(lm_increase_norm)/ static_cast<TDataType>(lm_dof_num);
            
            TDataType residual_disp_ratio; 
            
            // We initialize the solution
            if (mInitialResidualIsSet == false)
            {
                mDispInitialResidualNorm = (disp_residual_solution_norm == 0.0) ? 1.0 : disp_residual_solution_norm;
                residual_disp_ratio = 1.0;
                mInitialResidualIsSet = true;
            }
            
            // We calculate the ratio of the displacements
            residual_disp_ratio = mDispCurrentResidualNorm/mDispInitialResidualNorm;

            // We calculate the absolute norms
            TDataType residual_disp_abs = mDispCurrentResidualNorm/disp_dof_num;

            // We print the results
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                if (mpTable != nullptr)
                {
                    std::cout.precision(4);
                    auto& Table = mpTable->GetTable();
                    Table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << lm_ratio << mLMRatioTolerance << lm_abs << mLMAbsTolerance;
                }
                else
                {
                    std::cout.precision(4);
                    if (mPrintingOutput == false)
                    {
                        std::cout << BOLDFONT("MIXED CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[STEP] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        std::cout << BOLDFONT("\tDISPLACEMENT: RATIO = ") << residual_disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << residual_disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        std::cout << BOLDFONT("\tLAGRANGE MUL: RATIO = ") << lm_ratio << BOLDFONT(" EXP.RATIO = ") << mLMRatioTolerance << BOLDFONT(" ABS = ") << lm_abs << BOLDFONT(" EXP.ABS = ") << mLMAbsTolerance << std::endl;
                    }
                    else
                    {
                        std::cout << "MIXED CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[STEP] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        std::cout << "\tDISPLACEMENT: RATIO = " << residual_disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << residual_disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        std::cout << "\tLAGRANGE MUL: RATIO = " << lm_ratio << " EXP.RATIO = " << mLMRatioTolerance << " ABS = " << lm_abs << " EXP.ABS = " << mLMAbsTolerance << std::endl;
                    }
                }
            }

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = (residual_disp_ratio > lm_ratio) ? residual_disp_ratio : lm_ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = (lm_abs > mLMAbsTolerance) ? lm_abs : mLMAbsTolerance;
            
            if ((residual_disp_ratio <= mDispRatioTolerance || residual_disp_abs <= mDispAbsTolerance) &&
                    (lm_ratio <= mLMRatioTolerance || lm_abs <= mLMAbsTolerance) )
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    if (mpTable != nullptr)
                    {
                        auto& table = mpTable->GetTable();
                        if (mPrintingOutput == false)
                            table << BOLDFONT(FGRN("       Achieved"));
                        else
                            table << "Achieved";
                    }
                    else
                    {
                        if (mPrintingOutput == false)
                            std::cout << BOLDFONT("\tConvergence") << " is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        else
                            std::cout << "\tConvergence is achieved" << std::endl;
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
                        auto& table = mpTable->GetTable();
                        if (mPrintingOutput == false)
                            table << BOLDFONT(FRED("   Not achieved"));
                        else
                            table << "Not achieved";
                    }
                    else
                    {
                        if (mPrintingOutput == false)
                            std::cout << BOLDFONT("\tConvergence") << " is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        else
                            std::cout << "\tConvergence is not achieved" << std::endl;
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
            auto& table = mpTable->GetTable();
            table.AddColumn("DP RATIO", 10);
            table.AddColumn("EXP. RAT", 10);
            table.AddColumn("ABS", 10);
            table.AddColumn("EXP. ABS", 10);
            table.AddColumn("LM RATIO", 10);
            table.AddColumn("EXP. RAT", 10);
            table.AddColumn("ABS", 10);
            table.AddColumn("EXP. ABS", 10);
            table.AddColumn("CONVERGENCE", 15);
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
    bool mPrintingOutput;            // If the colors and bold are printed
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

