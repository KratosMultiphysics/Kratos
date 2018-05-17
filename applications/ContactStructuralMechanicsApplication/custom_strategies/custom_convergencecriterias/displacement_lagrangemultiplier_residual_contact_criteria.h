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

#if !defined(KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_CONTACT_CRITERIA_H)
#define KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_CONTACT_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
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
 * @class DisplacementLagrangeMultiplierResidualContactCriteria 
 * @ingroup ContactStructuralMechanicsApplication 
 * @brief Convergence criteria for contact problems
 * This class implements a convergence control based on nodal displacement and
 * lagrange multiplier values. The error is evaluated separately for each of them, and
 * relative and absolute tolerances for both must be specified.
 * @author Vicente Mataix Ferrandiz
 */
template<   class TSparseSpace,
            class TDenseSpace >
class DisplacementLagrangeMultiplierResidualContactCriteria 
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementLagrangeMultiplierResidualContactCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >     BaseType;

    typedef TSparseSpace                                  SparseSpaceType;

    typedef typename BaseType::TDataType                        TDataType;

    typedef typename BaseType::DofsArrayType                DofsArrayType;

    typedef typename BaseType::TSystemMatrixType        TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType        TSystemVectorType;

    typedef OpenMPUtils::PartitionVector                  PartitionVector;

    typedef std::size_t                                           KeyType;
    
    typedef TableStreamUtility::Pointer           TablePrinterPointerType;
    
    typedef std::size_t                                         IndexType;

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
    
    DisplacementLagrangeMultiplierResidualContactCriteria(  
        TDataType DispRatioTolerance,
        TDataType DispAbsTolerance,
        TDataType LMRatioTolerance,
        TDataType LMAbsTolerance,
        bool EnsureContact = false,
        const bool PrintingOutput = false
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mEnsureContact(EnsureContact),
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
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
    {
    }
    
    /// Destructor.
    ~DisplacementLagrangeMultiplierResidualContactCriteria() override = default;

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
        if (SparseSpaceType::Size(b) != 0) { //if we are solving for something
            // Initialize
            TDataType disp_residual_solution_norm = 0.0, lm_residual_solution_norm = 0.0;
            IndexType disp_dof_num(0),lm_dof_num(0);

            // Loop over Dofs
            #pragma omp parallel for reduction(+:disp_residual_solution_norm,lm_residual_solution_norm,disp_dof_num,lm_dof_num)
            for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
                auto it_dof = rDofSet.begin() + i;

                std::size_t dof_id;
                TDataType residual_dof_value;
                
                if (it_dof->IsFree()) {
                    dof_id = it_dof->EquationId();
                    residual_dof_value = b[dof_id];

                    const auto curr_var = it_dof->GetVariable();
                    if ((curr_var == VECTOR_LAGRANGE_MULTIPLIER_X) || (curr_var == VECTOR_LAGRANGE_MULTIPLIER_Y) || (curr_var == VECTOR_LAGRANGE_MULTIPLIER_Z) || (curr_var == LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)) {
                        lm_residual_solution_norm += residual_dof_value * residual_dof_value;
                        lm_dof_num++;
                    } else {
                        disp_residual_solution_norm += residual_dof_value * residual_dof_value;
                        disp_dof_num++;
                    }
                }
            }

            mDispCurrentResidualNorm = disp_residual_solution_norm;
            mLMCurrentResidualNorm = lm_residual_solution_norm;
            
            TDataType residual_disp_ratio = 1.0;
            TDataType residual_lm_ratio = 1.0;
            
            // We initialize the solution
            if (mInitialResidualIsSet == false) {
                mDispInitialResidualNorm = (disp_residual_solution_norm == 0.0) ? 1.0 : disp_residual_solution_norm;
                mLMInitialResidualNorm = (lm_residual_solution_norm == 0.0) ? 1.0 : lm_residual_solution_norm;
                residual_disp_ratio = 1.0;
                residual_lm_ratio = 1.0;
                mInitialResidualIsSet = true;
            }
            
            // We calculate the ratio of the displacements
            residual_disp_ratio = mDispCurrentResidualNorm/mDispInitialResidualNorm;
            
            // We calculate the ratio of the LM
            residual_lm_ratio = mLMCurrentResidualNorm/mLMInitialResidualNorm;

            KRATOS_ERROR_IF(mEnsureContact && residual_lm_ratio == 0.0) << "ERROR::CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;
            
            // We calculate the absolute norms
            const TDataType residual_disp_abs = mDispCurrentResidualNorm/disp_dof_num;
            const TDataType residual_lm_abs = mLMCurrentResidualNorm/lm_dof_num;

            // The process info of the model part
            ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
            
            // We print the results // TODO: Replace for the new log
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                if (r_process_info.Has(TABLE_UTILITY)) {
                    std::cout.precision(4);
                    TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                    auto& Table = p_table->GetTable();
                    Table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << residual_lm_ratio << mLMRatioTolerance << residual_lm_abs << mLMAbsTolerance;
                } else {
                    std::cout.precision(4);
                    if (mPrintingOutput == false) {
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << BOLDFONT("RESIDUAL CONVERGENCE CHECK") << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << BOLDFONT("\tDISPLACEMENT: RATIO = ") << residual_disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << residual_disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << BOLDFONT("\tLAGRANGE MUL: RATIO = ") << residual_lm_ratio << BOLDFONT(" EXP.RATIO = ") << mLMRatioTolerance << BOLDFONT(" ABS = ") << residual_lm_abs << BOLDFONT(" EXP.ABS = ") << mLMAbsTolerance << std::endl;
                    } else {
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << "RESIDUAL CONVERGENCE CHECK" << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << "\tDISPLACEMENT: RATIO = " << residual_disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << residual_disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << "\tLAGRANGE MUL: RATIO = " << residual_lm_ratio << " EXP.RATIO = " << mLMRatioTolerance << " ABS = " << residual_lm_abs << " EXP.ABS = " << mLMAbsTolerance << std::endl;
                    }
                }
            }

            r_process_info[CONVERGENCE_RATIO] = (residual_disp_ratio > residual_lm_ratio) ? residual_disp_ratio : residual_lm_ratio;
            r_process_info[RESIDUAL_NORM] = (residual_lm_abs > mLMAbsTolerance) ? residual_lm_abs : mLMAbsTolerance;
            
            // We check if converged
            const bool disp_converged = (residual_disp_ratio <= mDispRatioTolerance || residual_disp_abs <= mDispAbsTolerance);
            const bool lm_converged = (!mEnsureContact && residual_lm_ratio == 0.0) ? true : (residual_lm_ratio <= mLMRatioTolerance || residual_lm_abs <= mLMAbsTolerance);
            
            if (disp_converged && lm_converged ) {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& Table = p_table->GetTable();
                        if (mPrintingOutput == false)
                            Table << BOLDFONT(FGRN("       Achieved"));
                        else
                            Table << "Achieved";
                    } else {
                        if (mPrintingOutput == false)
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << "\tResidual convergence is achieved" << std::endl;
                    }
                }
                return true;
            } else {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& table = p_table->GetTable();
                        if (mPrintingOutput == false)
                            table << BOLDFONT(FRED("   Not achieved"));
                        else
                            table << "Not achieved";
                    } else {
                        if (mPrintingOutput == false)
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualContactCriteria") << "\tResidual convergence is not achieved" << std::endl;
                    }
                }
                return false;
            }
        } else // In this case all the displacements are imposed!
            return true;
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the contact problem. (unused)
     */
    
    void Initialize( ModelPart& rModelPart) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
        
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && mTableIsInitialized == false) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& table = p_table->GetTable();
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
    
    bool mInitialResidualIsSet; /// This "flag" is set in order to set that the initial residual is already computed
    
    const bool mEnsureContact; /// This "flag" is used to check that the norm of the LM is always greater than 0 (no contact)
    
    bool mPrintingOutput;      /// If the colors and bold are printed
    bool mTableIsInitialized;  /// If the table is already initialized
    
    TDataType mDispRatioTolerance;      /// The ratio threshold for the norm of the displacement residual
    TDataType mDispAbsTolerance;        /// The absolute value threshold for the norm of the displacement residual 
    TDataType mDispInitialResidualNorm; /// The reference norm of the displacement residual
    TDataType mDispCurrentResidualNorm; /// The current norm of the displacement residual
    
    TDataType mLMRatioTolerance;      /// The ratio threshold for the norm of the LM  residual
    TDataType mLMAbsTolerance;        /// The absolute value threshold for the norm of the LM  residual
    TDataType mLMInitialResidualNorm; /// The reference norm of the LM residual
    TDataType mLMCurrentResidualNorm; /// The current norm of the LM residual
    
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

#endif /* KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_CONTACT_CRITERIA_H */

