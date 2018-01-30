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

#if !defined(KRATOS_MORTAR_AND_CRITERIA_H)
#define  KRATOS_MORTAR_AND_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/table_stream_utility.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "utilities/color_utilities.h"
#include "utilities/condition_number_utility.h"

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

/** Short class definition.
Detail class definition.

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}

*/

template<class TSparseSpace,
         class TDenseSpace
         >
class MortarAndConvergenceCriteria 
    : public And_Criteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /** Counted pointer of MortarAndConvergenceCriteria */

    KRATOS_CLASS_POINTER_DEFINITION(MortarAndConvergenceCriteria );

    typedef And_Criteria< TSparseSpace, TDenseSpace >                BaseType;

    typedef TSparseSpace                                      SparseSpaceType;

    typedef typename TSparseSpace::MatrixType                SparseMatrixType;

    typedef typename TSparseSpace::VectorType                SparseVectorType;

    typedef typename TDenseSpace::MatrixType                  DenseMatrixType;

    typedef typename TDenseSpace::VectorType                  DenseVectorType;
    
    typedef typename BaseType::TDataType                            TDataType;

    typedef typename BaseType::DofsArrayType                    DofsArrayType;

    typedef typename BaseType::TSystemMatrixType            TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType            TSystemVectorType;
    
    typedef TableStreamUtility::Pointer               TablePrinterPointerType;
    
    typedef ConditionNumberUtility::Pointer ConditionNumberUtilityPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** 
     * Constructor.
     */
    MortarAndConvergenceCriteria(
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pFirstCriterion,
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pSecondCriterion,
        TablePrinterPointerType pTable = nullptr,
        const bool PrintingOutput = false,
        ConditionNumberUtilityPointerType pConditionNumberUtility = nullptr
        )
        :And_Criteria< TSparseSpace, TDenseSpace >(pFirstCriterion, pSecondCriterion),
        mpTable(pTable),
        mPrintingOutput(PrintingOutput),
        mpConditionNumberUtility(pConditionNumberUtility),
        mTableIsInitialized(false)
    {
    }

    /**
     * Copy constructor.
     */
    MortarAndConvergenceCriteria(MortarAndConvergenceCriteria const& rOther)
      :BaseType(rOther)
      ,mpTable(rOther.mpTable)
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
      ,mpConditionNumberUtility(rOther.mpConditionNumberUtility)
     {
         BaseType::mpFirstCriterion  = rOther.mpFirstCriterion;
         BaseType::mpSecondCriterion = rOther.mpSecondCriterion;      
     }

    /** Destructor.
    */
    ~MortarAndConvergenceCriteria () override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Criteria that need to be called after getting the solution
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
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            if (mpTable != nullptr)
                mpTable->AddToRow<unsigned int>(rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER]);
        
        bool criterion_result = BaseType::PostCriteria(rModelPart, rDofSet, A, Dx, b);
        
        if (mpConditionNumberUtility != nullptr)
        {
            TSystemMatrixType copy_A; // NOTE: Can not be const, TODO: Change the solvers to const
            const double condition_number = mpConditionNumberUtility->GetConditionNumber(copy_A);
            
            if (mpTable != nullptr)
            {
                std::cout.precision(4);
                auto& Table = mpTable->GetTable();
                Table  << condition_number;
            }
            else
            {
                if (mPrintingOutput == false)
                    std::cout << "\n" << BOLDFONT("CONDITION NUMBER:") << "\t " << std::scientific << condition_number << std::endl;
                else
                    std::cout << "\n" << "CONDITION NUMBER:" << "\t" << std::scientific << condition_number << std::endl;
            }
        }
        
        if (criterion_result == true && rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            if (mpTable != nullptr)
                mpTable->PrintFooter();
        
        return criterion_result;
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        if (mpTable != nullptr && mTableIsInitialized == false)
        {
            (mpTable->GetTable()).SetBold(!mPrintingOutput);
            (mpTable->GetTable()).AddColumn("ITER", 4);
        }
        
        mTableIsInitialized = true;
        BaseType::Initialize(rModelPart);
         
        if (mpTable != nullptr && mpConditionNumberUtility != nullptr)
            (mpTable->GetTable()).AddColumn("COND.NUM.", 10);
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
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            std::cout.precision(4);
            if (mPrintingOutput == false)
                std::cout << "\n\n" << BOLDFONT("CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[STEP] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            else
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[STEP] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
                
            if (mpTable != nullptr)
                mpTable->PrintHeader();
        }
        
        BaseType::InitializeSolutionStep(rModelPart,rDofSet,A,Dx,b);
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
        BaseType::FinalizeSolutionStep(rModelPart,rDofSet,A,Dx,b);
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
    
    TablePrinterPointerType mpTable;                            // Pointer to the fancy table 
    bool mPrintingOutput;                                       // If the colors and bold are printed
    ConditionNumberUtilityPointerType mpConditionNumberUtility; // The utility to compute the condition number
    bool mTableIsInitialized;                                   // If the table is already initialized
    
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

}; /* Class ClassName */

///@}

///@name Type Definitions */
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MORTAR_AND_CRITERIA_H  defined */

