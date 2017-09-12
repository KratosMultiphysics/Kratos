//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//  

#if !defined(KRATOS_FANCY_AND_CRITERIA_H)
#define  KRATOS_FANCY_AND_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/table_stream_utility.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#if !defined(_WIN32)
    #include "utilities/color_utilities.h"
#endif

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
class FancyAndConvergenceCriteria : public And_Criteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /** Counted pointer of FancyAndConvergenceCriteria */

    KRATOS_CLASS_POINTER_DEFINITION(FancyAndConvergenceCriteria );

    typedef And_Criteria< TSparseSpace, TDenseSpace >           BaseType;

    typedef TSparseSpace                                 SparseSpaceType;

    typedef typename BaseType::TDataType                       TDataType;

    typedef typename BaseType::DofsArrayType               DofsArrayType;

    typedef typename BaseType::TSystemMatrixType       TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType       TSystemVectorType;

    typedef boost::shared_ptr<TableStreamUtility> TableStreamPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** 
     * Constructor.
     */
    FancyAndConvergenceCriteria(
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pFirstCriterion,
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pSecondCriterion,
        TableStreamPointerType pTable = nullptr, 
        const bool PrintingOutput = false 
        )
        :And_Criteria< TSparseSpace, TDenseSpace >(pFirstCriterion, pSecondCriterion),
        mpTable(pTable),
        mPrintingOutput(PrintingOutput),
        mTableIsInitialized(false)
    {
    }

    /**
     * Copy constructor.
     */
    FancyAndConvergenceCriteria(FancyAndConvergenceCriteria const& rOther)
      :BaseType(rOther)
     {
         BaseType::mpFirstCriterion  = rOther.mpFirstCriterion;
         BaseType::mpSecondCriterion = rOther.mpSecondCriterion;
         mpTable = rOther.mpTable;
         mPrintingOutput = rOther.mPrintingOutput;
         mTableIsInitialized = rOther.mTableIsInitialized;
     }

    /** Destructor.
    */
    ~FancyAndConvergenceCriteria () override = default;

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
        {
            if (mpTable != nullptr)
            {
                const unsigned int nl_iteration = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
                mpTable->AddToRow<unsigned int>(nl_iteration);
            }
        }
        
        const bool criterion_result = BaseType::PostCriteria(rModelPart,rDofSet,A,Dx,b);
        
        if (criterion_result == true && rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            if (mpTable != nullptr)
            {
                mpTable->PrintFooter();
            }
        }
        
        return criterion_result;
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        if (mpTable != nullptr && mTableIsInitialized == false)
        {
            auto& Table = mpTable->GetTable();
            Table.AddColumn("ITER", 4);
            mTableIsInitialized = true;
        }
        
         BaseType::Initialize(rModelPart);
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
            {
            #if !defined(_WIN32)
                std::cout << "\n\n" << BOLDFONT("CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            #else
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            #endif
            }
            else
            {
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            }
                
            if (mpTable != nullptr)
            {
                mpTable->PrintHeader();
            }
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
    
    TableStreamPointerType mpTable; // Pointer to the fancy table 
    
    bool mPrintingOutput;            // If the colors and bold are printed
    
    bool mTableIsInitialized;        // If the table is already initialized
    
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

}; /* Class FancyAndConvergenceCriteria */

///@}

///@name Type Definitions */
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_FANCY_AND_CRITERIA_H  defined */

