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

#if !defined(KRATOS_MORTAR_AND_CRITERIA_H)
#define  KRATOS_MORTAR_AND_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/color_utilities.h"
#include "custom_utilities/bprinter_utility.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"

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
class MortarAndConvergenceCriteria : public And_Criteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /** Counted pointer of MortarAndConvergenceCriteria */

    KRATOS_CLASS_POINTER_DEFINITION(MortarAndConvergenceCriteria );

    typedef And_Criteria< TSparseSpace, TDenseSpace >         BaseType;

    typedef TSparseSpace                               SparseSpaceType;

    typedef typename BaseType::TDataType                     TDataType;

    typedef typename BaseType::DofsArrayType             DofsArrayType;

    typedef typename BaseType::TSystemMatrixType     TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType     TSystemVectorType;

    typedef boost::shared_ptr<BprinterUtility> TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** 
     * Constructor.
     */
    MortarAndConvergenceCriteria(
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pFirstCriterion,
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pSecondCriterion,
        TablePrinterPointerType pTable = nullptr
        )
        :And_Criteria< TSparseSpace, TDenseSpace >(pFirstCriterion, pSecondCriterion),
        mpTable(pTable)
    {
    }

    /**
     * Copy constructor.
     */
    MortarAndConvergenceCriteria(MortarAndConvergenceCriteria const& rOther)
      :BaseType(rOther)
     {
         BaseType::mpFirstCriterion   =  rOther.mpFirstCriterion;
         BaseType::mpSecondCriterion  =  rOther.mpSecondCriterion;      
     }

    /** Destructor.
    */
    virtual ~MortarAndConvergenceCriteria () {}

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
                const unsigned int NLIteration = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
                mpTable->AddToRow<unsigned int>(NLIteration);
            }
        }
        
        bool CriterionResult = BaseType::PostCriteria(rModelPart,rDofSet,A,Dx,b);
        
        if (CriterionResult == true && rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            if (mpTable != nullptr)
            {
                mpTable->PrintFooter();
            }
        }
        
        return CriterionResult;
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        if (mpTable != nullptr)
        {
            auto& Table = mpTable->GetTable();
            Table.AddColumn("ITER", 4);
//             mpTable->AddColumn("ITER", 4);
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
            std::cout << "\n\n" << BOLDFONT("CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl << std::scientific;
                
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
    
    TablePrinterPointerType mpTable;

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

