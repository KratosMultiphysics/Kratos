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
//

#if !defined(KRATOS_FANCY_DISPLACEMENT_CRITERIA )
#define  KRATOS_FANCY_DISPLACEMENT_CRITERIA

/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "utilities/bprinter_utility.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
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
class FancyDisplacementCriteria : virtual public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FancyDisplacementCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >  BaseType;

    typedef TSparseSpace                               SparseSpaceType;

    typedef typename BaseType::TDataType                     TDataType;

    typedef typename BaseType::DofsArrayType             DofsArrayType;

    typedef typename BaseType::TSystemMatrixType     TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType     TSystemVectorType;
    
    typedef std::size_t                                        KeyType;
    
    typedef boost::shared_ptr<BprinterUtility> TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
    */
    FancyDisplacementCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm,
        TablePrinterPointerType pTable = nullptr,
        const bool StandaloneTable = true 
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mpTable(pTable),
          mStandaloneTable(StandaloneTable),
          mTableIsInitialized(false)
    {
        mDispRatioTolerance = NewRatioTolerance;
        mDispAbsTolerance = AlwaysConvergedNorm;
    }

    /** Copy constructor.
    */
    FancyDisplacementCriteria( FancyDisplacementCriteria const& rOther )
      :BaseType(rOther)
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
      ,mpTable(rOther.mpTable)
      ,mStandaloneTable(rOther.mStandaloneTable)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
    {
    }

    /** Destructor.
    */
    ~FancyDisplacementCriteria() override {}


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
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0 && mStandaloneTable == true)
        {
            if (mpTable != nullptr)
            {
                const unsigned int nl_iteration = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
                mpTable->AddToRow<unsigned int>(nl_iteration);
            }
        }
        
        bool criterion_result;
        
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            TDataType final_correction_norm = TSparseSpace::TwoNorm(Dx);

            TDataType disp_ratio = 0.0;

            CalculateReferenceNorm(rDofSet);

            if(final_correction_norm == 0)
            {
                disp_ratio = 0.0;
            }
            else
            {
                if(mReferenceDispNorm == 0)
                {
                    KRATOS_ERROR << "NaN norm is detected" << std::endl;
                }
                
                disp_ratio = final_correction_norm/mReferenceDispNorm;
            }
            
            const double size_dx = SparseSpaceType::Size(Dx);

            const double disp_abs = (final_correction_norm/std::sqrt(size_dx));

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                if (mpTable != nullptr)
                {
                    std::cout.precision(4);
                    auto& Table = mpTable->GetTable();
                    Table  << disp_ratio  << mDispRatioTolerance  << disp_abs  << mDispAbsTolerance;
                }
                else
                {
                    std::cout.precision(4);
                #if !defined(_WIN32)
                    std::cout << BOLDFONT("DISPLACEMENT CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                    std::cout << BOLDFONT("\tDISPLACEMENT: RATIO = ") << disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                #else
                    std::cout << "DISPLACEMENT CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                    std::cout << "\tDISPLACEMENT: RATIO = " << disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                #endif
                }
            }

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = disp_ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = disp_abs;

            if ( disp_ratio <= mDispRatioTolerance  ||  disp_abs < mDispAbsTolerance )  //  || (final_correction_norm/x.size())<=1e-7)
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
                        std::cout << BOLDFONT("\tDisplacement") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    #else
                        std::cout << "\tDisplacement convergence is achieved" << std::endl;
                    #endif
                    }
                }

                criterion_result = true;
            }
            else
            {
                criterion_result = false;
            }
        }
        else //in this case all the displacements are imposed!
        {
            criterion_result = true;
        }
        
        if (criterion_result == true && rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0 && mStandaloneTable == true)
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
     * @param rModelPart Reference to the ModelPart containing the contact problem. (unused)
     */
        
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
        
        if (mpTable != nullptr && mTableIsInitialized == false)
        {
            auto& Table = mpTable->GetTable();
            if (mStandaloneTable == true) Table.AddColumn("ITER", 4);
            Table.AddColumn("DP RATIO", 10);
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

    TDataType mDispRatioTolerance;   // The tolerance ratio

    TDataType mDispAbsTolerance;     // The absolute converged norm

    TDataType mReferenceDispNorm;    // The reference displacement norm
    
    TablePrinterPointerType mpTable; // Pointer to the fancy table 
    
    bool mStandaloneTable;           // If the table is not appended to any other table
    
    bool mTableIsInitialized;        // If the table is already initialized

    ///@}
    ///@name Private Operators
    ///@{

    void CalculateReferenceNorm(DofsArrayType& rDofSet)
    {
        mReferenceDispNorm = TDataType();
        TDataType temp;

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                temp = i_dof->GetSolutionStepValue();
                mReferenceDispNorm += temp*temp;
            }
        }
        mReferenceDispNorm = std::sqrt(mReferenceDispNorm);
    }

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

///@name Type Definitions
///@{


///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_FANCY_DISPLACEMENT_CRITERIA  defined */

