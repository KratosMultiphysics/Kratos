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

#if !defined(KRATOS_FANCY_RESIDUAL_CRITERIA )
#define  KRATOS_FANCY_RESIDUAL_CRITERIA

// System includes 

// External includes 

// Project includes 
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
class FancyResidualCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions 
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FancyResidualCriteria );

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

    //* Constructor.
    
    FancyResidualCriteria(
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
        mResRatioTolerance = NewRatioTolerance;
        mResAbsTolerance = AlwaysConvergedNorm;
        mInitialResidualIsSet = false;
    }

    //* Copy constructor.
    
    FancyResidualCriteria( FancyResidualCriteria const& rOther )
      :BaseType(rOther) 
      ,mInitialResidualIsSet(rOther.mInitialResidualIsSet)
      ,mResRatioTolerance(rOther.mResRatioTolerance)
      ,mInitialResidualNorm(rOther.mInitialResidualNorm)
      ,mCurrentResidualNorm(rOther.mCurrentResidualNorm)
      ,mResAbsTolerance(rOther.mResAbsTolerance)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
      ,mpTable(rOther.mpTable)
      ,mStandaloneTable(rOther.mStandaloneTable)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
    {
    }

    //* Destructor.
    
    ~FancyResidualCriteria() override {}

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
        
        if (TSparseSpace::Size(b) != 0) //if we are solving for something
        {
            if (mInitialResidualIsSet == false)
            {
                mInitialResidualNorm = TSparseSpace::TwoNorm(b);
                mInitialResidualIsSet = true;
            }

            TDataType res_ratio;
            mCurrentResidualNorm = TSparseSpace::TwoNorm(b);

            if(mInitialResidualNorm == 0.0)
            {
                res_ratio = 0.0;
            }
            else
            {
                res_ratio = mCurrentResidualNorm/mInitialResidualNorm;
            }

            const double b_size = TSparseSpace::Size(b);
            TDataType res_abs = (mCurrentResidualNorm/b_size);
            
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {                
                if (mpTable != nullptr)
                {
                    std::cout.precision(4);
                    auto& Table = mpTable->GetTable();
                    Table  << res_ratio  << mResRatioTolerance  << res_abs  << mResAbsTolerance;
                }
                else
                {
                    std::cout.precision(4);
                #if !defined(_WIN32)
                    std::cout << BOLDFONT("RESIDUAL CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                    std::cout << BOLDFONT("\tRESIDUAL: RATIO = ") << res_ratio << BOLDFONT(" EXP.RATIO = ") << mResRatioTolerance << BOLDFONT(" ABS = ") << res_abs << BOLDFONT(" EXP.ABS = ") << mResAbsTolerance << std::endl;
                #else
                    std::cout << "RESIDUAL CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                    std::cout << "\tRESIDUAL: RATIO = " << res_ratio << " EXP.RATIO = " << mResRatioTolerance << " ABS = " << res_abs << " EXP.ABS = " << mResAbsTolerance << std::endl;
                #endif
                }
            }

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = res_ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = res_abs;

            if (res_ratio <= mResRatioTolerance || res_abs < mResAbsTolerance)
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
                        std::cout << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    #else
                        std::cout << "\tResidual convergence is achieved" << std::endl;
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
        else
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
            Table.AddColumn("RES. RAT", 10);
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
        
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0 && mStandaloneTable == true)
        {
            std::cout.precision(4);
        #if !defined(_WIN32)
            std::cout << "\n\n" << BOLDFONT("CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
        #else
            std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
        #endif
                
            if (mpTable != nullptr)
            {
                mpTable->PrintHeader();
            }
        }
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

    bool mInitialResidualIsSet;      // The initial residual is set

    TDataType mResRatioTolerance;    // The tolerance ratio

    TDataType mInitialResidualNorm;  // The initial residual norm

    TDataType mCurrentResidualNorm;  // The current residual norm

    TDataType mResAbsTolerance;      // The absolute converged norm

    TDataType mReferenceDispNorm;    // The reference displacement norm

    TablePrinterPointerType mpTable; // Pointer to the fancy table 
    
    bool mStandaloneTable;           // If the table is not appended to any other table
    
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

}; // Class ClassName 

///@} 

///@name Type Definitions 
///@{


///@} 

}  // namespace Kratos.

#endif // KRATOS_NEW_DISPLACEMENT_CRITERIA  defined 

