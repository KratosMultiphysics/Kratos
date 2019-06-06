//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_LINE_SEARCH_STRATEGY )
#define  KRATOS_LINE_SEARCH_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"


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

/// Short class definition.

/// Detail class definition.

//URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

//URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

//URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

//URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


//URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

//URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

//URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

//URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}



template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class LineSearchStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(LineSearchStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Constructor.
     */

    LineSearchStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ): ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag),
       mSettings(Parameters(R"({})"))
    {
        CheckDefaultsAndAssignSettings();
    }

    LineSearchStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ): ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag),
       mSettings(Parameters(R"({})"))
    {
        CheckDefaultsAndAssignSettings();
    }

    LineSearchStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters Settings,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ): ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag),
       mSettings(Settings)
    {
        CheckDefaultsAndAssignSettings();
    }

    LineSearchStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters Settings,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ): ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag),
       mSettings(Settings)
    {
        CheckDefaultsAndAssignSettings();
    }

    /**
     * Destructor.
     */

    ~LineSearchStrategy() override
    {
    }



    ///@}
    ///@name Operators

    ///@{

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
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "LineSearchStrategy";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    Parameters mSettings;
    int mMaxLineSearchIterations;
    double mFirstAlphaValue;
    double mSecondAlphaValue;
    double mMinAlpha;
    double mMaxAlpha;
    double mLineSearchTolerance;


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

protected:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{
    void CheckDefaultsAndAssignSettings()
    {
        Parameters default_settings(R"({
            "max_line_search_iterations" : 5,
            "first_alpha_value"          : 0.5,
            "second_alpha_value"         : 1.0,
            "min_alpha"                  : 0.1,
            "max_alpha"                  : 2.0,
            "line_search_tolerance"      : 0.5
        })");
        mSettings.ValidateAndAssignDefaults(default_settings);
        mMaxLineSearchIterations = mSettings["max_line_search_iterations"].GetInt();
        mFirstAlphaValue = mSettings["first_alpha_value"].GetDouble();
        mSecondAlphaValue = mSettings["second_alpha_value"].GetDouble();
        mMinAlpha = mSettings["min_alpha"].GetDouble();
        mMaxAlpha = mSettings["max_alpha"].GetDouble();
        mLineSearchTolerance = mSettings["line_search_tolerance"].GetDouble();
        
    }

    /**
     * Here the database is updated
     */

    void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh
    ) override
    {
        typename TSchemeType::Pointer pScheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = this->GetBuilderAndSolver();

        TSystemVectorType aux(TSparseSpace::Size(b));
        
        double x1 = mFirstAlphaValue;
        double x2 = mSecondAlphaValue;

        bool converged = false;
        int it = 0;
        double xprevious = 0.0;

        //Compute residual with 1 coefficient update (x1)
        //since no update was performed yet, this includes an increment wrt the previous
        //solution of x1*Dx
        TSparseSpace::Assign(aux,x1-xprevious, Dx);
        xprevious = x1;
        BaseType::UpdateDatabase(A,aux,b,MoveMesh); 
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
        double r1 = TSparseSpace::Dot(aux,b);
        
        double rmax = std::abs(r1);
        while(!converged && it < mMaxLineSearchIterations) {

            //Compute residual with 2 coefficient update (x2)
            //since the database was initialized with x1*Dx
            //we need to apply ONLY THE INCREMENT, that is (x2-xprevious)*Dx
            TSparseSpace::Assign(aux,x2-xprevious, Dx);
            xprevious = x2;
            BaseType::UpdateDatabase(A,aux,b,MoveMesh);
            TSparseSpace::SetToZero(b);
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
            double r2 = TSparseSpace::Dot(aux,b);

            if(it == 0) {
                rmax = std::max(rmax,std::abs(r2));
            }
            double rmin = std::min(std::abs(r1),std::abs(r2));

            //Find optimum
            double x = 1.0;
            if(std::abs(r1 - r2) > 1e-10)
                x =  (r1*x2 - r2*x1)/(r1 - r2);
            
            if(x < mMinAlpha) {
                x = mMinAlpha;
            } else if(x > mMaxAlpha) {
                x = mMaxAlpha;
            }                

            //Perform final update
            TSparseSpace::Assign(aux,x-xprevious, Dx);
            xprevious = x;
            BaseType::UpdateDatabase(A,aux,b,MoveMesh);
            if(rmin < mLineSearchTolerance*rmax) {
                std::cout << "LINE SEARCH it " << it << " coeff = " << x <<  " r1 = " << r1 << " r2 = " << r2 << std::endl;
                converged = true;
                TSparseSpace::Assign(aux,x, Dx);
                break;
            }

            //note that we compute the next residual only if it is strictly needed (we break on the line before if it is not needed)
            TSparseSpace::SetToZero(b);
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
            double rf = TSparseSpace::Dot(aux,b);

            std::cout << "LINE SEARCH it " << it << " coeff = " << x << " rf = " << rf << " r1 = " << r1 << " r2 = " << r2 << std::endl;


            if(std::abs(rf) < rmax*mLineSearchTolerance) {
                converged = true;
                TSparseSpace::Assign(aux,x, Dx);
            } else {
                if(std::abs(r1)>std::abs(r2)) {
                    r1 = rf;
                    x1 = x;
                } else {
                    r2 = r1;
                    x2 = x1;
                    r1 = rf;
                    x1 = x;
                }
                converged = false;
            }


            it++;
        }
        TSparseSpace::SetToZero(b);
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

    /**
     * Copy constructor.
     */

    LineSearchStrategy(const LineSearchStrategy& Other)
    {
    };


    ///@}

}; /* Class LineSearchStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos. */

#endif /* KRATOS_LINE_SEARCH_STRATEGY  defined */
