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
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit LineSearchStrategy(ModelPart& rModelPart, Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters["move_mesh_flag"].GetBool())
    {
    }

    /**
     * Constructor.
     */
    explicit LineSearchStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ): ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag)
    {}

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
    ): ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag)
    {}

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
    ///@name Friends
    ///@{


    ///@}

private:
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

protected:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

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

        TSystemVectorType aux(b.size()); //TODO: do it by using the space
        TSparseSpace::Assign(aux,0.5, Dx);

        //compute residual without update
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
        double ro = TSparseSpace::TwoNorm(b);

        //compute half step residual
        BaseType::UpdateDatabase(A,aux,b,MoveMesh);
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
        double rh = TSparseSpace::TwoNorm(b);

        //compute full step residual (add another half Dx to the previous half)
        BaseType::UpdateDatabase(A,aux,b,MoveMesh);
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
        double rf = TSparseSpace::TwoNorm(b);

        //compute optimal (limited to the range 0-1)
        //parabola is y = a*x^2 + b*x + c -> min/max for
        //x=0   --> r=ro
        //x=1/2 --> r=rh
        //x=1   --> r =
        //c= ro,     b= 4*rh -rf -3*ro,  a= 2*rf - 4*rh + 2*ro
        //max found if a>0 at the position  xmax = (rf/4 - rh)/(rf - 2*rh);
        double parabola_a = 2*rf + 2*ro - 4*rh;
        double parabola_b = 4*rh - rf - 3*ro;
        double xmin = 1e-3;
        double xmax = 1.0;
        if( parabola_a > 0) //if parabola has a local minima
        {
            xmax = -0.5 * parabola_b/parabola_a; // -b / 2a
            if( xmax > 1.0)
                xmax = 1.0;
            else if(xmax < -1.0)
                xmax = -1.0;
        }
        else //parabola degenerates to either a line or to have a local max. best solution on either extreme
        {
            if(rf < ro)
                xmax = 1.0;
            else
                xmax = xmin; //should be zero, but otherwise it will stagnate
        }

        //perform final update
        TSparseSpace::Assign(aux,-(1.0-xmax), Dx);
        BaseType::UpdateDatabase(A,aux,b,MoveMesh);
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
