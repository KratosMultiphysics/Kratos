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


#if !defined(KRATOS_UPWIND_LINE_SEARCH_STRATEGY )
#define  KRATOS_UPWIND_LINE_SEARCH_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "upwind_residualbased_newton_raphson_strategy.hpp"


namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class UpwindLineSearchStrategy
    : public UpwindResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(UpwindLineSearchStrategy);

    typedef UpwindResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

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
     * Default Constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    UpwindLineSearchStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ): UpwindResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag)
    {
        Parameters default_settings = this->GetDefaultSettings();
        OverrideDefaultSettingsWithParameters(default_settings, MaxIterations, ReformDofSetAtEachStep, CalculateReactions);
        this->AssignSettings(default_settings);
    }

    /**
     * Destructor.
     */

    ~UpwindLineSearchStrategy() override
    {
    }
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "UpwindLineSearchStrategy";
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


private:
    ///@name Protected member Variables
    ///@{
    int mMaxLineSearchIterations;   //Maximum number of iterations line search do
    double mFirstAlphaValue;        //First alpha guess value used for the first iteration
    double mSecondAlphaValue;       //Second alpha guess value used for the first iteration
    double mMinAlpha;               //Minimum possible alpha value at the end of the algorithm
    double mMaxAlpha;               //Maximum possible alpha value at the end of the algorithm
    double mLineSearchTolerance;    //Tolerance of the line search algorithm, defined as the ratio between
                                    //maximum residual*alpha*dx and current iteration residual*alpha*dx
    ///@}

protected:
    ///@name Private Operators
    ///@{

    /**
     * Here the database is updated
     * @param A The LHS matrix of the system of equations
     * @param Dx The incremement in the solution
     * @param b The RHS vector of the system of equations
     * @param MoveMesh The flag that allows to move the mesh
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

        TSystemVectorType aux(Dx);

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
                KRATOS_INFO("UpwindLineSearchStrategy") << "LINE SEARCH it " << it << " coeff = " << x <<  " r1 = " << r1 << " r2 = " << r2 << std::endl;
                converged = true;
                TSparseSpace::Assign(aux,x, Dx);
                break;
            }

            //note that we compute the next residual only if it is strictly needed (we break on the line before if it is not needed)
            TSparseSpace::SetToZero(b);
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
            double rf = TSparseSpace::Dot(aux,b);

            KRATOS_INFO("UpwindLineSearchStrategy") << "LINE SEARCH it " << it << " coeff = " << x << " rf = " << rf << " r1 = " << r1 << " r2 = " << r2 << std::endl;


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

    /**
     * @brief This method returns the default settings
     */
    Parameters GetDefaultSettings() override
    {
        Parameters base_default_settings = BaseType::GetDefaultSettings();
        Parameters default_settings(R"({
            "max_line_search_iterations" : 5,
            "first_alpha_value"          : 0.5,
            "second_alpha_value"         : 1.0,
            "min_alpha"                  : 0.1,
            "max_alpha"                  : 2.0,
            "line_search_tolerance"      : 0.5
        })");
        default_settings.AddMissingParameters(base_default_settings);
        return default_settings;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param Settings Parameters that are assigned to the member variables
     */
    void AssignSettings(Parameters Settings) override
    {
        BaseType::AssignSettings(Settings);
        mMaxLineSearchIterations = Settings["max_line_search_iterations"].GetInt();
        mFirstAlphaValue = Settings["first_alpha_value"].GetDouble();
        mSecondAlphaValue = Settings["second_alpha_value"].GetDouble();
        mMinAlpha = Settings["min_alpha"].GetDouble();
        mMaxAlpha = Settings["max_alpha"].GetDouble();
        mLineSearchTolerance = Settings["line_search_tolerance"].GetDouble();
    }

    /**
     * @brief This method overrides the default settings wiht user provided parameters
     * @param DefaultSettings Parameters with default settings
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     */
    void OverrideDefaultSettingsWithParameters(
        Parameters DefaultSettings,
        const double MaxIterations,
        const bool ReformDofSetAtEachStep,
        const bool CalculateReactions
    )
    {
        DefaultSettings["max_iterations"].SetInt(MaxIterations);
        DefaultSettings["reform_dofs_at_each_step"].SetBool(ReformDofSetAtEachStep);
        DefaultSettings["calculate_reactions"].SetBool(CalculateReactions);
    }



    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * Copy constructor.
     */

    UpwindLineSearchStrategy(const UpwindLineSearchStrategy& Other)
    {
    };


    ///@}

}; /* Class UpwindLineSearchStrategy */

///@}

} /* namespace Kratos. */

#endif /* KRATOS_UPWIND_LINE_SEARCH_STRATEGY  defined */
