//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

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

/**
 * @class LineSearchStrategy
 * @ingroup KratosCore
 * @brief This strategy extends the standard Newton-Raphson method by incorporating a line search algorithm to improve convergence.
 * @details The Line Search strategy is a variant of the Newton-Raphson method that aims to enhance convergence for highly non-linear problems.
 * In each non-linear iteration, after computing the search direction (the solution increment) $\Delta\mathbf{u}$, this strategy finds an optimal scaling factor $\alpha$ to update the solution. The update is performed as:
 * $$
 * \mathbf{u}_{k+1} = \mathbf{u}_k + \alpha \cdot \Delta\mathbf{u}_k
 * $$
 * The goal of the line search is to determine an $\alpha$ that minimizes the energy potential (or equivalently, the component of the residual) along the search direction $\Delta\mathbf{u}_k$. This is achieved by an iterative procedure that approximates the minimum.
 * This approach can prevent divergence in cases where the full Newton step ($\alpha=1$) would be too large, significantly improving the robustness and stability of the non-linear solution process. The behavior of the line search is controlled by several parameters, such as the maximum number of line search iterations, initial guesses for alpha, and the desired tolerance.
 * @see ResidualBasedNewtonRaphsonStrategy
 * @author Riccardo Rossi
 * @tparam TSparseSpace The sparse space definition
 * @tparam TDenseSpace The dense space definition
 * @tparam TLinearSolver The linear solver definition
 */
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

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef LineSearchStrategy<TSparseSpace,TDenseSpace,TLinearSolver> ClassType;

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
     * @brief Default constructor
     */
    explicit LineSearchStrategy() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit LineSearchStrategy(ModelPart& rModelPart, Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

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
    explicit LineSearchStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : BaseType(rModelPart, pScheme, pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag)
    {
        Parameters default_settings = this->GetDefaultParameters();
        OverrideDefaultSettingsWithParameters(default_settings, MaxIterations, ReformDofSetAtEachStep, CalculateReactions);
        this->AssignSettings(default_settings);
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit LineSearchStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        Parameters default_settings = this->GetDefaultParameters();
        OverrideDefaultSettingsWithParameters(default_settings, MaxIterations, ReformDofSetAtEachStep, CalculateReactions);
        this->AssignSettings(default_settings);
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    KRATOS_DEPRECATED_MESSAGE("Constructor deprecated, please use the constructor without linear solver")
    explicit LineSearchStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : BaseType(rModelPart, pScheme, pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag)
    {
        Parameters default_settings = this->GetDefaultParameters();
        OverrideDefaultSettingsWithParameters(default_settings, MaxIterations, ReformDofSetAtEachStep, CalculateReactions);
        this->AssignSettings(default_settings);
    }

    /**
     * Constructor with Settings
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param ThisParameters Settings used in the strategy
     */
    LineSearchStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters ThisParameters
        ): BaseType(rModelPart, pScheme, pNewLinearSolver,pNewConvergenceCriteria, BaseType::GetDefaultParameters())
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * Constructor with Settings and pointer to BuilderAndSolver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param ThisParameters Settings used in the strategy
     */
    LineSearchStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters ThisParameters
        ): BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, BaseType::GetDefaultParameters())
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
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

    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                       : "line_search_strategy",
            "max_line_search_iterations" : 5,
            "first_alpha_value"          : 0.5,
            "second_alpha_value"         : 1.0,
            "min_alpha"                  : 0.1,
            "max_alpha"                  : 2.0,
            "line_search_tolerance"      : 0.5
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "line_search_strategy";
    }

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

    int mMaxLineSearchIterations;   //Maximum number of iterations line search do
    double mFirstAlphaValue;        //First alpha guess value used for the first iteration
    double mSecondAlphaValue;       //Second alpha guess value used for the first iteration
    double mMinAlpha;               //Minimum possible alpha value at the end of the algorithm
    double mMaxAlpha;               //Maximum possible alpha value at the end of the algorithm
    double mLineSearchTolerance;    //Tolerance of the line search algorithm, defined as the ratio between
                                    //maximum residual*alpha*dx and current iteration residual*alpha*dx
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
                KRATOS_INFO("LineSearchStrategy") << "LINE SEARCH it " << it << " coeff = " << x <<  " r1 = " << r1 << " r2 = " << r2 << std::endl;
                converged = true;
                TSparseSpace::Assign(aux,x, Dx);
                break;
            }

            //note that we compute the next residual only if it is strictly needed (we break on the line before if it is not needed)
            TSparseSpace::SetToZero(b);
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
            double rf = TSparseSpace::Dot(aux,b);

            KRATOS_INFO("LineSearchStrategy") << "LINE SEARCH it " << it << " coeff = " << x << " rf = " << rf << " r1 = " << r1 << " r2 = " << r2 << std::endl;


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
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
        mMaxLineSearchIterations = ThisParameters["max_line_search_iterations"].GetInt();
        mFirstAlphaValue = ThisParameters["first_alpha_value"].GetDouble();
        mSecondAlphaValue = ThisParameters["second_alpha_value"].GetDouble();
        mMinAlpha = ThisParameters["min_alpha"].GetDouble();
        mMaxAlpha = ThisParameters["max_alpha"].GetDouble();
        mLineSearchTolerance = ThisParameters["line_search_tolerance"].GetDouble();
    }

    /**
     * @brief This method overrides the default settings with user provided parameters
     * @param rDefaultSettings Parameters with default settings
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     */
    void OverrideDefaultSettingsWithParameters(
        Parameters& rDefaultSettings,
        const double MaxIterations,
        const bool ReformDofSetAtEachStep,
        const bool CalculateReactions
        )
    {
        rDefaultSettings["max_iteration"].SetInt(MaxIterations);
        rDefaultSettings["reform_dofs_at_each_step"].SetBool(ReformDofSetAtEachStep);
        rDefaultSettings["compute_reactions"].SetBool(CalculateReactions);
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
