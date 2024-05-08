//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolas Sibuet
//


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
class AnnPromLineSearchStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(AnnPromLineSearchStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef AnnPromLineSearchStrategy<TSparseSpace,TDenseSpace,TLinearSolver> ClassType;

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
    explicit AnnPromLineSearchStrategy() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit AnnPromLineSearchStrategy(ModelPart& rModelPart, Parameters ThisParameters)
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
    explicit AnnPromLineSearchStrategy(
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
    explicit AnnPromLineSearchStrategy(
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
    explicit AnnPromLineSearchStrategy(
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
    AnnPromLineSearchStrategy(
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
    AnnPromLineSearchStrategy(
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

    ~AnnPromLineSearchStrategy() override
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
            "name"                       : "ann_prom_wolfe_line_search_strategy",
            "max_line_search_expand_iterations" : 10,
            "max_line_search_zoom_iterations" : 100,
            "coefficient_1_value"          : 1e-4,
            "second_alpha_value"         : 0.9,
            "max_alpha"                  : 10.0
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
        return "ann_prom_line_search_strategy";
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
        return "AnnPromLineSearchStrategy";
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
        TSystemVectorType& Dq_orig = BaseType::GetModelPart().GetRootModelPart().GetValue(ROM_SOLUTION_INCREMENT);
        TSystemVectorType& q_base = BaseType::GetModelPart().GetRootModelPart().GetValue(ROM_SOLUTION_BASE);
        TSystemVectorType& q_total = BaseType::GetModelPart().GetRootModelPart().GetValue(ROM_SOLUTION_TOTAL);

        TSystemVectorType& x_old = BaseType::GetModelPart().GetRootModelPart().GetValue(SOLUTION_BASE);
        x_old -= Dx;

        // KRATOS_INFO("AnnPromLineSearchStrategy") << "ROM_SOLUTION_INCREMENT " << DxRom << std::endl;
        // KRATOS_INFO("AnnPromLineSearchStrategy") << "ROM_SOLUTION_BASE " << XRomPrevious << std::endl;
        // KRATOS_INFO("AnnPromLineSearchStrategy") << "ROM_SOLUTION_TOTAL " << XRomTotal << std::endl;
        // KRATOS_INFO("AnnPromLineSearchStrategy") << "SOLUTION_BASE " << UPrevious[8174] << std::endl;

        TSystemVectorType p = Dq_orig;
        // TSystemVectorType p(Dq_orig);
        // TSparseSpace::Assign(p, 1.0/TSparseSpace::TwoNorm(Dq_orig) ,Dq_orig);

        TDataType alpha_old = 0.0;
        
        TDataType metric_base;
        TSystemVectorType grad_base(q_base);
        ComputeMetricAndGradient(q_base, p, x_old, alpha_old, globalPhiEffective, A, b, MoveMesh, metric_base, grad_base);
        TDataType proj_base = TSparseSpace::Dot(grad_base, p);

        TDataType metric_old = metric_base;
        TDataType proj_old = proj_base;

        TDataType alpha = mInitialAlpha;

        TDataType metric;
        TSystemVectorType grad(grad_base);
        TDataType proj;

        bool converged = false;
        int it = 0;
        while(!converged && it < mMaxExpandIterations) {

            ComputeMetricAndGradient(q_base, p, x_old, alpha, globalPhiEffective, A, b, MoveMesh, metric, grad);
            proj = TSparseSpace::Dot(grad,p)
            
            // If new point violates sufficient decrease criteria or is higher than last one, it denotes a valid search interval
            if((metric > metric_base + mCoef1*alpha*proj_base) || (it>0 && metric > metric_old)) {
                KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH found valid interval. x_low: " << alpha_old << "(f: " << metric_old << "), x_high: " << alpha << "(f: " << metric << ")" << std::endl;
                alpha = Zoom(metric_old, proj_old, alpha_old, metric, proj, alpha, q_base, metric_base, grad_base, p, x_old, globalPhiEffective, A, b, MoveMesh);
                converged = true;
                continue;
            }

            // If new point satisfies strong Wolfe conditions, we use it as our result
            if(std::abs(proj) <= mCoef2*std::abs(proj_base)) {
                KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH found valid point. Returning x: " << alpha << std::endl;
                converged = true;
                continue;
            }

            // If the new point has positive slope, then it denotes a valid search interval
            if(proj >= 0.0) {
                KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH found valid interval. x_low: " << alpha << "(f: " << metric << "), x_high: " << alpha_old << "(f: " << metric_old << ")" << std::endl;
                alpha = Zoom(metric, proj, alpha, metric_old, proj_old, alpha_old, q_base, metric_base, grad_base, p, x_old, globalPhiEffective, A, b, MoveMesh);
                converged = true;
                continue;
            }

            metric_old = metric;
            proj_old = proj;
            alpha_old = alpha;

            alpha = 2.0*alpha;
            if (alpha >= mAlphaMax){
                KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH went over the maximum x value. Returning x: " << mApphaMax << std::endl;
                alpha = mAlphaMax;
                converged = true;
                continue;
            }

            it++;
        }

        ComputeMetricAndGradient(q_base, p, x_old, alpha, globalPhiEffective, A, b, MoveMesh, metric, grad);
        KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH final value x: " << alpha << "(f: " << metric << ")" << std::endl;

        TSparseSpace::Assign(Dq_orig, alpha, p);
        q_total += Dq_orig;

        TSparseSpace::SetToZero(A);
        TSparseSpace::SetToZero(b);

    }

   /**
    * @brief This method searches for a point that satisfies the Strong Wolfe Conditions, within the specified interval.
    * @param metric_low Initial low metric value
    * @param proj_low Initial projection value at the low point
    * @param alpha_low Initial alpha value at the low point
    * @param metric_high Initial high metric value
    * @param proj_high Initial projection value at the high point
    * @param alpha_high Initial alpha value at the high point
    * @param q_base The original reduced snapshot vector
    * @param metric_base Metric function evaluated at q_base
    * @param grad_base Gradient of the metric function evaluated at q_base (in reduced space)
    * @param p The search direction vector (in reduced space)
    * @param x_old The values in the current system snapshot
    * @param globalPhiEffective System's global phi effective
    * @param A The LHS matrix of the system of equations
    * @param b The RHS vector of the system of equations
    * @param MoveMesh The flag that allows to move the mesh
   */
    TDataType Zoom(
        TDataType metric_low,
        TDataType proj_low,
        TDataType alpha_low,
        TDataType metric_high,
        TDataType proj_high,
        TDataType alpha_high,
        const TSystemVectorType& q_base,
        const TDataType& metric_base,
        const TSystemVectorType& grad_base,
        const TSystemVectorType& p,
        const TSystemVectorType& x_old,
        const TSystemMatrixType& globalPhiEffective,
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const bool MoveMesh
    )
    {
        TDataType proj_base = TSparseSpace::Dot(grad_base,p);

        TDataType metric;
        TSystemVectorType grad(grad_base);
        TDataType proj;
        
        bool converged = false;
        int it = 0;
        while(!converged && it < mMaxZoomIterations) {
            
            alpha = 0.5*(alpha_low+alpha_high);

            ComputeMetricAndGradient(q_base, p, x_old, alpha, globalPhiEffective, A, b, MoveMesh, metric, grad);
            proj =  TSparseSpace::Dot(grad,p);

            // If new point violates sufficient decrease condition or is higher than metric_low, set it as new high point
            if((metric > metric_base + mCoef1*alpha*proj_base) || (metric >= metric_low)) {
                KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH reduced interval to: x_low: " << alpha_low << "(f: " << metric_low << "), x_high: " << alpha << "(f: " << metric << ")" << std::endl;
                alpha_high = alpha;
                metric_high = metric;
                proj_high = proj;

            } else {

                // If the new point satisfies the curvature condition, then return it as our result
                if(std::abs(proj) <= mCoef2*std::abs(proj_base)) {
                    KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH found valid point. Returning x: " << alpha << std::endl;
                    return alpha;
                }
                
                // If proj has same sign as alpha_high - alpha_low, consider the low point as the new high point
                if(proj*(alpha_high-alpha_low) >= 0) {
                    alpha_high = alpha_low;
                    metric_high = metric_low;
                    proj_high = proj_low;
                }

                alpha_low = alpha;
                metric_low = metric;
                proj_low = proj;
            }

            it++;
        }

        KRATOS_INFO("AnnPromLineSearchStrategy") << "LINE SEARCH passed the iterations limit. Returning x: " << alpha << std::endl;
        return alpha;
    }
                
    /**
     * @brief This method assigns the value of the metric at the specific point and its gradient. IMPORTANT-> It updates the database and x_old
     * @
     * @param q_base The original reduced snapshot vector
     * @param p The search direction vector (in reduced space)
     * @param x_old The values in the current system snapshot
     * @param alpha Stepsize to apply
     * @param globalPhiEffective System's global phi effective
     * @param A The LHS matrix of the system of equations
     * @param b The RHS vector of the system of equations
     * @param MoveMesh The flag that allows to move the mesh
     * @param metric Output: Metric function evaluated at q_base
     * @param grad Output: Gradient of the metric function evaluated at q_base (in reduced space)
     * 
     */
    void ComputeMetricAndGradient(
        const TSystemVectorType& q_base,
        const TSystemVectorType& p,
        const TSystemVectorType& x_old,
        const TDataType& alpha,
        const TSystemMatrixType& globalPhiEffective,
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const bool MoveMesh,
        TDataType& metric,
        TSystemVectorType& grad
    )
    {
        typename TSchemeType::Pointer pScheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = this->GetBuilderAndSolver();
        TSparseSpace::SetToZero(A);
        TSparseSpace::SetToZero(b);

        TSystemVectorType Dq(p);
        TSystemVectorType x_total(x_old);
        TSystemMatrixType phi(globalPhiEffective);

        TSparseSpace::Assign(Dq, alpha, p);
        pBuilderAndSolver->GetXAndDecoderGradient(q_base+Dq, x_total, phi);
        const TSystemVectorType Dx = x_total-x_old;
        BaseType::UpdateDatabase(A, Dx, b, MoveMesh);
        x_old = x_total;
        pBuilderAndSolver->Build(pScheme, BaseType::GetModelPart(), A, b);

        TSystemVectorType aux(x_old);
        TSparseSpace::TransposeMult(A, b, aux);
        TSparseSpace::TransposeMult(phi, aux, grad);
        TSparseSpace::Assign(grad, -1.0, grad);
        metric = TSparseSpace::Dot(b, b);
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

    AnnPromLineSearchStrategy(const AnnPromLineSearchStrategy& Other)
    {
    };


    ///@}

}; /* Class AnnPromLineSearchStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos. */
