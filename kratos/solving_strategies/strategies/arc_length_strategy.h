//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Ignasi Pouplana
//

#if !defined(KRATOS_ARC_LENGTH_STRATEGY)
#define KRATOS_ARC_LENGTH_STRATEGY

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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
 * @class ArcLengthStrategy
 * @ingroup KratosCore
 * @brief This is the base ArcLengthStrategy
 * @details The theoretical details can be found in "Geometrical interpretation of the ARC-LENGTH method", M. Fafard and B. Massicotte,
 * Computers and Structures Vol 46 pp 603-615 (1993). (Ramm arc-length method)
 * @author Alejandro Cornejo
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ArcLengthStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
    public:
    ///@name Type Definitions
    ///@{

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(ArcLengthStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
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

    // Constructor
    ArcLengthStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            mDesiredIterations = rParameters["desired_iterations"].GetInt();
            mMaxRadiusFactor   = rParameters["max_radius_factor"].GetDouble();
            mMinRadiusFactor   = rParameters["min_radius_factor"].GetDouble();
            mInitializeArcLengthWasPerformed = false;

            // we initialize the list of load processes to be taken into account
            if (rParameters["loads_sub_model_part_list"].size() > 0) {
                mSubModelPartList.resize(rParameters["loads_sub_model_part_list"].size());
                mVariableNames.resize(rParameters["loads_variable_list"].size());

                if (mSubModelPartList.size() != mVariableNames.size())
                    KRATOS_THROW_ERROR( std::logic_error, "For each SubModelPart there must be a corresponding nodal Variable", "")

                for (unsigned int i = 0; i < mVariableNames.size(); i++) {
                    mSubModelPartList[i] = &( model_part.GetSubModelPart(rParameters["loads_sub_model_part_list"][i].GetString()));
                    mVariableNames[i] = rParameters["loads_variable_list"][i].GetString();
                }
            }
        }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~ArcLengthStrategy() override
    {
        BaseType::~BaseType();
    }


    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (!BaseType::mInitializeWasPerformed)
        {
            BaseType::Initialize();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        SparseSpaceType::Clear(mpf);
        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);

        TSystemVectorType& mf      = *mpf;
        TSystemVectorType& mDxf    = *mpDxf;
        TSystemVectorType& mDxb    = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;

        SparseSpaceType::Resize(mf, 0);
        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        BaseType::Clear();

        KRATOS_CATCH("");
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (!mInitializeArcLengthWasPerformed) {
            //set up the system
            if (!mpBuilderAndSolver->GetDofSetIsInitializedFlag()) {
                //setting up the list of the DOFs to be solved
                mpBuilderAndSolver->SetUpDofSet(mpScheme, BaseType::GetModelPart());

                //shaping correctly the system
                mpBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
            }

            // Compute initial radius (mRadius_0)
            mpBuilderAndSolver->ResizeAndInitializeVectors(mpScheme, mpA, mpDx, mpb, BaseType::GetModelPart());
            TSystemMatrixType& mA = *mpA;
            TSystemVectorType& mDx = *mpDx;
            TSystemVectorType& mb = *mpb;
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

            mRadius_0 = TSparseSpace::TwoNorm(mDx);
            mRadius = mRadius_0;

            // Compute vector of reference external force (mf)
            this->InitializeSystemVector(mpf);
            TSystemVectorType& mf = *mpf;
            TSparseSpace::SetToZero(mf);

            // We build it now to only include external loads
            mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mf);

            //Initialize the loading factor Lambda
            mLambda = 0.0;
            mLambda_old = 1.0;

            // Initialize Norm of solution
            mNormxEquilibrium = 0.0;

            mInitializeArcLengthWasPerformed = true;

            KRATOS_INFO("Arc Length Strategy") << "Strategy Initialized" << std::endl;
        }

        BaseType::InitializeSolutionStep();

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;


        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_INFO("Arc-Length Strategy") << "INITIAL ARC-LENGTH RADIUS: " << mRadius_0 << std::endl;
        KRATOS_INFO("Arc-Length Strategy") << "ARC-LENGTH RADIUS: " << mRadius/mRadius_0 << " X initial radius" << std::endl;
        KRATOS_INFO("Arc-Length Strategy") << "ARC-LENGTH LAMBDA: " << mLambda << std::endl;

        ModelPart& r_model_part = BaseType::GetModelPart();

        // Initialize variables
		DofsArrayType& r_dof_set    = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& r_A      = *mpA;
        TSystemVectorType& r_Dx     = *mpDx;
        TSystemVectorType& r_b      = *mpb;
        TSystemVectorType& r_f      = *mpf;
        TSystemVectorType& r_Dxb    = *mpDxb;
        TSystemVectorType& r_Dxf    = *mpDxf;
        TSystemVectorType& r_DxPred = *mpDxPred;
        TSystemVectorType& r_DxStep = *mpDxStep;

        // Initialize iterations info
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), r_A, r_Dx, r_b);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, r_A, r_Dx, r_b);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, r_A, r_Dx, r_b);






        return true;
    }

    void UpdateExternalLoads()
    {
        // Update External Loads
        for (unsigned int i = 0; i < mVariableNames.size(); i++) {
            ModelPart& r_sub_model_part = *(mSubModelPartList[i]);
            const std::string& r_variable_name = mVariableNames[i];

            if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
                const Variable<double>& var = KratosComponents< Variable<double> >::Get( r_variable_name );
                block_for_each(model_part.Nodes(), [&](auto& r_node){
                    r_node.FastGetSolutionStepValue(var) *= (mLambda/mLambda_old);
                });
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(r_variable_name)) {
                typedef Variable<array_1d<double,3>> array_type;
                const array_type& r_var = KratosComponents<array_type>::Get(r_variable_name);

                block_for_each(model_part.Conditions(), [&](auto& r_condition){
                    r_condition.GetValue(r_var) *= (mLambda/mLambda_old);
                });

                // TODO-> add for node loads

            } else {
                KRATOS_THROW_ERROR( std::logic_error, "One variable of the applied loads has a non supported type. Variable: ", r_variable_name )
            }
        }

        // Save the applied Lambda factor
        mLambda_old = mLambda;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                : "arc_length_strategy",
            "use_old_stiffness_in_first_iteration": false,
            "max_iteration"                       : 10,
            "reform_dofs_at_each_step"            : false,
            "compute_reactions"                   : false,
            "desired_iterations"                  : 4,
            "max_radius_factor"                   : 10.0,
            "min_radius_factor"                   : 0.1,
            "builder_and_solver_settings"         : {},
            "convergence_criteria_settings"       : {},
            "linear_solver_settings"              : {},
            "scheme_settings"                     : {}
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
        return "arc_length_strategy";
    }

    /**
     * @brief It resizes and initializes a system vector
     */
    void InitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (pv == NULL)
        {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), false);
    }

    /**
     * @brief It saves system vector pointer
     */
    void SaveInitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (pv == NULL)
        {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), true);
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
        return "ArcLengthStrategy";
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

    TSystemVectorPointerType mpf;      /// Vector of reference external forces
    TSystemVectorPointerType mpDxf;    /// Delta x of A*Dxf=f
    TSystemVectorPointerType mpDxb;    /// Delta x of A*Dxb=b
    TSystemVectorPointerType mpDxPred; /// Delta x of prediction phase
    TSystemVectorPointerType mpDxStep; /// Delta x of the current step

    unsigned int mDesiredIterations;   /// This is used to calculate the radius of the next step

    bool mInitializeArcLengthWasPerformed;

    double mMaxRadiusFactor, mMinRadiusFactor; /// Used to limit the radius of the arc length strategy
    double mRadius, mRadius_0;                 /// Radius of the arc length strategy
    double mLambda, mLambda_old;               /// current and old loading factor
    double mNormxEquilibrium;                  /// Norm of the solution vector in equilibrium
    double mDLambdaStep;                       /// Delta lambda of the current step

    std::vector<ModelPart*> mSubModelPartList; /// List of  SubModelParts associated to an external load
    std::vector<std::string> mVariableNames;   /// Name of the nodal variable associated to each SubModelPart


    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief Here the database is updated
     * @param A The LHS matrix of the system of equations
     * @param Dx The incremement in the solution
     * @param b The RHS vector of the system of equations
     * @param MoveMesh The flag that allows to move the mesh
     */

    virtual void UpdateDatabase(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh)
    {
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        p_scheme->Update(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        // Move the mesh if needed
        if (MoveMesh == true)
            BaseType::MoveMesh();
        

        // TODO
    }


    /**
     * @brief This method prints information after reach the max number of iterations
     */

    virtual void MaxIterationsExceeded()
    {
        KRATOS_INFO_IF("ARC-LENGTH Strategy", this->GetEchoLevel() > 0)
            << "ATTENTION: max iterations ( " << mMaxIterationNumber
            << " ) exceeded!" << std::endl;
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

    ArcLengthStrategy(const ArcLengthStrategy &Other){};

    ///@}

}; /* Class ArcLengthStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_ARC_LENGTH_STRATEGY  defined */