//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
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
 * @author Alejandro Cornejo and Ignasi de Pouplana
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
    typedef typename BaseType::TConvergenceCriteriaType TConvergenceCriteriaType;
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
        Parameters ThisParameters
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewConvergenceCriteria, pNewBuilderAndSolver, ThisParameters["max_iteration"].GetInt(),
                ThisParameters["compute_reactions"].GetBool(), ThisParameters["reform_dofs_at_each_step"].GetBool(),
                ThisParameters["move_mesh_flag"].GetBool())
        {
            ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
            AssignSettings(ThisParameters);
        }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~ArcLengthStrategy() override = default;


    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        BaseType::AssignSettings(ThisParameters);
        mDesiredIterations = ThisParameters["desired_iterations"].GetInt();
        mMaxRadiusFactor   = ThisParameters["max_radius_factor"].GetDouble();
        mMinRadiusFactor   = ThisParameters["min_radius_factor"].GetDouble();
        mInitializeArcLengthWasPerformed = false;

        // we initialize the list of load processes to be taken into account
        if (ThisParameters["loads_sub_model_part_list"].size() > 0) {
            mSubModelPartList.resize(ThisParameters["loads_sub_model_part_list"].size());
            mVariableNames.resize(ThisParameters["loads_variable_list"].size());
        
            KRATOS_ERROR_IF(mSubModelPartList.size() != mVariableNames.size()) << "For each SubModelPart there must be a corresponding nodal Variable" << std::endl;
        
            const auto& r_sub_model_parts_names = ThisParameters["loads_sub_model_part_list"].GetStringArray();
            for (std::size_t i = 0; i < mVariableNames.size(); i++) {
                mSubModelPartList[i] = &( r_model_part.GetSubModelPart(r_sub_model_parts_names[i]));
                mVariableNames[i] = ThisParameters["loads_variable_list"][i].GetString();
            }
        }
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;
        BaseType::Initialize();
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
            ModelPart& r_model_part = BaseType::GetModelPart();
            //set up the system
            if (!this->mpBuilderAndSolver->GetDofSetIsInitializedFlag()) {
                // Setting up the list of the DOFs to be solved
                this->mpBuilderAndSolver->SetUpDofSet(this->mpScheme, r_model_part);

                // Shaping correctly the system
                this->mpBuilderAndSolver->SetUpSystem(r_model_part);
            }

            // Compute initial radius (mRadius_0)
            this->mpBuilderAndSolver->ResizeAndInitializeVectors(this->mpScheme, this->mpA, this->mpDx, this->mpb, r_model_part);
            TSystemMatrixType& rA  = *(this->mpA);
            TSystemVectorType& rDx = *(this->mpDx);
            TSystemVectorType& rb  = *(this->mpb);
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            this->mpBuilderAndSolver->BuildAndSolve(this->mpScheme, r_model_part, rA, rDx, rb);

            mRadius_0 = TSparseSpace::TwoNorm(rDx);
            mRadius = mRadius_0;

            // Compute vector of reference external force (mf)
            this->InitializeSystemVector(mpf);
            TSystemVectorType& rf = *mpf;
            TSparseSpace::SetToZero(rf);

            // We build it now to only include external loads
            this->mpBuilderAndSolver->BuildRHS(this->mpScheme, r_model_part, rf);

            //Initialize the loading factor Lambda
            mLambda = 0.0;
            mLambda_old = 1.0;

            // Initialize Norm of solution
            mNormxEquilibrium = 0.0;

            mInitializeArcLengthWasPerformed = true;

            KRATOS_INFO_IF("ArcLengthStrategy", BaseType::GetEchoLevel() > 0) << "Strategy Initialized" << std::endl;
        }
    
        BaseType::InitializeSolutionStep();
        SaveInitializeSystemVector(mpf);
        InitializeSystemVector(mpDxf);
        InitializeSystemVector(mpDxb);
        InitializeSystemVector(mpDxPred);
        InitializeSystemVector(mpDxStep);

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& r_model_part  = BaseType::GetModelPart();

        const std::size_t iteration_number = r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER];
    
        // Update the radius
        mRadius = mRadius * std::sqrt(double(mDesiredIterations) / double(iteration_number));
        if (mRadius > mMaxRadiusFactor*mRadius_0) {
            mRadius = mMaxRadiusFactor*mRadius_0;
        } else if (mRadius < mMinRadiusFactor*mRadius_0) {
            mRadius = mMinRadiusFactor*mRadius_0;
        }

        BaseType::FinalizeSolutionStep();

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_INFO_IF("ArcLengthStrategy", BaseType::GetEchoLevel() > 0) << "INITIAL ARC-LENGTH RADIUS: " << mRadius_0 << std::endl;
        KRATOS_INFO_IF("ArcLengthStrategy", BaseType::GetEchoLevel() > 0) << "ARC-LENGTH RADIUS: " << mRadius/mRadius_0 << " X initial radius" << std::endl;
        mInsideIterationLoop = false;

        ModelPart& r_model_part = BaseType::GetModelPart();

        // Initialize variables
        DofsArrayType& r_dof_set    = this->mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& r_A      = *(this->mpA);
        TSystemVectorType& r_Dx     = *(this->mpDx);
        TSystemVectorType& r_b      = *(this->mpb);
        TSystemVectorType& r_f      = *mpf;
        TSystemVectorType& r_Dxb    = *mpDxb;
        TSystemVectorType& r_Dxf    = *mpDxf;
        TSystemVectorType& r_DxPred = *mpDxPred;
        TSystemVectorType& r_DxStep = *mpDxStep;

        // Initialize iterations info
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

        this->mpScheme->InitializeNonLinIteration(r_model_part, r_A, r_Dx, r_b);
        this->mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, r_A, r_Dx, r_b);
        bool is_converged = this->mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, r_A, r_Dx, r_b);

        TSparseSpace::SetToZero(r_A);
        TSparseSpace::SetToZero(r_b);
        TSparseSpace::SetToZero(r_Dxf);

        // Now we compute Dxf
        this->mpBuilderAndSolver->Build(this->mpScheme, r_model_part, r_A, r_b);
        this->mpBuilderAndSolver->ApplyDirichletConditions(this->mpScheme, r_model_part, r_A, r_Dx, r_b);
        TSparseSpace::Assign(r_b, 1.0, r_f);
        this->mpBuilderAndSolver->SystemSolve(r_A, r_Dxf, r_b);
        double lambda_increment = mRadius / TSparseSpace::TwoNorm(r_Dxf);
        mDLambdaStep = lambda_increment;
        mLambda += lambda_increment;

        KRATOS_INFO_IF("ArcLengthStrategy", BaseType::GetEchoLevel() > 0) << "ARC-LENGTH LAMBDA: " << mLambda << std::endl;

        TSparseSpace::InplaceMult(r_Dxf, lambda_increment);
        TSparseSpace::Assign(r_DxPred, 1.0, r_Dxf);
        TSparseSpace::Assign(r_DxStep, 1.0, r_DxPred);
        TSparseSpace::InplaceMult(r_Dxf, 1.0 / lambda_increment);
        UpdateDatabase(r_A, r_DxPred, r_b, BaseType::MoveMeshFlag());

        this->mpScheme->FinalizeNonLinIteration(r_model_part, r_A, r_DxPred, r_b);
        this->mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, r_A, r_Dx, r_b);

        if (is_converged) {
            if (this->mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(r_b);
                this->mpBuilderAndSolver->BuildRHS(this->mpScheme, r_model_part, r_b);
            }

            is_converged = this->mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, r_A, r_Dxf, r_b);
        }

        while (!is_converged && iteration_number++ < BaseType::mMaxIterationNumber) {
            mInsideIterationLoop = true;
            //setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            this->mpScheme->InitializeNonLinIteration(r_model_part, r_A, r_Dx, r_b);
            this->mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, r_A, r_Dx, r_b);
            is_converged = this->mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, r_A, r_Dx, r_b);

            TSparseSpace::SetToZero(r_A);
            TSparseSpace::SetToZero(r_b);
            TSparseSpace::SetToZero(r_Dxf);

            // We compute r_Dxf 
            this->mpBuilderAndSolver->Build(this->mpScheme, r_model_part, r_A, r_b);
            this->mpBuilderAndSolver->ApplyDirichletConditions(this->mpScheme, r_model_part, r_A, r_Dxf, r_b);
            TSparseSpace::Assign(r_b, 1.0, r_f);
            this->mpBuilderAndSolver->SystemSolve(r_A, r_Dxf, r_b);

            TSparseSpace::SetToZero(r_A);
            TSparseSpace::SetToZero(r_b);
            TSparseSpace::SetToZero(r_Dxb);

            this->mpBuilderAndSolver->BuildAndSolve(this->mpScheme, r_model_part, r_A, r_Dxb, r_b);
            lambda_increment = -TSparseSpace::Dot(r_DxPred, r_Dxb) / TSparseSpace::Dot(r_DxPred, r_Dxf);
            TSparseSpace::Assign(r_Dx, 1.0, r_Dxb + lambda_increment*r_Dxf);

            // Update results
            mDLambdaStep += lambda_increment;
            mLambda += lambda_increment;
            TSparseSpace::UnaliasedAdd(r_DxStep, 1.0, r_Dx);
            UpdateDatabase(r_A, r_Dx, r_b, BaseType::MoveMeshFlag());

            this->mpScheme->FinalizeNonLinIteration(r_model_part, r_A, r_Dx, r_b);
            this->mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, r_A, r_Dx, r_b);

            if (is_converged) {
                if (this->mpConvergenceCriteria->GetActualizeRHSflag()) {
                    TSparseSpace::SetToZero(r_b);
                    this->mpBuilderAndSolver->BuildRHS(this->mpScheme, r_model_part, r_b);
                }
                is_converged = this->mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, r_A, r_Dx, r_b);
            }
        }
        // Prints a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber) {
            MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("Arc-Length Strategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << BaseType::mMaxIterationNumber << " iterations" << std::endl;
        }
        //calculate reactions if required
        if (BaseType::mCalculateReactionsFlag)
            this->mpBuilderAndSolver->CalculateReactions(this->mpScheme, r_model_part, r_A, r_Dx, r_b);
        return is_converged;
    }

    /**
     * @brief This method updates the value of the external load according to the new load factor Lambda
     */
    void UpdateExternalLoads()
    {
        // Update External Loads
        for (unsigned int i = 0; i < mVariableNames.size(); i++) {
            ModelPart& r_sub_model_part = *(mSubModelPartList[i]);
            const std::string& r_variable_name = mVariableNames[i];

            if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
                const Variable<double>& var = KratosComponents<Variable<double>>::Get(r_variable_name);
                block_for_each(r_sub_model_part.Nodes(), [&](Node& r_node){
                    r_node.FastGetSolutionStepValue(var) *= (mLambda/mLambda_old);
                });
            } else if (KratosComponents<Variable<array_1d<double,3>>>::Has(r_variable_name)) {
                typedef Variable<array_1d<double,3>> array_type;
                const array_type& r_var = KratosComponents<array_type>::Get(r_variable_name);

                block_for_each(r_sub_model_part.Conditions(), [&](Condition& r_condition) {
                    /* we are multipying by Lambda instead of Lambda/Lambda_old because
                    the load processes reset the loads to its initial values
                    when the InitSolStep is called */
                    if (mInsideIterationLoop)
                        r_condition.GetValue(r_var) *= mLambda / mLambda_old;
                    else
                        r_condition.GetValue(r_var) *= mLambda;
                });

                // TODO-> add for node loads

            } else {
                KRATOS_ERROR << "One variable of the applied loads has a non supported type. Variable: " << r_variable_name << std::endl;
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
            "name"                                 : "arc_length_strategy",
            "desired_iterations"                   : 4,
            "max_radius_factor"                    : 10.0,
            "min_radius_factor"                    : 0.1,
            "loads_sub_model_part_list"            : [],
            "loads_variable_list"                  : []
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
        if (!pv) {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != this->mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(this->mpBuilderAndSolver->GetEquationSystemSize(), false);
    }

    /**
     * @brief It saves system vector pointer
     */
    void SaveInitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (!pv) {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }
        TSystemVectorType& v = *pv;
        if (v.size() != this->mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(this->mpBuilderAndSolver->GetEquationSystemSize(), true);
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
    bool mInsideIterationLoop;

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
    void UpdateDatabase(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh) override
    {
        BaseType::UpdateDatabase(rA, rDx, rb, MoveMesh);
        UpdateExternalLoads();
    }


    /**
     * @brief This method prints information after reach the max number of iterations
     */

    void MaxIterationsExceeded() override
    {
        KRATOS_INFO_IF("ARC-LENGTH Strategy", this->GetEchoLevel() > 0) << "ATTENTION: max iterations ("<< this->mMaxIterationNumber <<") exceeded!" << std::endl;
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
