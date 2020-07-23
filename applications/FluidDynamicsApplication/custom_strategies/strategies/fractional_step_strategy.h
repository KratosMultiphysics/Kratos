//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_FRACTIONAL_STEP_STRATEGY
#define KRATOS_FRACTIONAL_STEP_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/solver_settings.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

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
 * @brief Fractional-step strategy for incompressible Navier-Stokes formulation
 * This strategy implements a splitting scheme for the incompressible Navier-Stokes equations.
 * It is intended to be used in combination with the FractionalStep element in the FluidDynamicsApplication.
 * The fractional step index, which is stored in the ProcessInfo, takes the values
 * 1 : Momentum step (calculate fractional step velocity)
 * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
 * 4 : Pressure step
 * 5 : Computation of projections
 * 6 : End of step velocity
 * @tparam TSparseSpace Sparse space template type
 * @tparam TDenseSpace Dense space template type
 * @tparam TLinearSolver Linear solver template type
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class FractionalStepStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FractionalStepStrategy
    KRATOS_CLASS_POINTER_DEFINITION(FractionalStepStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

    typedef SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver> SolverSettingsType;

    ///@}
    ///@name Life Cycle
    ///@{

    FractionalStepStrategy(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector):
        BaseType(rModelPart,false),
        mCalculateReactionsFlag(false),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        KRATOS_WARNING("FractionalStepStrategy") << "This constructor is deprecated. Use the one with the \'CalculateReactionsFlag\' instead." << std::endl;
        InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FractionalStepStrategy(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector,
               const Kratos::Variable<int>& PeriodicVar):
        BaseType(rModelPart,false),
        mCalculateReactionsFlag(false),
        mrPeriodicIdVar(PeriodicVar)
    {
        KRATOS_WARNING("FractionalStepStrategy") << "This constructor is deprecated. Use the one with the \'CalculateReactionsFlag\' instead." << std::endl;
        InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FractionalStepStrategy(
        ModelPart& rModelPart,
        SolverSettingsType& rSolverConfig,
        bool PredictorCorrector,
        bool CalculateReactionsFlag)
        : BaseType(rModelPart,false)
        , mCalculateReactionsFlag(CalculateReactionsFlag)
        , mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FractionalStepStrategy(
        ModelPart& rModelPart,
        SolverSettingsType& rSolverConfig,
        bool PredictorCorrector,
        bool CalculateReactionsFlag,
        const Kratos::Variable<int>& PeriodicVar)
        : BaseType(rModelPart,false)
        , mCalculateReactionsFlag(CalculateReactionsFlag)
        , mrPeriodicIdVar(PeriodicVar)
    {
        InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    /// Destructor.
    ~FractionalStepStrategy() override{}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        // Set up nodes to use slip conditions if needed.
        if (mUseSlipConditions) {
            auto& r_model_part = BaseType::GetModelPart();
            const int n_conds = r_model_part.NumberOfConditions();
#pragma omp parallel for
            for (int i_cond = 0; i_cond < n_conds; ++i_cond) {
                auto it_cond = r_model_part.ConditionsBegin() + i_cond;
                if (it_cond->Is(SLIP)) {
                    auto& r_geom = it_cond->GetGeometry();
                    for (auto& r_node : r_geom) {
                        r_node.SetLock();
                        r_node.Set(SLIP, true);
                        r_node.UnSetLock();
                    }
                }
            }
        }
    }

    int Check() override
    {
        KRATOS_TRY;

        // Base strategy check
        int ierr = BaseType::Check();
        if (ierr != 0) {
            return ierr;
        }

        // Check time order and buffer size
        const auto& r_model_part = BaseType::GetModelPart();
        KRATOS_ERROR_IF(mTimeOrder == 2 && r_model_part.GetBufferSize() < 3)
            << "Buffer size too small for fractional step strategy (BDF2), needed 3, got " << r_model_part.GetBufferSize() << std::endl;
        KRATOS_ERROR_IF(mTimeOrder == 1 && r_model_part.GetBufferSize() < 2)
            << "Buffer size too small for fractional step strategy (Backward Euler), needed 2, got " << r_model_part.GetBufferSize() << std::endl;

        // Check elements and conditions
        const auto &r_current_process_info = r_model_part.GetProcessInfo();
        for (const auto& r_element : r_model_part.Elements()) {
            ierr = r_element.Check(r_current_process_info);
            if (ierr != 0) {
               break;
            }
        }

        for (const auto& r_condition : r_model_part.Conditions()) {
            ierr = r_condition.Check(r_current_process_info);
            if (ierr != 0) {
                break;
            }
        }

        return ierr;

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override
    {
        // Initialize BDF2 coefficients
        SetTimeCoefficients();
    }

    bool SolveSolutionStep() override
    {
        bool converged = false;
        if (mPredictorCorrector) {
            const unsigned int echo_level = BaseType::GetEchoLevel();
            // Iterative solution for pressure
            for (unsigned int it = 0; it < mMaxPressureIter; ++it) {
                KRATOS_INFO_IF("FractionalStepStrategy", echo_level > 1) << "Pressure iteration " << it << std::endl;
                const auto convergence_output = this->SolveStep();
                converged = this->CheckPressureConvergence(std::get<1>(convergence_output));
                if (converged) {
                    KRATOS_INFO_IF("FractionalStepStrategy", echo_level > 0) << "Predictor-corrector converged in " << it + 1 << " iterations." << std::endl;
                    break;
                }
            }
            KRATOS_WARNING_IF("FractionalStepStrategy", !converged && echo_level > 0) << "Predictor-corrector iterations did not converge." << std::endl;
        } else {
            // Solve for fractional step velocity, then update pressure once
            const auto convergence_output = this->SolveStep();
            // If not doing predictor corrector iterations, norm_dp will
            // typically be "large" since we are not iterating on pressure.
            // It makes no sense to report that the iteration didn't converge
            // based on this. Hence, what we report is the convergence of the
            // fractional step velocity.
            converged = std::get<0>(convergence_output);
        }

        // Calculate reactions
        if (mCalculateReactionsFlag) {
            CalculateReactions();
        }

        return converged;
    }

    void FinalizeSolutionStep() override
    {
        if (mReformDofSet) {
            this->Clear();
        }
    }

    //TODO: Move to private section as soon as we remove the Python exposure
    /**
     * @brief Calculates the reactions
     * This methods calculates the reactions of the momentum equation.
     * These are computed as minus the RHS and saved in the REACTION variable
     */
    virtual void CalculateReactions()
    {
        auto &r_model_part = BaseType::GetModelPart();
        auto &r_process_info = r_model_part.GetProcessInfo();
        const int n_elems = r_model_part.NumberOfElements();

        // Set fractional step index to the momentum equation step
        const int original_step = r_process_info[FRACTIONAL_STEP];
        r_process_info.SetValue(FRACTIONAL_STEP, 1);

        // Allocate and initialize values for REACTION calculation
        LocalSystemVectorType RHS_Contribution;
        LocalSystemMatrixType LHS_Contribution;
        const auto &r_const_process_info = r_process_info;
        VariableUtils().SetHistoricalVariableToZero(REACTION, r_model_part.Nodes());

#pragma omp parallel for private(RHS_Contribution, LHS_Contribution)
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            // Build local system
            auto it_elem = r_model_part.ElementsBegin() + i_elem;
            it_elem->CalculateLocalSystem(
                LHS_Contribution,
                RHS_Contribution,
                r_const_process_info);

            // Accumulate minus the RHS as the reaction
            unsigned int index = 0;
            auto& r_geom = it_elem->GetGeometry();
            const unsigned int n_nodes = r_geom.PointsNumber();
            for (unsigned int i = 0; i < n_nodes; ++i) {
                r_geom[i].SetLock();
                auto& r_reaction = r_geom[i].FastGetSolutionStepValue(REACTION);
                for (unsigned int d = 0; d < mDomainSize; ++d) {
                    r_reaction[d] -= RHS_Contribution[index++];
                }
                r_geom[i].UnSetLock();
            }
        }

        // Synchronize the local REACTION values
        r_model_part.GetCommunicator().AssembleCurrentData(REACTION);

        // Reset original fractional step index
        r_process_info.SetValue(FRACTIONAL_STEP, original_step);
    }

    virtual void AddIterationStep(Process::Pointer pNewStep)
    {
        mExtraIterationSteps.push_back(pNewStep);
    }

    virtual void ClearExtraIterationSteps()
    {
        mExtraIterationSteps.clear();
    }

    void Clear() override
    {
        mpMomentumStrategy->Clear();
        mpPressureStrategy->Clear();
    }

    ///@}
    ///@name Access
    ///@{

    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        int StrategyLevel = Level > 0 ? Level - 1 : 0;
        mpMomentumStrategy->SetEchoLevel(StrategyLevel);
        mpPressureStrategy->SetEchoLevel(StrategyLevel);
    }

    /**
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "FractionalStepStrategy" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected Life Cycle
    ///@{


    ///@}
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    double mVelocityTolerance;

    double mPressureTolerance;

    unsigned int mMaxVelocityIter;

    unsigned int mMaxPressureIter;

    unsigned int mDomainSize;

    unsigned int mTimeOrder;

    bool mPredictorCorrector;

    bool mUseSlipConditions;

    bool mReformDofSet;

    bool mCalculateReactionsFlag;

    /// Scheme for the solution of the momentum equation
    StrategyPointerType mpMomentumStrategy;

    /// Scheme for the solution of the mass equation
    StrategyPointerType mpPressureStrategy;

    std::vector< Process::Pointer > mExtraIterationSteps;

    const Kratos::Variable<int>& mrPeriodicIdVar;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Set the Time Coefficients object
     * Calculate the coefficients for the BDF2 time iteration.
     * These are stored in the BDF_COEFFICIENTS variable of the ProcessInfo container.
     */
    void SetTimeCoefficients()
    {
        KRATOS_TRY;

        auto &r_process_info = (BaseType::GetModelPart()).GetProcessInfo();

        if (mTimeOrder == 2)
        {
            //calculate the BDF coefficients
            double Dt = r_process_info[DELTA_TIME];
            double OldDt = r_process_info.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            double Rho = OldDt / Dt;
            double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

            Vector& BDFcoeffs = r_process_info[BDF_COEFFICIENTS];
            BDFcoeffs.resize(3, false);

            BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
            BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
            BDFcoeffs[2] = TimeCoeff; //coefficient for step n-1 (1/2Dt if Dt is constant)
        }
        else if (mTimeOrder == 1)
        {
            double Dt = r_process_info[DELTA_TIME];
            double TimeCoeff = 1.0 / Dt;

            Vector& BDFcoeffs = r_process_info[BDF_COEFFICIENTS];
            BDFcoeffs.resize(2, false);

            BDFcoeffs[0] = TimeCoeff; //coefficient for step n+1 (1/Dt)
            BDFcoeffs[1] = -TimeCoeff; //coefficient for step n (-1/Dt)
        }

        KRATOS_CATCH("");
    }

    virtual std::tuple<bool,double> SolveStep()
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        const int n_nodes = rModelPart.NumberOfNodes();

        // 1. Fractional step momentum iteration
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);

        bool Converged = false;
        for(unsigned int it = 0; it < mMaxVelocityIter; ++it)
        {
            KRATOS_INFO_IF("FractionalStepStrategy", BaseType::GetEchoLevel() > 1) << "Momentum iteration " << it << std::endl;

            // build momentum system and solve for fractional step velocity increment
            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);
            double NormDv = mpMomentumStrategy->Solve();

//            // Compute projections (for stabilization)
//            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
//            this->ComputeSplitOssProjections(rModelPart);

//            // Additional steps // Moved to end of step
//            for (std::vector<Process::Pointer>::iterator iExtraSteps = mExtraIterationSteps.begin();
//                 iExtraSteps != mExtraIterationSteps.end(); ++iExtraSteps)
//                (*iExtraSteps)->Execute();

            // Check convergence
            Converged = this->CheckFractionalStepConvergence(NormDv);

            if (Converged)
            {
                KRATOS_INFO_IF("FractionalStepStrategy", BaseType::GetEchoLevel() > 0) << "Fractional velocity converged in " << it + 1 << " iterations." << std::endl;
                break;
            }
        }

        KRATOS_INFO_IF("FractionalStepStrategy", !Converged && BaseType::GetEchoLevel() > 0) << "Fractional velocity iterations did not converge." << std::endl;

        // Compute projections (for stabilization)
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
        this->ComputeSplitOssProjections(rModelPart);

        // 2. Pressure solution (store pressure variation in PRESSURE_OLD_IT)
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,5);

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = rModelPart.NodesBegin() + i_node;
            const double old_press = it_node->FastGetSolutionStepValue(PRESSURE);
            it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) = -old_press;
        }

        KRATOS_INFO_IF("FractionalStepStrategy", BaseType::GetEchoLevel() > 0) << "Calculating Pressure." << std::endl;
        double NormDp = mpPressureStrategy->Solve();

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = rModelPart.NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) += it_node->FastGetSolutionStepValue(PRESSURE);
        }

        // 3. Compute end-of-step velocity
        KRATOS_INFO_IF("FractionalStepStrategy", BaseType::GetEchoLevel() > 0) << "Updating Velocity." << std::endl;
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,6);

        this->CalculateEndOfStepVelocity();

        /*
        mpPressureStrategy->Clear();
        double NormDu = mpPressureStrategy->Solve();
        mpPressureStrategy->Clear();
        */

        // Additional steps
        for (std::vector<Process::Pointer>::iterator iExtraSteps = mExtraIterationSteps.begin();
             iExtraSteps != mExtraIterationSteps.end(); ++iExtraSteps)
            (*iExtraSteps)->Execute();

        // Set the output tuple as the fractional velocity convergence and pressure norm
        return std::make_tuple(Converged, NormDp);
    }

    bool CheckFractionalStepConvergence(const double NormDv)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        const int n_nodes = rModelPart.NumberOfNodes();

        double NormV = 0.00;
#pragma omp parallel for reduction(+:NormV)
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            const auto it_node = rModelPart.NodesBegin() + i_node;
            const auto &r_vel = it_node->FastGetSolutionStepValue(VELOCITY);
            for (unsigned int d = 0; d < 3; ++d) {
                NormV += r_vel[d] * r_vel[d];
            }
        }
        NormV = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormV);
        NormV = sqrt(NormV);

        const double zero_tol = 1.0e-12;
        const double Ratio = (NormV < zero_tol) ? NormDv : NormDv / NormV;

        KRATOS_INFO_IF("FractionalStepStrategy", BaseType::GetEchoLevel() > 0) << "CONVERGENCE CHECK:" << std::endl;
        KRATOS_INFO_IF("FractionalStepStrategy", BaseType::GetEchoLevel() > 0)
            <<  std::scientific << std::setprecision(8) << "FRAC VEL.: ratio = " << Ratio <<"; exp.ratio = " << mVelocityTolerance << " abs = " << NormDv << std::endl;

        if (Ratio < mVelocityTolerance)
            return true;
        else
            return false;
    }

    bool CheckPressureConvergence(const double NormDp)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        const int n_nodes = rModelPart.NumberOfNodes();

        double NormP = 0.00;
#pragma omp parallel for reduction(+:NormP)
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            const auto it_node = rModelPart.NodesBegin() + i_node;
            const double Pr = it_node->FastGetSolutionStepValue(PRESSURE);
            NormP += Pr * Pr;
        }
        NormP = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormP);
        NormP = sqrt(NormP);

        const double zero_tol = 1.0e-12;
        const double Ratio = (NormP < zero_tol) ? NormDp : NormDp / NormP;

        KRATOS_INFO_IF("FractionalStepStrategy", BaseType::GetEchoLevel() > 0) << "Pressure relative error: " << Ratio << std::endl;

        if (Ratio < mPressureTolerance)
        {
            return true;
        }
        else
            return false;
    }


    virtual void ComputeSplitOssProjections(ModelPart& rModelPart)
    {
        array_1d<double,3> Out = ZeroVector(3);
        const int n_nodes = rModelPart.NumberOfNodes();
        const int n_elems = rModelPart.NumberOfElements();

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = rModelPart.NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(CONV_PROJ) = CONV_PROJ.Zero();
            it_node->FastGetSolutionStepValue(PRESS_PROJ) = PRESS_PROJ.Zero();
            it_node->FastGetSolutionStepValue(DIVPROJ) = 0.0;
            it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        }

#pragma omp parallel for
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            const auto it_elem = rModelPart.ElementsBegin() + i_elem;
            it_elem->Calculate(CONV_PROJ, Out, rModelPart.GetProcessInfo());
        }

        rModelPart.GetCommunicator().AssembleCurrentData(CONV_PROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(PRESS_PROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);

        // If there are periodic conditions, add contributions from both sides to the periodic nodes
        this->PeriodicConditionProjectionCorrection(rModelPart);

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = rModelPart.NodesBegin() + i_node;
            const double NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);
            it_node->FastGetSolutionStepValue(CONV_PROJ) /= NodalArea;
            it_node->FastGetSolutionStepValue(PRESS_PROJ) /= NodalArea;
            it_node->FastGetSolutionStepValue(DIVPROJ) /= NodalArea;
        }
    }

    virtual void CalculateEndOfStepVelocity()
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        const int n_nodes = rModelPart.NumberOfNodes();
        const int n_elems = rModelPart.NumberOfElements();

        array_1d<double,3> Out = ZeroVector(3);
        VariableUtils().SetHistoricalVariableToZero(FRACT_VEL, rModelPart.Nodes());

#pragma omp parallel for
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            const auto it_elem = rModelPart.ElementsBegin() + i_elem;
            it_elem->Calculate(VELOCITY, Out, rModelPart.GetProcessInfo());
        }

        rModelPart.GetCommunicator().AssembleCurrentData(FRACT_VEL);
        this->PeriodicConditionVelocityCorrection(rModelPart);

        // Force the end of step velocity to verify slip conditions in the model
        if (mUseSlipConditions)
            this->EnforceSlipCondition(SLIP);

        if (mDomainSize > 2)
        {
#pragma omp parallel for
            for (int i_node = 0; i_node < n_nodes; ++i_node) {
                auto it_node = rModelPart.NodesBegin() + i_node;
                const double NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);
                if ( ! it_node->IsFixed(VELOCITY_X) )
                    it_node->FastGetSolutionStepValue(VELOCITY_X) += it_node->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                if ( ! it_node->IsFixed(VELOCITY_Y) )
                    it_node->FastGetSolutionStepValue(VELOCITY_Y) += it_node->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
                if ( ! it_node->IsFixed(VELOCITY_Z) )
                    it_node->FastGetSolutionStepValue(VELOCITY_Z) += it_node->FastGetSolutionStepValue(FRACT_VEL_Z) / NodalArea;
            }
        }
        else
        {
#pragma omp parallel for
            for (int i_node = 0; i_node < n_nodes; ++i_node) {
                auto it_node = rModelPart.NodesBegin() + i_node;
                const double NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);
                if ( ! it_node->IsFixed(VELOCITY_X) )
                    it_node->FastGetSolutionStepValue(VELOCITY_X) += it_node->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                if ( ! it_node->IsFixed(VELOCITY_Y) )
                    it_node->FastGetSolutionStepValue(VELOCITY_Y) += it_node->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
            }
        }
    }

    /**
     * @brief Substract wall-normal component of velocity update to ensure that the final velocity satisfies slip conditions.
     * @param rSlipWallFlag If Node.Is(rSlipWallFlag) == true, the node is in the wall.
     */
    void EnforceSlipCondition(const Kratos::Flags& rSlipWallFlag)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        const int num_nodes_in_model_part = rModelPart.NumberOfNodes();

        #pragma omp parallel for
        for (int i = 0; i < num_nodes_in_model_part; i++)
        {
            ModelPart::NodeIterator itNode = rModelPart.NodesBegin() + i;
            const Node<3>& r_const_node = *itNode;

            if ( r_const_node.Is(rSlipWallFlag) )
            {
                const array_1d<double,3>& rNormal = itNode->FastGetSolutionStepValue(NORMAL);
                array_1d<double,3>& rDeltaVelocity = itNode->FastGetSolutionStepValue(FRACT_VEL);

                double Proj = rNormal[0] * rDeltaVelocity[0];
                double Norm = rNormal[0] * rNormal[0];

                for (unsigned int d = 1; d < mDomainSize; ++d)
                {
                    Proj += rNormal[d] * rDeltaVelocity[d];
                    Norm += rNormal[d] * rNormal[d];
                }

                Proj /= Norm;
                rDeltaVelocity -= Proj * rNormal;
            }
        }
    }

    /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
     * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
     * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
     * 2- The non-historical containers are added across processes, transmiting the right value from the condition owner to all partitions.\n
     * 3- The value on all periodic nodes is replaced by the one received in step 2.
     */
     void PeriodicConditionProjectionCorrection(ModelPart& rModelPart)
     {
         Communicator& r_comm = rModelPart.GetCommunicator();
         if (mrPeriodicIdVar.Key() != Kratos::Variable<int>::StaticObject().Key())
         {
             int GlobalNodesNum = r_comm.LocalMesh().Nodes().size();
             GlobalNodesNum = r_comm.GetDataCommunicator().SumAll(GlobalNodesNum);

             for (typename ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); itCond++ )
             {
                 ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
                 if (rGeom.PointsNumber() == 2)
                 {
                     Node<3>& rNode0 = rGeom[0];
                     int Node0Pair = rNode0.FastGetSolutionStepValue(mrPeriodicIdVar);

                     Node<3>& rNode1 = rGeom[1];
                     int Node1Pair = rNode1.FastGetSolutionStepValue(mrPeriodicIdVar);

                     // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                     if ( ( static_cast<int>(rNode0.Id()) == Node1Pair ) && (static_cast<int>(rNode1.Id()) == Node0Pair ) )
                     {
                         double NodalArea = rNode0.FastGetSolutionStepValue(NODAL_AREA) + rNode1.FastGetSolutionStepValue(NODAL_AREA);
                         array_1d<double,3> ConvProj = rNode0.FastGetSolutionStepValue(CONV_PROJ) + rNode1.FastGetSolutionStepValue(CONV_PROJ);
                         array_1d<double,3> PressProj = rNode0.FastGetSolutionStepValue(PRESS_PROJ) + rNode1.FastGetSolutionStepValue(PRESS_PROJ);
                         double DivProj = rNode0.FastGetSolutionStepValue(DIVPROJ) + rNode1.FastGetSolutionStepValue(DIVPROJ);

                         rNode0.GetValue(NODAL_AREA) = NodalArea;
                         rNode0.GetValue(CONV_PROJ) = ConvProj;
                         rNode0.GetValue(PRESS_PROJ) = PressProj;
                         rNode0.GetValue(DIVPROJ) = DivProj;
                         rNode1.GetValue(NODAL_AREA) = NodalArea;
                         rNode1.GetValue(CONV_PROJ) = ConvProj;
                         rNode1.GetValue(PRESS_PROJ) = PressProj;
                         rNode1.GetValue(DIVPROJ) = DivProj;
                     }
                 }
                 else if (rGeom.PointsNumber() == 4 && rGeom[0].FastGetSolutionStepValue(mrPeriodicIdVar) > GlobalNodesNum)
                 {
                     double NodalArea = rGeom[0].FastGetSolutionStepValue(NODAL_AREA);
                     array_1d<double,3> ConvProj = rGeom[0].FastGetSolutionStepValue(CONV_PROJ);
                     array_1d<double,3> PressProj = rGeom[0].FastGetSolutionStepValue(PRESS_PROJ);
                     double DivProj = rGeom[0].FastGetSolutionStepValue(DIVPROJ);

                     for (unsigned int i = 1; i < 4; i++)
                     {
                         NodalArea += rGeom[i].FastGetSolutionStepValue(NODAL_AREA);
                         ConvProj += rGeom[i].FastGetSolutionStepValue(CONV_PROJ);
                         PressProj += rGeom[i].FastGetSolutionStepValue(PRESS_PROJ);
                         DivProj += rGeom[i].FastGetSolutionStepValue(DIVPROJ);
                     }

                     for (unsigned int i = 0; i < 4; i++)
                     {
                         rGeom[i].GetValue(NODAL_AREA) = NodalArea;
                         rGeom[i].GetValue(CONV_PROJ) = ConvProj;
                         rGeom[i].GetValue(PRESS_PROJ) = PressProj;
                         rGeom[i].GetValue(DIVPROJ) = DivProj;
                     }
                 }
             }

             rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
             rModelPart.GetCommunicator().AssembleNonHistoricalData(CONV_PROJ);
             rModelPart.GetCommunicator().AssembleNonHistoricalData(PRESS_PROJ);
             rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

             for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
             {
                 if (itNode->GetValue(NODAL_AREA) != 0.0)
                 {
                     itNode->FastGetSolutionStepValue(NODAL_AREA) = itNode->GetValue(NODAL_AREA);
                     itNode->FastGetSolutionStepValue(CONV_PROJ) = itNode->GetValue(CONV_PROJ);
                     itNode->FastGetSolutionStepValue(PRESS_PROJ) = itNode->GetValue(PRESS_PROJ);
                     itNode->FastGetSolutionStepValue(DIVPROJ) = itNode->GetValue(DIVPROJ);

                     // reset for next iteration
                     itNode->GetValue(NODAL_AREA) = 0.0;
                     itNode->GetValue(CONV_PROJ) = CONV_PROJ.Zero();
                     itNode->GetValue(PRESS_PROJ) = PRESS_PROJ.Zero();
                     itNode->GetValue(DIVPROJ) = 0.0;
                 }
             }
         }
     }

     void PeriodicConditionVelocityCorrection(ModelPart& rModelPart)
     {
         Communicator& r_comm = rModelPart.GetCommunicator();
         if (mrPeriodicIdVar.Key() != Kratos::Variable<int>::StaticObject().Key())
         {
             int GlobalNodesNum = r_comm.LocalMesh().Nodes().size();
             GlobalNodesNum = r_comm.GetDataCommunicator().SumAll(GlobalNodesNum);

             for (typename ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); itCond++ )
             {
                 ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
                 if (rGeom.PointsNumber() == 2)
                 {
                     Node<3>& rNode0 = rGeom[0];
                     int Node0Pair = rNode0.FastGetSolutionStepValue(mrPeriodicIdVar);

                     Node<3>& rNode1 = rGeom[1];
                     int Node1Pair = rNode1.FastGetSolutionStepValue(mrPeriodicIdVar);

                     // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                     if ( ( static_cast<int>(rNode0.Id()) == Node1Pair ) && (static_cast<int>(rNode1.Id()) == Node0Pair ) )
                     {
                         array_1d<double,3> DeltaVel = rNode0.FastGetSolutionStepValue(FRACT_VEL) + rNode1.FastGetSolutionStepValue(FRACT_VEL);

                         rNode0.GetValue(FRACT_VEL) = DeltaVel;
                         rNode1.GetValue(FRACT_VEL) = DeltaVel;
                     }
                 }
                 else if (rGeom.PointsNumber() == 4 && rGeom[0].FastGetSolutionStepValue(mrPeriodicIdVar) > GlobalNodesNum)
                 {
                     array_1d<double,3> DeltaVel = rGeom[0].FastGetSolutionStepValue(FRACT_VEL);
                     for (unsigned int i = 1; i < 4; i++)
                     {
                         DeltaVel += rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
                     }

                     for (unsigned int i = 0; i < 4; i++)
                     {
                         rGeom[i].GetValue(FRACT_VEL) = DeltaVel;
                     }
                 }
             }

             rModelPart.GetCommunicator().AssembleNonHistoricalData(FRACT_VEL);

             for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
             {
                 array_1d<double,3>& rDeltaVel = itNode->GetValue(FRACT_VEL);
                 if ( rDeltaVel[0]*rDeltaVel[0] + rDeltaVel[1]*rDeltaVel[1] + rDeltaVel[2]*rDeltaVel[2] != 0.0)
                 {
                     itNode->FastGetSolutionStepValue(FRACT_VEL) = itNode->GetValue(FRACT_VEL);
                     rDeltaVel = ZeroVector(3);
                 }
             }
         }
     }


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void InitializeStrategy(
        SolverSettingsType& rSolverConfig,
        bool PredictorCorrector)
    {
        KRATOS_TRY;

        mTimeOrder = rSolverConfig.GetTimeOrder();

        // Check that input parameters are reasonable and sufficient.
        this->Check();

        mDomainSize = rSolverConfig.GetDomainSize();

        mPredictorCorrector = PredictorCorrector;

        mUseSlipConditions = rSolverConfig.UseSlipConditions();

        mReformDofSet = rSolverConfig.GetReformDofSet();

        BaseType::SetEchoLevel(rSolverConfig.GetEchoLevel());

        // Initialize strategies for each step
        bool HaveVelStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Velocity,mpMomentumStrategy);

        if (HaveVelStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Velocity,mVelocityTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Velocity,mMaxVelocityIter);
        }
        else
        {
            KRATOS_ERROR << "FractionalStepStrategy error: No Velocity strategy defined in FractionalStepSettings" << std::endl;
        }

        bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure,mpPressureStrategy);

        if (HavePressStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Pressure,mPressureTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Pressure,mMaxPressureIter);
        }
        else
        {
            KRATOS_ERROR << "FractionalStepStrategy error: No Pressure strategy defined in FractionalStepSettings" << std::endl;
        }

        Process::Pointer pTurbulenceProcess;
        bool HaveTurbulence = rSolverConfig.GetTurbulenceModel(pTurbulenceProcess);

        if (HaveTurbulence)
            mExtraIterationSteps.push_back(pTurbulenceProcess);


        // Check input parameters
        this->Check();

        KRATOS_CATCH("");
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FractionalStepStrategy& operator=(FractionalStepStrategy const& rOther){}

    /// Copy constructor.
    FractionalStepStrategy(FractionalStepStrategy const& rOther){}


    ///@}

}; /// Class FStepStrategy

///@}
///@name Type Definitions
///@{


///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_FRACTIONAL_STEP_STRATEGY
