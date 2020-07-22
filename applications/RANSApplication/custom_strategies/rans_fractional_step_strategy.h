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
//  Extended by :    Suneth Warnakulasuriya

#if !defined(KRATOS_RANS_FRACTIONAL_STEP_STRATEGY)
#define KRATOS_RANS_FRACTIONAL_STEP_STRATEGY

// System includes

// External includes

// Project includes
#include "custom_strategies/strategies/fractional_step_strategy.h"

// Application includes

namespace Kratos
{
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
 * This strategy implements a splitting scheme for the incompressible
 * Navier-Stokes equations. It is intended to be used in combination with the
 * FractionalStep element in the FluidDynamicsApplicatoin. The fractional step
 * index, which is stored in the ProcessInfo, takes the values 1 : Momentum step
 * (calculate fractional step velocity) 2-3 : Unused (reserved for componentwise
 * calculation of frac step velocity) 4 : Pressure step 5 : Computation of
 * projections 6 : End of step velocity
 * @tparam TSparseSpace Sparse space template type
 * @tparam TDenseSpace Dense space template type
 * @tparam TLinearSolver Linear solver template type
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class RansFractionalStepStrategy
    : public FractionalStepStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of RansFractionalStepStrategy
    KRATOS_CLASS_POINTER_DEFINITION(RansFractionalStepStrategy);

    typedef FractionalStepStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef SolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> SolverSettingsType;

    ///@}
    ///@name Life Cycle
    ///@{

    RansFractionalStepStrategy(
        ModelPart& rModelPart,
        SolverSettingsType& rSolverConfig,
        bool PredictorCorrector,
        bool CalculateReactionsFlag)
    : BaseType(rModelPart, rSolverConfig, PredictorCorrector, CalculateReactionsFlag)
    {
        KRATOS_INFO(this->Info()) << "Created fractional step strategy." << std::endl;
    }

    RansFractionalStepStrategy(
        ModelPart& rModelPart,
        SolverSettingsType& rSolverConfig,
        bool PredictorCorrector,
        bool CalculateReactionsFlag,
        const Kratos::Variable<int>& PeriodicVar)
    : BaseType(rModelPart, rSolverConfig, PredictorCorrector, CalculateReactionsFlag, PeriodicVar)
    {
        KRATOS_INFO(this->Info())
            << "Created periodic fractional step strategy." << std::endl;
    }

    /// Destructor.
    ~RansFractionalStepStrategy() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
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
        std::stringstream buffer;
        buffer << "RansFractionalStepStrategy";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    std::tuple<bool, double> SolveStep() override
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        // 1. Fractional step momentum iteration
        r_process_info.SetValue(FRACTIONAL_STEP, 1);

        bool converged = false;
        for (unsigned int it = 0; it < this->mMaxVelocityIter; ++it) {
            KRATOS_INFO_IF("RansFractionalStepStrategy", BaseType::GetEchoLevel() > 1)
                << "Momentum iteration " << it << std::endl;

            // build momentum system and solve for fractional step velocity increment
            r_process_info.SetValue(FRACTIONAL_STEP, 1);
            double norm_dv = this->mpMomentumStrategy->Solve();

            // Check convergence
            converged = this->CheckFractionalStepConvergence(norm_dv);

            if (converged) {
                KRATOS_INFO_IF("RansFractionalStepStrategy", BaseType::GetEchoLevel() > 0)
                    << "Fractional velocity converged in " << it + 1
                    << " iterations." << std::endl;
                break;
            }
        }

        KRATOS_INFO_IF("RansFractionalStepStrategy", !converged && BaseType::GetEchoLevel() > 0)
            << "Fractional velocity iterations did not converge." << std::endl;

        // Compute projections (for stabilization)
        r_process_info.SetValue(FRACTIONAL_STEP, 4);
        this->ComputeSplitOssProjections(r_model_part);

        // 2. Pressure solution (store pressure variation in PRESSURE_OLD_IT)
        r_process_info.SetValue(FRACTIONAL_STEP, 5);

        const double eta = (r_process_info.Has(PRESSURE_COEFFICIENT))
                               ? r_process_info[PRESSURE_COEFFICIENT]
                               : 1.0;

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            const double old_press = it_node->FastGetSolutionStepValue(PRESSURE);
            it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) = -eta * old_press;
        }

        KRATOS_INFO_IF("RansFractionalStepStrategy", BaseType::GetEchoLevel() > 0)
            << "Calculating Pressure." << std::endl;
        double NormDp = this->mpPressureStrategy->Solve();

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) +=
                it_node->FastGetSolutionStepValue(PRESSURE);
        }

        // 3. Compute end-of-step velocity
        KRATOS_INFO_IF("RansFractionalStepStrategy", BaseType::GetEchoLevel() > 0)
            << "Updating Velocity." << std::endl;
        r_process_info.SetValue(FRACTIONAL_STEP, 6);

        this->CalculateEndOfStepVelocity();

        // Set the output tuple as the fractional velocity convergence and pressure norm
        return std::make_tuple(converged, NormDp);
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
    RansFractionalStepStrategy& operator=(RansFractionalStepStrategy const& rOther)
    {
    }

    /// Copy constructor.
    RansFractionalStepStrategy(RansFractionalStepStrategy const& rOther)
    {
    }

    ///@}

}; /// Class FStepStrategy

///@}
///@name Type Definitions
///@{

///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_RANS_FRACTIONAL_STEP_STRATEGY
