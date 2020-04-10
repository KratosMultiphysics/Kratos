//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_COUPLED_STRATEGY_H)
#define KRATOS_RANS_COUPLED_STRATEGY_H

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "coupled_strategy_item.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@addtogroup RANSApplication
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

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class CoupledStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
    struct ConvergenceVariableSettings
    {
    public:
        std::string Name;
        double RelativeTolerance;
        double AbsoluteTolerance;
        double RelativeError = 1.0;
        double AbsoluteError = 1.0;
        bool IsConverged;
    };

public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of CoupledStrategy
    KRATOS_CLASS_POINTER_DEFINITION(CoupledStrategy);

    using BaseType = SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using CouplingStrategyItemType =
        CoupledStrategyItem<TSparseSpace, TDenseSpace, TLinearSolver>;

    ///@}
    ///@name Life Cycle
    ///@{

    CoupledStrategy(ModelPart& rModelPart,
                    const bool CheckConvergenceCriteriaConvergence,
                    const bool CheckVariableTransientConvergence,
                    const bool MoveMesh = false,
                    const int MaxIterations = 10)
        : BaseType(rModelPart, MoveMesh),
          mCheckConvergenceCriteriaConvergence(CheckConvergenceCriteriaConvergence),
          mCheckVariableTransientConvergence(CheckVariableTransientConvergence),
          mMaxIterations(MaxIterations)
    {
        if (mCheckVariableTransientConvergence)
        {
            const auto buffer_size = rModelPart.GetBufferSize();
            KRATOS_ERROR_IF(buffer_size < 2)
                << "CoupledStrategy with transient variable convergence check "
                   "requires buffer size greater than 1 in "
                << rModelPart.Name() << ". [ buffer_size = " << buffer_size << " ].\n";
        }
    }

    /// Destructor.
    ~CoupledStrategy() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddStrategyItem(typename CouplingStrategyItemType::Pointer pStrategyItem)
    {
        mStrategiesList.push_back(pStrategyItem);
    }

    void AddConvergenceCheckVariable(const std::string& rVariableName,
                                     const double RelativeTolerance = 1e-3,
                                     const double AbsoluteTolerance = 1e-5)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF((!KratosComponents<Variable<double>>::Has(rVariableName) &&
                         !KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)))
            << "Only double or 3d variables are supported as convergence "
               "variables. "
            << rVariableName << " is not found in either lists.\n";

        ConvergenceVariableSettings current_variable;
        current_variable.Name = rVariableName;
        current_variable.RelativeTolerance = RelativeTolerance;
        current_variable.AbsoluteTolerance = AbsoluteTolerance;

        mConvergenceVariableSettingsList.push_back(current_variable);

        KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 0)
            << "Added " << rVariableName << " transient convergence variable\n";

        KRATOS_CATCH("");
    }

    void Initialize() override
    {
        KRATOS_TRY

        for (auto& p_strategy : mStrategiesList)
            for (auto& p_process : p_strategy->GetAuxiliaryProcessList())
                p_process->ExecuteInitialize();

        std::stringstream buffer;
        buffer << "Strategies list: \n";
        for (auto& strategy : mStrategiesList)
        {
            strategy->GetStrategy().Initialize();
            buffer << strategy->GetStrategyInfo();
        }
        buffer << "  Convergence settings:\n";
        buffer << "      Check convergence criteria: "
               << (mCheckConvergenceCriteriaConvergence ? "yes" : "no") << "\n";
        buffer << "      Check transinet variables : "
               << (mCheckVariableTransientConvergence ? "yes" : "no") << "\n";

        if (mCheckVariableTransientConvergence)
        {
            buffer << "      Variables:\n";
            for (const auto& variable_settings : mConvergenceVariableSettingsList)
            {
                buffer << "          Name: " << variable_settings.Name << "\n";
                buffer << "             relative tolerance: "
                       << variable_settings.RelativeTolerance << "\n";
                buffer << "             absolute tolerance: "
                       << variable_settings.AbsoluteTolerance << "\n";
            }
        }

        KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 0) << buffer.str();

        KRATOS_CATCH("");
    }

    int Check() override
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(mStrategiesList.size() == 0)
            << "No strategies are found for CoupledStrategy.";

        KRATOS_ERROR_IF((mCheckVariableTransientConvergence &&
                         mConvergenceVariableSettingsList.size() == 0))
            << "Transient variable convergence check requires non-empty "
               "convergence check variable list.\n";

        for (auto& p_strategy : mStrategiesList)
            for (auto& p_process : p_strategy->GetAuxiliaryProcessList())
                p_process->Check();

        for (auto& strategy_item : mStrategiesList)
            strategy_item->GetStrategy().Check();

        return 0;

        KRATOS_CATCH("");
    }

    bool IsConverged() override
    {
        KRATOS_TRY

        bool is_converged = true;

        if (mCheckConvergenceCriteriaConvergence)
        {
            for (auto& p_strategy_item : mStrategiesList)
            {
                const bool is_current_strategy_converged =
                    p_strategy_item->GetStrategy().IsConverged();
                is_converged = (is_converged) ? is_current_strategy_converged : is_converged;

                KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 1 && is_current_strategy_converged)
                    << "*** CONVERGENCE ACHIEVED FOR COUPLED STRATEGY *** : "
                    << p_strategy_item->GetName() << std::endl;
            }
        }

        if (mCheckVariableTransientConvergence)
        {
            const ModelPart& r_model_part = this->GetModelPart();
            for (ConvergenceVariableSettings& r_variable_settings : mConvergenceVariableSettingsList)
            {
                double relative_error, absolute_error;
                if (KratosComponents<Variable<double>>::Has(r_variable_settings.Name))
                {
                    const Variable<double>& r_variable =
                        KratosComponents<Variable<double>>::Get(
                            r_variable_settings.Name);

                    RansVariableUtilities::CalculateTransientVariableConvergence(
                        relative_error, absolute_error, r_model_part, r_variable);
                }
                else
                {
                    const Variable<array_1d<double, 3>>& r_variable =
                        KratosComponents<Variable<array_1d<double, 3>>>::Get(
                            r_variable_settings.Name);

                    RansVariableUtilities::CalculateTransientVariableConvergence(
                        relative_error, absolute_error, r_model_part, r_variable);
                }

                if (relative_error != 0.0 && absolute_error != 0.0)
                {
                    r_variable_settings.RelativeError = relative_error;
                    r_variable_settings.AbsoluteError = absolute_error;
                }

                r_variable_settings.IsConverged =
                    (r_variable_settings.RelativeError < r_variable_settings.RelativeTolerance ||
                     r_variable_settings.AbsoluteError < r_variable_settings.AbsoluteTolerance);

                is_converged = (is_converged) ? r_variable_settings.IsConverged : is_converged;

                if (this->GetEchoLevel() > 1)
                {
                    std::stringstream buffer;
                    buffer << std::scientific << std::setprecision(6)
                           << "[ Obtained ratio: " << r_variable_settings.RelativeError
                           << "; Expected ratio: " << r_variable_settings.RelativeTolerance
                           << "; Absolute norm: " << r_variable_settings.AbsoluteError
                           << "; Expected norm: " << r_variable_settings.AbsoluteTolerance
                           << " ] - " << r_variable_settings.Name << std::endl;

                    KRATOS_INFO(this->Info()) << buffer.str();
                }
            }
            std::stringstream buffer;
            const auto& r_variable_settings =
                *(mConvergenceVariableSettingsList.begin());
            buffer << (r_variable_settings.IsConverged ? r_variable_settings.Name : "");
            for (int i = 1;
                 i < static_cast<int>(mConvergenceVariableSettingsList.size()); ++i)
            {
                const auto& r_variable_settings =
                    *(mConvergenceVariableSettingsList.begin() + i);
                buffer << (r_variable_settings.IsConverged
                               ? (", " + r_variable_settings.Name)
                               : "");
            }
            KRATOS_INFO_IF(this->Info(), (this->GetEchoLevel() > 1 && buffer.str() != ""))
                << buffer.str() << " *** TRANSIENT CONVERGENCE ACHIEVED ***"
                << std::endl;
        }

        KRATOS_INFO_IF(this->Info(), (is_converged && this->GetEchoLevel() > 1))
            << "*** CONVERGENCE ACHIEVED ***" << std::endl;

        return is_converged;

        KRATOS_CATCH("");
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

        const int number_of_solving_strategies = this->mStrategiesList.size();

        bool is_converged = false;
        int& iteration = r_current_process_info[COUPLING_ITERATION];
        iteration = 1;

        while (!is_converged && iteration <= this->mMaxIterations)
        {
            for (int i = 0; i < number_of_solving_strategies; ++i)
            {
                auto& p_solving_strategy = this->mStrategiesList[i];
                if (p_solving_strategy->IsStrategySolvable())
                {
                    p_solving_strategy->GetStrategy().SolveSolutionStep();

                    // execute update processes
                    for (auto& p_process : p_solving_strategy->GetAuxiliaryProcessList())
                        p_process->Execute();

                    KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 1)
                        << "Executed update processes for "
                        << p_solving_strategy->GetName() << ".\n";

                    const unsigned int iterations =
                        r_current_process_info[NL_ITERATION_NUMBER];
                    KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 1)
                        << "Solving " << p_solving_strategy->GetName()
                        << " used " << iterations << " non-linear iterations.\n";
                }
            }

            is_converged = this->IsConverged();
            // Check for all the strategy convergence
            if (is_converged)
            {
                for (int i = 0; i < number_of_solving_strategies; ++i)
                {
                    auto& p_solving_strategy = this->mStrategiesList[i];
                    p_solving_strategy->GetStrategy().SolveSolutionStep();

                    // execute update processes
                    for (auto& p_process : p_solving_strategy->GetAuxiliaryProcessList())
                        p_process->Execute();
                }
                is_converged = this->IsConverged();
            }
            ++iteration;
        }

        return is_converged;

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        for (auto& p_strategy : mStrategiesList)
            for (auto& p_process : p_strategy->GetAuxiliaryProcessList())
                p_process->ExecuteInitializeSolutionStep();

        for (auto& p_solving_strategy : this->mStrategiesList)
            p_solving_strategy->GetStrategy().InitializeSolutionStep();

        KRATOS_CATCH("");
    }

    void Predict() override
    {
        KRATOS_TRY

        for (auto& p_solving_strategy : this->mStrategiesList)
            p_solving_strategy->GetStrategy().Predict();

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        for (auto& p_solving_strategy : this->mStrategiesList)
            p_solving_strategy->GetStrategy().FinalizeSolutionStep();

        for (auto& p_strategy : mStrategiesList)
            for (auto& p_process : p_strategy->GetAuxiliaryProcessList())
                p_process->ExecuteFinalizeSolutionStep();

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        KRATOS_TRY

        for (auto& p_solving_strategy : this->mStrategiesList)
            p_solving_strategy->GetStrategy().Clear();

        KRATOS_CATCH("");
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
        std::stringstream buffer;
        buffer << "CoupledStrategy";
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

    bool mCheckConvergenceCriteriaConvergence;
    bool mCheckVariableTransientConvergence;
    int mMaxIterations;

    std::vector<typename CouplingStrategyItemType::Pointer> mStrategiesList;
    std::vector<ConvergenceVariableSettings> mConvergenceVariableSettingsList;

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
    CoupledStrategy& operator=(CoupledStrategy const& rOther)
    {
    }

    /// Copy constructor.
    CoupledStrategy(CoupledStrategy const& rOther)
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

#endif // KRATOS_RANS_COUPLED_STRATEGY_H
