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

#if !defined(SCALAR_CO_SOLVING_PROCESS_H_INCLUDED)
#define SCALAR_CO_SOLVING_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "custom_utilities/rans_variable_utils.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// The base class for all ScalarCoSolvingProcesses in Kratos.
/** The ScalarCoSolvingProcess is the base class for all ScalarCoSolvingProcesses and defines a simple interface for them.
    Execute method is used to execute the ScalarCoSolvingProcess algorithms. While the parameters of this method
  can be very different from one ScalarCoSolvingProcess to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all ScalarCoSolvingProcess parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other ScalarCoSolvingProcesses or the base ScalarCoSolvingProcess class.
*/
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ScalarCoSolvingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> SolvingStrategyType;

    /// Pointer definition of ScalarCoSolvingProcess
    KRATOS_CLASS_POINTER_DEFINITION(ScalarCoSolvingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ScalarCoSolvingProcess(ModelPart& rModelPart,
                           Parameters& rParameters,
                           Variable<double>& rConvergenceVariable)
        : Process(), mrModelPart(rModelPart), mrConvergenceVariable(rConvergenceVariable)
    {
        Parameters default_parameters = Parameters(R"(
        {
            "relative_tolerance"    : 1e-3,
            "absolute_tolerance"    : 1e-5,
            "max_iterations"        : 10,
            "echo_level"            : 0
        })");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = rParameters["echo_level"].GetInt();
        mConvergenceRelativeTolerance = rParameters["relative_tolerance"].GetDouble();
        mConvergenceAbsoluteTolerance = rParameters["absolute_tolerance"].GetDouble();
        mMaxIterations = rParameters["max_iterations"].GetInt();
    }

    /// Destructor.
    ~ScalarCoSolvingProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddStrategy(SolvingStrategyType* pStrategy)
    {
        mrSolvingStrategiesList.push_back(pStrategy);
    }

    /// Execute method is used to execute the ScalarCoSolvingProcess algorithms.
    void Execute() override
    {
        if (mrModelPart.GetProcessInfo()[IS_CO_SOLVING_PROCESS_ACTIVE])
            SolveSolutionStep();
    }

    virtual int Check() override
    {
        KRATOS_CHECK_VARIABLE_KEY(IS_CO_SOLVING_PROCESS_ACTIVE);
        KRATOS_CHECK_VARIABLE_KEY(OLD_CONVERGENCE_VARIABLE);
        KRATOS_CHECK_VARIABLE_KEY(this->mrConvergenceVariable);

        for (ModelPart::NodeType& r_node : mrModelPart.Nodes())
        {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(OLD_CONVERGENCE_VARIABLE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(this->mrConvergenceVariable, r_node);
        }

        for (SolvingStrategyType* strategy : mrSolvingStrategiesList)
            strategy->Check();

        KRATOS_ERROR_IF(mrSolvingStrategiesList.size() == 0)
            << "No strategies are found for ScalarCoSolvingProcess.";

        return 0;
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
    virtual std::string Info() const override
    {
        return "ScalarCoSolvingProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ScalarCoSolvingProcess";
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
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    int mEchoLevel;
    bool mIsCoSolvingProcessActive;
    ///@}
    ///@name Operations
    ///@{

    virtual void UpdateBeforeSolveSolutionStep()
    {
        KRATOS_ERROR << "Calling the base class "
                        "ScalarCoSolvingProcess::UpdateBeforeSolveSolutionStep."
                        " Please override it in derrived class.";
    }

    virtual void UpdateAfterSolveSolutionStep()
    {
        KRATOS_ERROR << "Calling the base class "
                        "ScalarCoSolvingProcess::UpdateAfterSolveSolutionStep. "
                        "Please override it in derrived class.";
    }

    virtual void UpdateConvergenceVariable()
    {
        KRATOS_ERROR << "Calling the base class "
                        "ScalarCoSolvingProcess::UpdateConvergenceVariable. "
                        "Please override it in derrived class.";
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::vector<SolvingStrategyType*> mrSolvingStrategiesList;
    Variable<double>& mrConvergenceVariable;

    int mMaxIterations;
    double mConvergenceAbsoluteTolerance;
    double mConvergenceRelativeTolerance;

    ///@}
    ///@name Operations
    ///@{

    void SolveSolutionStep()
    {
        this->UpdateBeforeSolveSolutionStep();
        this->UpdateConvergenceVariable();

        for (SolvingStrategyType* p_solving_strategy : this->mrSolvingStrategiesList)
        {
            p_solving_strategy->InitializeSolutionStep();
            p_solving_strategy->Predict();
        }

        bool is_converged = false;
        int iteration = 1;

        RansVariableUtils rans_variable_utils;

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();

        int iteration_format_length =
            static_cast<int>(std::log10(this->mMaxIterations)) + 1;

        while (!is_converged && iteration <= this->mMaxIterations)
        {
            // Copy the existing convergence variable to the old_convergence_variable
            rans_variable_utils.CopyScalarVar(
                this->mrConvergenceVariable, OLD_CONVERGENCE_VARIABLE, r_nodes);

            for (SolvingStrategyType* p_solving_strategy : this->mrSolvingStrategiesList)
                p_solving_strategy->SolveSolutionStep();

            this->UpdateConvergenceVariable();

            const double increase_norm = rans_variable_utils.GetScalarVariableDifferenceNormSquare(
                r_nodes, OLD_CONVERGENCE_VARIABLE, this->mrConvergenceVariable);
            double solution_norm = rans_variable_utils.GetScalarVariableSolutionNormSquare(
                r_nodes, this->mrConvergenceVariable);

            if (solution_norm <= std::numeric_limits<double>::epsilon())
                solution_norm = 1.0;

            double convergence_relative = increase_norm / solution_norm;
            double convergence_absolute = std::sqrt(increase_norm) / r_nodes.size();

            is_converged = (convergence_relative < this->mConvergenceRelativeTolerance ||
                            convergence_absolute < this->mConvergenceAbsoluteTolerance);

            if (this->mEchoLevel > 1)
            {
                std::stringstream conv_check_msg;
                conv_check_msg
                    << "[Itr.#" << std::setw(iteration_format_length) << iteration
                    << "] CONVERGENCE CHECK: " << mrConvergenceVariable.Name()
                    << " ratio = " << std::setprecision(3) << std::scientific << convergence_relative
                    << "; exp. ratio = " << this->mConvergenceRelativeTolerance
                    << "; abs = " << convergence_absolute
                    << "; exp.abs = " << this->mConvergenceAbsoluteTolerance << "\n";
                KRATOS_INFO(this->Info()) << conv_check_msg.str();

                if (is_converged)
                {
                    std::stringstream conv_msg;
                    conv_msg
                        << "[Itr.#" << std::setw(iteration_format_length) << iteration
                        << "] CONVERGENCE CHECK: " << mrConvergenceVariable.Name()
                        << " *** CONVERGENCE IS ACHIEVED ***\n";
                    KRATOS_INFO(this->Info()) << conv_msg.str();
                }
            }

            iteration++;
        }

        this->UpdateAfterSolveSolutionStep();

        KRATOS_WARNING_IF(this->Info(), !is_converged)
            << "\n-------------------------------------------------------"
            << "\n    WARNING: Max coupling iterations reached.          "
            << "\n             Please increase coupling max_iterations   "
            << "\n             or decrease coupling                      "
            << "\n             relative_tolerance/absolute tolerance     "
            << "\n-------------------------------------------------------"
            << "\n";

        for (SolvingStrategyType* p_solving_strategy : this->mrSolvingStrategiesList)
            p_solving_strategy->FinalizeSolutionStep();
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(
        ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    // ScalarCoSolvingProcess(ScalarCoSolvingProcess const& rOther);

    ///@}

}; // Class ScalarCoSolvingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::istream& operator>>(std::istream& rIStream,
                                ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::ostream& operator<<(std::ostream& rOStream,
                                const ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // SCALAR_CO_SOLVING_PROCESS_H_INCLUDED defined
