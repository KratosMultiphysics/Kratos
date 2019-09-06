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
            "echo_level"            : 0,
            "relaxation_factor"     : 1.0
        })");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = rParameters["echo_level"].GetInt();
        mConvergenceRelativeTolerance = rParameters["relative_tolerance"].GetDouble();
        mConvergenceAbsoluteTolerance = rParameters["absolute_tolerance"].GetDouble();
        mRelaxationFactor = rParameters["relaxation_factor"].GetDouble();
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

    void AddStrategy(SolvingStrategyType* pStrategy, Variable<double>* pScalarVariable)
    {
        mrSolvingStrategiesList.push_back(pStrategy);
        mrSolvingVariablesList.push_back(pScalarVariable);
    }

    void AddAuxiliaryProcess(Process::Pointer pAuxiliaryProcess)
    {
        mAuxiliaryProcessList.push_back(pAuxiliaryProcess);
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
        KRATOS_CHECK_VARIABLE_KEY(this->mrConvergenceVariable);

        KRATOS_CHECK_IS_FALSE(!mrModelPart.HasNodalSolutionStepVariable(this->mrConvergenceVariable));

        for (SolvingStrategyType* strategy : mrSolvingStrategiesList)
            strategy->Check();

        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->Check();

        KRATOS_ERROR_IF(mrSolvingStrategiesList.size() == 0)
            << "No strategies are found for ScalarCoSolvingProcess.";

        return 0;
    }

    virtual void ExecuteInitialize() override
    {
        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->ExecuteInitialize();
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

    void ExecuteAuxiliaryProcesses()
    {
        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->Execute();
    }

    void ExecuteAuxiliaryProcessesInitializeSolutionStep()
    {
        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->ExecuteInitializeSolutionStep();
    }

    void SolveSolutionStep()
    {
        this->UpdateBeforeSolveSolutionStep();

        for (SolvingStrategyType* p_solving_strategy : this->mrSolvingStrategiesList)
        {
            p_solving_strategy->InitializeSolutionStep();
            p_solving_strategy->Predict();
        }

        bool is_converged = false;
        int iteration = 1;

        RansVariableUtils rans_variable_utils;

        Communicator& r_communicator = mrModelPart.GetCommunicator();

        ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();
        const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

        Vector old_values(r_nodes.size());
        Vector new_values(r_nodes.size());
        Vector delta_values(r_nodes.size());

        int iteration_format_length =
            static_cast<int>(std::log10(this->mMaxIterations)) + 1;

        while (!is_converged && iteration <= this->mMaxIterations)
        {
            rans_variable_utils.GetNodalVariablesVector(
                old_values, r_nodes, this->mrConvergenceVariable);

            for (int i = 0;
                 i < static_cast<int>(this->mrSolvingStrategiesList.size()); ++i)
            {
                auto p_solving_strategy = this->mrSolvingStrategiesList[i];
                auto p_scalar_variable = this->mrSolvingVariablesList[i];

                p_solving_strategy->SolveSolutionStep();
                const unsigned int iterations = r_current_process_info[NL_ITERATION_NUMBER];
                KRATOS_INFO_IF(this->Info(), this->mEchoLevel > 0)
                    << "Solving " << p_scalar_variable->Name() << " used "
                    << iterations << " iterations.\n";
            }

            this->UpdateConvergenceVariable();

            rans_variable_utils.GetNodalVariablesVector(
                new_values, r_nodes, this->mrConvergenceVariable);
            noalias(delta_values) = new_values - old_values;

            // This vector stores norms of the residual
            // index - 0 : increase_norm
            // index - 1 : solution_norm
            // index - 3 : number of nodes
            std::vector<double> residual_norms(3);
            residual_norms[0] = std::pow(norm_2(delta_values), 2);
            residual_norms[1] = std::pow(norm_2(new_values), 2);
            residual_norms[2] = static_cast<double>(r_nodes.size());
            r_communicator.GetDataCommunicator().SumAll(residual_norms);

            noalias(new_values) = old_values + delta_values * mRelaxationFactor;
            rans_variable_utils.SetNodalVariables(new_values, r_nodes,
                                                  this->mrConvergenceVariable);
            r_communicator.SynchronizeVariable(this->mrConvergenceVariable);

            if (residual_norms[1] <= std::numeric_limits<double>::epsilon())
                residual_norms[1] = 1.0;

            double convergence_relative = residual_norms[0] / residual_norms[1];
            double convergence_absolute =
                std::sqrt(residual_norms[0]) / residual_norms[2];

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
                    conv_msg << "[Itr.#" << std::setw(iteration_format_length) << iteration
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::vector<SolvingStrategyType*> mrSolvingStrategiesList;
    std::vector<Variable<double>*> mrSolvingVariablesList;
    std::vector<Process::Pointer> mAuxiliaryProcessList;
    Variable<double>& mrConvergenceVariable;

    int mMaxIterations;
    double mConvergenceAbsoluteTolerance;
    double mConvergenceRelativeTolerance;
    double mRelaxationFactor;

    ///@}
    ///@name Operations
    ///@{

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
