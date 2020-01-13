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
#include <limits>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "input_output/vtk_output.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "custom_utilities/rans_variable_utilities.h"

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

    using SolvingStrategyType = SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Pointer definition of ScalarCoSolvingProcess
    KRATOS_CLASS_POINTER_DEFINITION(ScalarCoSolvingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ScalarCoSolvingProcess(ModelPart& rModelPart,
                           Parameters rParameters,
                           Variable<double>& rConvergenceVariable)
        : Process(), mrModelPart(rModelPart), mrConvergenceVariable(rConvergenceVariable)
    {
        Parameters default_parameters = Parameters(R"(
        {
            "relative_tolerance"                : 1e-3,
            "absolute_tolerance"                : 1e-5,
            "max_iterations"                    : 10,
            "echo_level"                        : 0,
            "relaxation_factor"                 : 1.0,
            "number_of_parent_solve_iterations" : 0,
            "vtk_output_settings"               : {},
            "vtk_output_frequency"              : 1,
            "vtk_output_prefix"                 : ""
        })");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = rParameters["echo_level"].GetInt();
        mConvergenceRelativeTolerance = rParameters["relative_tolerance"].GetDouble();
        mConvergenceAbsoluteTolerance = rParameters["absolute_tolerance"].GetDouble();
        mRelaxationFactor = rParameters["relaxation_factor"].GetDouble();
        mMaxIterations = rParameters["max_iterations"].GetInt();
        mSkipIterations = rParameters["number_of_parent_solve_iterations"].GetInt();
        mVtkOutputFrequency = rParameters["vtk_output_frequency"].GetInt();
        mVtkOutputPrefix = rParameters["vtk_output_prefix"].GetString();

        Parameters empty_parameters(R"({})");
        if (!rParameters["vtk_output_settings"].IsEquivalentTo(empty_parameters))
        {
            KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
                << "Adding VtkOutput.\n";
            mpVtkOutput = Kratos::make_unique<VtkOutput>(rModelPart, rParameters["vtk_output_settings"]);
        }

        mCurrentParentIteration = 0;
    }

    /// Destructor.
    ~ScalarCoSolvingProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddStrategy(typename SolvingStrategyType::Pointer pStrategy,
                     const Variable<double>& rScalarVariable)
    {
        mSolvingStrategiesList.push_back(pStrategy);
        mSolvingVariableNamesList.push_back(rScalarVariable.Name());
    }

    void AddAuxiliaryProcess(Process::Pointer pAuxiliaryProcess)
    {
        mAuxiliaryProcessList.push_back(pAuxiliaryProcess);
    }

    void SetParentSolvingStrategy(typename SolvingStrategyType::Pointer pParentSolvingStrategy)
    {
        mpParentSolvingStrategy = pParentSolvingStrategy;
    }

    /// Execute method is used to execute the ScalarCoSolvingProcess algorithms.
    void Execute() override
    {
        if (mrModelPart.GetProcessInfo()[IS_CO_SOLVING_PROCESS_ACTIVE])
            SolveEquations();
    }

    virtual int Check() override
    {
        KRATOS_CHECK_IS_FALSE(!mrModelPart.HasNodalSolutionStepVariable(this->mrConvergenceVariable));

        for (auto strategy : mSolvingStrategiesList)
            strategy->Check();

        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->Check();

        KRATOS_ERROR_IF(mSolvingStrategiesList.size() == 0)
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

    virtual void UpdateBeforeSolveEquations()
    {
        KRATOS_ERROR << "Calling the base class "
                        "ScalarCoSolvingProcess::UpdateBeforeSolveEquations."
                        " Please override it in derrived class.";
    }

    virtual void UpdateAfterSolveEquations()
    {
        KRATOS_ERROR << "Calling the base class "
                        "ScalarCoSolvingProcess::UpdateAfterSolveEquations. "
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

    void SolveEquations()
    {
        ++mCurrentParentIteration;

        if (mCurrentParentIteration > mSkipIterations ||
            mpParentSolvingStrategy->IsConverged())
        {
            mCurrentParentIteration = 0;
            this->UpdateBeforeSolveEquations();

            for (auto p_solving_strategy : this->mSolvingStrategiesList)
            {
                p_solving_strategy->InitializeSolutionStep();
                p_solving_strategy->Predict();
            }

            bool is_converged = false;
            int iteration = 1;

            Communicator& r_communicator = mrModelPart.GetCommunicator();

            ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();
            ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();
            const int parent_solve_iteration = r_current_process_info[NL_ITERATION_NUMBER];

            Vector old_values(r_nodes.size());
            Vector new_values(r_nodes.size());
            Vector delta_values(r_nodes.size());

            int iteration_format_length =
                static_cast<int>(std::log10(this->mMaxIterations)) + 1;

            while (!is_converged && iteration <= this->mMaxIterations)
            {
                RansVariableUtilities::GetNodalVariablesVector(
                    old_values, r_nodes, this->mrConvergenceVariable);

                for (int i = 0;
                     i < static_cast<int>(this->mSolvingStrategiesList.size()); ++i)
                {
                    auto p_solving_strategy = this->mSolvingStrategiesList[i];
                    auto scalar_variable_name = this->mSolvingVariableNamesList[i];

                    p_solving_strategy->SolveSolutionStep();
                    const unsigned int iterations =
                        r_current_process_info[NL_ITERATION_NUMBER];
                    KRATOS_INFO_IF(this->Info(), this->mEchoLevel > 0)
                        << "Solving " << scalar_variable_name << " used "
                        << iterations << " iterations.\n";
                }

                ++mVtkOutputStep;
                if (mpVtkOutput && mVtkOutputStep >= mVtkOutputFrequency)
                {
                    mVtkOutputStep = 0;
                    std::stringstream s_label;
                    s_label << mVtkOutputPrefix << "_Rank_"
                            << mrModelPart.GetCommunicator().MyPID() << std::fixed
                            << "_Step_" << r_current_process_info[STEP]
                            << "_ParentItr_" << parent_solve_iteration
                            << std::setw(iteration_format_length) << "_CoupleItr_"
                            << std::setfill('0') << iteration;
                    mpVtkOutput->PrintOutput(s_label.str());
                    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1) << "Writing Vtk output to " << s_label.str() << "\n";
                }

                this->UpdateConvergenceVariable();

                RansVariableUtilities::GetNodalVariablesVector(
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
                const std::vector<double>& total_residual_norms =
                    r_communicator.GetDataCommunicator().SumAll(residual_norms);

                noalias(new_values) = old_values + delta_values * mRelaxationFactor;
                RansVariableUtilities::SetNodalVariables(
                    r_nodes, new_values, this->mrConvergenceVariable);
                r_communicator.SynchronizeVariable(this->mrConvergenceVariable);

                double convergence_relative =
                    total_residual_norms[0] /
                    (total_residual_norms[1] <= std::numeric_limits<double>::epsilon()
                         ? 1.0
                         : total_residual_norms[1]);
                double convergence_absolute =
                    std::sqrt(total_residual_norms[0]) / total_residual_norms[2];

                is_converged =
                    (convergence_relative < this->mConvergenceRelativeTolerance ||
                     convergence_absolute < this->mConvergenceAbsoluteTolerance);

                if (this->mEchoLevel > 1)
                {
                    std::stringstream conv_check_msg;
                    conv_check_msg
                        << "[Itr.#" << std::setw(iteration_format_length)
                        << iteration << "/" << this->mMaxIterations
                        << "] CONVERGENCE CHECK: " << mrConvergenceVariable.Name()
                        << " ratio = " << std::setprecision(3)
                        << std::scientific << convergence_relative
                        << "; exp. ratio = " << this->mConvergenceRelativeTolerance
                        << "; abs = " << convergence_absolute
                        << "; exp.abs = " << this->mConvergenceAbsoluteTolerance << "\n";
                    KRATOS_INFO(this->Info()) << conv_check_msg.str();

                    if (is_converged)
                    {
                        std::stringstream conv_msg;
                        conv_msg
                            << "[Itr.#" << std::setw(iteration_format_length)
                            << iteration << "/" << this->mMaxIterations
                            << "] CONVERGENCE CHECK: " << mrConvergenceVariable.Name()
                            << " *** CONVERGENCE IS ACHIEVED ***\n";
                        KRATOS_INFO(this->Info()) << conv_msg.str();
                    }
                }

                iteration++;
            }

            this->UpdateAfterSolveEquations();

            KRATOS_INFO_IF(this->Info(), !is_converged && this->mEchoLevel > 2)
                << "\n-------------------------------------------------------"
                << "\n    INFO: Max coupling iterations reached.             "
                << "\n          Please increase coupling max_iterations      "
                << "\n          or decrease coupling                         "
                << "\n          relative_tolerance/absolute tolerance        "
                << "\n-------------------------------------------------------"
                << "\n";

            for (auto p_solving_strategy : this->mSolvingStrategiesList)
                p_solving_strategy->FinalizeSolutionStep();
        }
        else
        {
            KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
                << "Skipping co-solving process for parent solve to continue, "
                   "since parent solve itertions are less than "
                   "\"number_of_parent_solve_iterations\" [ "
                << mCurrentParentIteration << " <= " << mSkipIterations << " ].\n";
        }
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::vector<typename SolvingStrategyType::Pointer> mSolvingStrategiesList;
    std::vector<std::string> mSolvingVariableNamesList;
    std::vector<Process::Pointer> mAuxiliaryProcessList;
    std::string mVtkOutputPrefix;
    Variable<double>& mrConvergenceVariable;

    VtkOutput::UniquePointer mpVtkOutput;

    typename SolvingStrategyType::Pointer mpParentSolvingStrategy;

    int mMaxIterations;
    int mSkipIterations;
    int mCurrentParentIteration;
    int mVtkOutputStep = 0;
    int mVtkOutputFrequency;
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
