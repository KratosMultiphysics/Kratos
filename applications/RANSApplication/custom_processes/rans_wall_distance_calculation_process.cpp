//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "factories/linear_solver_factory.h"
#include "includes/communicator.h"
#include "includes/variables.h"
#include "linear_solvers/linear_solver.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/variational_distance_calculation_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "rans_wall_distance_calculation_process.h"

namespace Kratos
{
RansWallDistanceCalculationProcess::RansWallDistanceCalculationProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mMaxIterations = rParameters["max_iterations"].GetInt();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mWallFlagVariableName = rParameters["wall_flag_variable_name"].GetString();
    mWallFlagVariableValue = rParameters["wall_flag_variable_value"].GetBool();
    mRecalculateAtEachTimeStep = rParameters["re_calculate_at_each_time_step"].GetBool();
    mCorrectDistancesUsingNeighbors =
        rParameters["correct_distances_using_neighbors"].GetBool();
    mLinearSolverParameters = rParameters["linear_solver_settings"];

    KRATOS_CATCH("");
}

Process::Pointer RansWallDistanceCalculationProcess::GetWallDistanceCalculationProcess(
    ModelPart& rModelPart,
    Parameters LinearSolverParameters,
    const int MaxIterations)
{
    KRATOS_TRY

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;

    auto p_linear_solver =
        LinearSolverFactory<SparseSpaceType, DenseSpaceType>()
            .Create(LinearSolverParameters);

    auto p_builder_and_solver =
        Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType>>(
            p_linear_solver);

    const int domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

    Process::Pointer p_process = nullptr;

    if (domain_size == 2) {
        p_process =
            Kratos::make_shared<VariationalDistanceCalculationProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType>>(
                rModelPart, p_linear_solver, p_builder_and_solver, MaxIterations);
    } else if (domain_size == 3) {
        p_process =
            Kratos::make_shared<VariationalDistanceCalculationProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType>>(
                rModelPart, p_linear_solver, p_builder_and_solver, MaxIterations);
    } else {
        KRATOS_ERROR << "Unknown domain size = " << domain_size;
    }

    return p_process;

    KRATOS_CATCH("");
}

int RansWallDistanceCalculationProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(DISTANCE))
        << "DISTANCE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0.0;

    KRATOS_CATCH("");
}

void RansWallDistanceCalculationProcess::ExecuteInitialize()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    mpVariationalDistanceCalculationProcess = this->GetWallDistanceCalculationProcess(
        r_model_part, mLinearSolverParameters, mMaxIterations);

    CalculateWallDistances();

    KRATOS_CATCH("");
}

void RansWallDistanceCalculationProcess::ExecuteInitializeSolutionStep()
{
    if (mRecalculateAtEachTimeStep) {
        CalculateWallDistances();
    }
}

std::string RansWallDistanceCalculationProcess::Info() const
{
    return std::string("RansWallDistanceCalculationProcess");
}

void RansWallDistanceCalculationProcess::CalculateWallDistances()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const Flags& r_wall_flag = KratosComponents<Flags>::Get(mWallFlagVariableName);

    VariableUtils variable_utilities;

    variable_utilities.SetVariable(DISTANCE, 1.0, r_model_part.Nodes());
    variable_utilities.SetVariable(DISTANCE, 0.0, r_model_part.Nodes(),
                                   r_wall_flag, mWallFlagVariableValue);

    mpVariationalDistanceCalculationProcess->Execute();

    if (mCorrectDistancesUsingNeighbors) {
        CorrectWallDistances();
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Wall distances calculated in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

void RansWallDistanceCalculationProcess::CorrectWallDistances()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_communicator = r_model_part.GetCommunicator();
    auto& r_data_communicator = r_communicator.GetDataCommunicator();

    auto& r_nodes = r_communicator.LocalMesh().Nodes();

    VariableUtils().SetNonHistoricalVariableToZero(DISTANCE, r_model_part.Nodes());

    FindGlobalNodalNeighboursProcess find_nodal_neighbours_process(
        r_data_communicator, r_model_part);
    find_nodal_neighbours_process.Execute();

    class GlobalPointerAdder
    {
    public:
        typedef GlobalPointersVector<ModelPart::NodeType> value_type;
        value_type gp_vector;

        value_type GetValue()
        {
            gp_vector.Unique();
            return gp_vector;
        }

        void LocalReduce(const value_type& rGPVector)
        {
            for (auto& r_gp : rGPVector.GetContainer()) {
                this->gp_vector.push_back(r_gp);
            }
        }
        void ThreadSafeReduce(GlobalPointerAdder& rOther)
        {
#pragma omp critical
            {
                for (auto& r_gp : rOther.gp_vector.GetContainer()) {
                    this->gp_vector.push_back(r_gp);
                }
            }
        }
    };

    auto all_global_pointers =
        block_for_each<GlobalPointerAdder>(r_nodes, [](ModelPart::NodeType& rNode) {
            return rNode.GetValue(NEIGHBOUR_NODES);
        });

    GlobalPointerCommunicator<ModelPart::NodeType> pointer_comm(
        r_data_communicator, all_global_pointers);

    auto distance_proxy =
        pointer_comm.Apply([](const GlobalPointer<ModelPart::NodeType>& gp) -> double {
            return gp->FastGetSolutionStepValue(DISTANCE);
        });

    int number_of_modified_nodes = block_for_each<SumReduction<int>>(
        r_nodes, [&](ModelPart::NodeType& rNode) -> int {
            if (rNode.FastGetSolutionStepValue(DISTANCE) < 0.0) {
                const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);

                int count = 0;
                double average_value = 0.0;
                for (int j_node = 0;
                     j_node < static_cast<int>(r_neighbours.size()); ++j_node) {
                    const double j_distance = distance_proxy.Get(r_neighbours(j_node));
                    if (j_distance > 0.0) {
                        average_value += j_distance;
                        count++;
                    }
                }

                if (count > 0) {
                    rNode.SetValue(DISTANCE, average_value / static_cast<double>(count));
                } else {
                    KRATOS_ERROR
                        << "Node " << rNode.Id() << " at " << rNode.Coordinates()
                        << " didn't find any neighbour(" << r_neighbours.size()
                        << ") with positive DISTANCE. Please recheck "
                        << mModelPartName << ".\n";
                }

                return 1;
            }
            return 0;
        });

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        double& r_distance = rNode.FastGetSolutionStepValue(DISTANCE);
        const double avg_distance = rNode.GetValue(DISTANCE);
        if (r_distance < 0.0) {
            r_distance = avg_distance;
        }
    });

    r_communicator.SynchronizeVariable(DISTANCE);
    number_of_modified_nodes =
        r_communicator.GetDataCommunicator().SumAll(number_of_modified_nodes);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && number_of_modified_nodes > 0)
        << "Corrected " << number_of_modified_nodes
        << " nodal wall distances in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

const Parameters RansWallDistanceCalculationProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                  : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "max_iterations"                   : 10,
            "echo_level"                       : 0,
            "wall_flag_variable_name"          : "STRUCTURE",
            "wall_flag_variable_value"         : true,
            "re_calculate_at_each_time_step"   : false,
            "correct_distances_using_neighbors": true,
            "linear_solver_settings" : {
                "solver_type"     : "amgcl"
            }
        })");

    return default_parameters;
}

} // namespace Kratos.
