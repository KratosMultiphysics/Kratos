//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

#include <sstream>

#include "mpi.h"

#include "containers/model.h"
#include "includes/data_communicator.h"
#include "includes/model_part.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/debug_utilities.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckSingleHistoricalVariableValue, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const int world_size = r_comm.Size();

    Model model;
    ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Put some nodes in every partition
    for(int i = 0; i < world_size; i++) {
        auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

        node->FastGetSolutionStepValue(PRESSURE) = i%world_size;
        node->FastGetSolutionStepValue(TEMPERATURE) = (i%world_size)*10;
        node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
    }

    // Build the communicator
    ParallelFillCommunicator(model_part, r_comm).Execute();

    MpiDebugUtilities::CheckHistoricalVariable(model_part, PRESSURE);
    MpiDebugUtilities::CheckHistoricalVariable(model_part, TEMPERATURE);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckSingleHistoricalVariableFixity, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const int world_size = r_comm.Size();

    Model model;
    ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Put some nodes in every partition
    for(int i = 0; i < world_size; i++) {
        auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

        node->FastGetSolutionStepValue(PRESSURE) = 0;
        // Not optimal but should scramble the values enough
        if(i%2) {
            node->Fix(PRESSURE);
        }
        node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
    }

    // Build the communicator
    ParallelFillCommunicator(model_part, r_comm).Execute();

    MpiDebugUtilities::CheckHistoricalVariable(model_part, PRESSURE);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckSingleHistoricalVariableValueError, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const int world_rank = r_comm.Rank();
    const int world_size = r_comm.Size();

    Model model;
    ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Put some nodes in every partition
    for(int i = 0; i < world_size; i++) {
        auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

        node->FastGetSolutionStepValue(PRESSURE) = (world_rank == 0);
        node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
    }

    // Build the communicator
    ParallelFillCommunicator(model_part, r_comm).Execute();

    std::stringstream error_message;

    error_message << "Value error(s) found" << std::endl;

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        MpiDebugUtilities::CheckHistoricalVariable(model_part, PRESSURE),
        error_message.str()
    );
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckSingleHistoricalVariableFixityError, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const int world_rank = r_comm.Rank();
    const int world_size = r_comm.Size();

    Model model;
    ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Put some nodes in every partition
    for(int i = 0; i < world_size; i++) {
        auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

        node->FastGetSolutionStepValue(PRESSURE) = 0;
        if(world_rank == 0) {
            node->Fix(PRESSURE);
        }
        node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
    }

    // Build the communicator
    ParallelFillCommunicator(model_part, r_comm).Execute();

    std::stringstream error_message;

    error_message << "Fixity error(s) found" << std::endl;

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        MpiDebugUtilities::CheckHistoricalVariable(model_part, PRESSURE),
        error_message.str()
    );
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckSingleHistoricalVariableCombinedError, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const int world_rank = r_comm.Rank();
    const int world_size = r_comm.Size();

    Model model;
    ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Put some nodes in every partition
    for(int i = 0; i < world_size; i++) {
        auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

        node->FastGetSolutionStepValue(PRESSURE) = (world_rank == 0);
        if(world_rank == 0) {
            node->Fix(PRESSURE);
        }
        node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
    }

    // Build the communicator
    ParallelFillCommunicator(model_part, r_comm).Execute();

    std::stringstream error_message;

    error_message << "Value error(s) found" << std::endl;
    error_message << "Fixity error(s) found" << std::endl;

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        MpiDebugUtilities::CheckHistoricalVariable(model_part, PRESSURE),
        error_message.str()
    );
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckSingleNonHistoricalVariableValue, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const int world_size = r_comm.Size();

    Model model;
    ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Put some nodes in every partition
    for(int i = 0; i < world_size; i++) {
        auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

        node->SetValue(PRESSURE, i%world_size);
        node->SetValue(TEMPERATURE, (i%world_size)*10);
        node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
    }

    // Build the communicator
    ParallelFillCommunicator(model_part, r_comm).Execute();

    MpiDebugUtilities::CheckNonHistoricalVariable(model_part, model_part.Nodes(), PRESSURE);
    MpiDebugUtilities::CheckNonHistoricalVariable(model_part, model_part.Nodes(), TEMPERATURE);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckSingleNonHistoricalVariableValueError, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const int world_rank = r_comm.Rank();
    const int world_size = r_comm.Size();

    Model model;
    ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Put some nodes in every partition
    for(int i = 0; i < world_size; i++) {
        auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

        node->SetValue(PRESSURE, (world_rank == 0));
        node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
    }

    // Build the communicator
    ParallelFillCommunicator(model_part, r_comm).Execute();

    std::stringstream error_message;

    error_message << "Value error(s) found" << std::endl;

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        MpiDebugUtilities::CheckNonHistoricalVariable(model_part, model_part.Nodes(), PRESSURE),
        error_message.str()
    );
}

// This will work with #5091 or when we move to C++17
// KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DebugToolsCheckMultipleHistoricalVariablesValue, KratosMPICoreFastSuite)
// {
//     const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
//     const int world_rank = r_comm.Rank();
//     const int world_size = r_comm.Size();

//     Model model;
//     ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

//     model_part.AddNodalSolutionStepVariable(PRESSURE);
//     model_part.AddNodalSolutionStepVariable(TEMPERATURE);
//     model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

//     // Put some nodes in every partition
//     for(int i = 0; i < world_size; i++) {
//         auto node = model_part.CreateNewNode(i, 0.0, 0.0, 0.0);

//         node->FastGetSolutionStepValue(PRESSURE) = i%world_size;
//         node->FastGetSolutionStepValue(TEMPERATURE) = (i%world_size)*10;
//         node->FastGetSolutionStepValue(PARTITION_INDEX) = i%world_size;
//     }

//     // Build the communicator
//     ParallelFillCommunicator(model_part, r_comm).Execute();

//     MpiDebugUtilities::CheckNodalHistoricalDatabase(model_part);
// }

}
}