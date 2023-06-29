//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "includes/fill_communicator.h"
#include "includes/parallel_environment.h"

namespace Kratos
{

FillCommunicator::FillCommunicator(ModelPart& rModelPart)
    : FillCommunicator(rModelPart, ParallelEnvironment::GetDefaultDataCommunicator())
{}

FillCommunicator::FillCommunicator(
    ModelPart& rModelPart,
    const DataCommunicator& rDataCommunicator)
    : mrDataComm(rDataCommunicator), mrBaseModelPart(rModelPart)
{}

void FillCommunicator::Execute()
{
    KRATOS_TRY
    KRATOS_CATCH("");
}

void FillCommunicator::PrintDebugInfo()
{
    PrintModelPartDebugInfo(mrBaseModelPart);
}

void FillCommunicator::PrintModelPartDebugInfo(const ModelPart& rModelPart)
{
    KRATOS_TRY

    std::cout.flush();
    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_data_communicator = r_communicator.GetDataCommunicator();
    r_data_communicator.Barrier();

    // Check rank and number of processors
    const int rank = r_data_communicator.Rank();
    const int num_processors = r_data_communicator.Size();
    KRATOS_ERROR_IF_NOT(rank == 0) << "Serial FillCommunicator current rank is not 0." << std::endl;
    KRATOS_ERROR_IF_NOT(num_processors == 1) << "Serial FillCommunicator number of processors larger than 1." << std::endl;

    // Check local and ghost mesh
    KRATOS_ERROR_IF_NOT(r_communicator.NeighbourIndices().size() == 0) << "There are not expected neighbour indices" << std::endl;
    KRATOS_ERROR_IF_NOT(r_communicator.GhostMesh().NumberOfNodes() == 0) << "There are unexpected nodes in the ghost mesh" << std::endl;
    KRATOS_ERROR_IF_NOT(r_communicator.InterfaceMesh().NumberOfNodes() == 0) << "There are unexpected nodes in the interface mesh." << std::endl;

    KRATOS_CATCH("");
}

std::string FillCommunicator::Info() const
{
    std::stringstream buffer;
    buffer << "FillCommunicator";
    return buffer.str();
}

void FillCommunicator::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FillCommunicator" << std::endl;
}

void FillCommunicator::PrintData(std::ostream& rOStream) const
{
}

}
