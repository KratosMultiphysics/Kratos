#include "hdf5_model_part_io.h"

namespace Kratos
{

HDF5ModelPartIO::HDF5ModelPartIO(Parameters& rParams, HDF5File::Pointer pFile)
: mpFile(pFile)
{
    m_pid = pFile->GetPID();
    m_total_processes = pFile->GetTotalProcesses();
}

bool HDF5ModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    if (GetTotalProcesses() == 1)
        ReadNodesSerial(rNodes);
    else
        ReadNodesParallel(rNodes);
    
    return true;
    KRATOS_CATCH("");
}

std::size_t HDF5ModelPartIO::ReadNodesNumber()
{
    return 0;
}

void HDF5ModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    if (GetTotalProcesses() == 1)
        WriteNodesSerial(rNodes);
    else
        WriteNodesParallel(rNodes);

    KRATOS_CATCH("");
}

void HDF5ModelPartIO::ReadElements(NodesContainerType& rNodes,
                                   PropertiesContainerType& rProperties,
                                   ElementsContainerType& rElements)
{
}

std::size_t HDF5ModelPartIO::ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
{
    return 0;
}

void HDF5ModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
}

void HDF5ModelPartIO::ReadConditions(NodesContainerType& rNodes,
                                     PropertiesContainerType& rProperties,
                                     ConditionsContainerType& rConditions)
{
}

std::size_t HDF5ModelPartIO::ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
{
    return 0;
}

void HDF5ModelPartIO::ReadInitialValues(ModelPart& rModelPart)
{
}

void HDF5ModelPartIO::ReadInitialValues(NodesContainerType& rNodes,
                                        ElementsContainerType& rElements,
                                        ConditionsContainerType& rConditions)
{
}

void HDF5ModelPartIO::ReadModelPart(ModelPart& rModelPart)
{
}

void HDF5ModelPartIO::WriteModelPart(ModelPart& rModelPart)
{
}

unsigned HDF5ModelPartIO::GetPID() const
{
    return m_pid;
}

unsigned HDF5ModelPartIO::GetTotalProcesses() const
{
    return m_total_processes;
}

void HDF5ModelPartIO::WriteNodesSerial(NodesContainerType const& rNodes) const
{

}

void HDF5ModelPartIO::WriteNodesParallel(NodesContainerType const& rNodes) const
{

}

void HDF5ModelPartIO::ReadNodesSerial(NodesContainerType const& rNodes) const
{

}

void HDF5ModelPartIO::ReadNodesParallel(NodesContainerType const& rNodes) const
{

}
} // namespace Kratos.
