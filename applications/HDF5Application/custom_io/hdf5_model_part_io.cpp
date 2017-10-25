#include "hdf5_model_part_io.h"

namespace Kratos
{

HDF5ModelPartIO::HDF5ModelPartIO(Parameters& rParams, HDF5File::Pointer pFile)
: mpFile(pFile)
{
    m_pid = pFile->GetPID();
}

bool HDF5ModelPartIO::ReadNodes(NodesContainerType& rThisNodes)
{
    return true;
}

std::size_t HDF5ModelPartIO::ReadNodesNumber()
{
    return 0;
}

void HDF5ModelPartIO::WriteNodes(NodesContainerType const& rThisNodes)
{

}

void HDF5ModelPartIO::ReadElements(NodesContainerType& rThisNodes,
                                   PropertiesContainerType& rThisProperties,
                                   ElementsContainerType& rThisElements)
{
}

std::size_t HDF5ModelPartIO::ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
{
    return 0;
}

void HDF5ModelPartIO::WriteElements(ElementsContainerType const& rThisElements)
{
}

void HDF5ModelPartIO::ReadConditions(NodesContainerType& rThisNodes,
                                     PropertiesContainerType& rThisProperties,
                                     ConditionsContainerType& rThisConditions)
{
}

std::size_t HDF5ModelPartIO::ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
{
    return 0;
}

void HDF5ModelPartIO::ReadInitialValues(ModelPart& rThisModelPart)
{
}

void HDF5ModelPartIO::ReadInitialValues(NodesContainerType& rThisNodes,
                                        ElementsContainerType& rThisElements,
                                        ConditionsContainerType& rThisConditions)
{
}

void HDF5ModelPartIO::ReadModelPart(ModelPart& rThisModelPart)
{
}

void HDF5ModelPartIO::WriteModelPart(ModelPart& rThisModelPart)
{
}

unsigned HDF5ModelPartIO::GetPID() const
{
    return m_pid;
}
} // namespace Kratos.
