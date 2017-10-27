#include "hdf5_model_part_io.h"
#include <vector>

namespace Kratos
{

HDF5ModelPartIO::HDF5ModelPartIO(Parameters& rParams, HDF5File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;
    Check();
    KRATOS_CATCH("");
}

bool HDF5ModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    HDF5File::Vector<int> node_ids;
    HDF5File::Vector<array_1d<double, 3>> node_coords;
    rNodes.clear();

    const std::size_t num_nodes = ReadNodesNumber();
    GetFile().ReadDataSet("/Nodes/Local/Id", node_ids, 0, num_nodes);
    GetFile().ReadDataSet("/Nodes/Local/Coordinate", node_coords, 0, num_nodes);

    rNodes.reserve(num_nodes);
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        const array_1d<double, 3>& r_coord = node_coords[i];
        NodeType::Pointer p_node = boost::make_shared<NodeType>(
            node_ids[i], r_coord[0], r_coord[1], r_coord[2]);
        rNodes.push_back(p_node);
    }

    return true;
    KRATOS_CATCH("");
}

std::size_t HDF5ModelPartIO::ReadNodesNumber()
{
    const std::vector<unsigned> dims = GetFile().GetDataDimensions("/Nodes/Local/Id");
    return dims[0];
}

void HDF5ModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    HDF5File::Vector<int> node_ids(rNodes.size());
    HDF5File::Vector<array_1d<double, 3>> node_coords(rNodes.size());

    unsigned pos = 0;
    for (const auto& r_node : rNodes)
    {
        node_ids[pos] = r_node.Id();
        node_coords[pos] = r_node.Coordinates();
        ++pos;
    }

    GetFile().WriteDataSet("/Nodes/Local/Id", node_ids);
    GetFile().WriteDataSet("/Nodes/Local/Coordinate", node_coords);
    
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
    KRATOS_TRY;

    HDF5File::Vector<int> elem_ids(rElements.size());
    HDF5File::Vector<int> elem_property_ids(rElements.size());
    //HDF5File::Vector<int> elem_topology;
    //elem_topology.reserve(rElements.size() * (rElements.front().GetGeometry().size() + 2));

    unsigned pos = 0;
    for (const auto& r_elem : rElements)
    {
        elem_ids[pos] = r_elem.Id();
        elem_property_ids[pos] = r_elem.GetProperties().Id();
        ++pos;
    }

    GetFile().WriteDataSet("/Elements/Id", elem_ids);
    GetFile().WriteDataSet("/Elements/PropertyId", elem_property_ids);

    KRATOS_CATCH("");
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

HDF5File& HDF5ModelPartIO::GetFile() const
{
    return *mpFile;
}

void HDF5ModelPartIO::Check()
{
    KRATOS_ERROR_IF(mpFile->GetTotalProcesses() != 1)
        << "Serial IO expects file access by a single process only." << std::endl;
}

} // namespace Kratos.
