#include "hdf5_model_part_io.h"
#include <vector>
#include <string>
#include <sstream>

namespace Kratos
{

HDF5ModelPartIO::HDF5ModelPartIO(Parameters& rParams, HDF5File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "list_of_elements": [],
            "list_of_conditions": []
        })");

    rParams.ValidateAndAssignDefaults(default_params);

    mElementNames.resize(rParams["list_of_elements"].size());
    for (unsigned i = 0; i < mElementNames.size(); ++i)
        mElementNames[i] = rParams["list_of_elements"].GetArrayItem(i).GetString();

    mConditionNames.resize(rParams["list_of_conditions"].size());
    for (unsigned i = 0; i < mConditionNames.size(); ++i)
        mConditionNames[i] = rParams["list_of_conditions"].GetArrayItem(i).GetString();

    Check();

    for (const std::string& r_elem_name : mElementNames)
    {
        const Element& r_elem = KratosComponents<Element>::Get(r_elem_name);
        std::stringstream elem_info;
        r_elem.PrintInfo(elem_info);
        mElementNamesToPrintInfos[r_elem_name] = elem_info.str();
        mElementPrintInfosToNames[elem_info.str()] = r_elem_name;
    }

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

    // Write the headers.
    for (const std::string& r_elem_name : mElementNames)
        GetFile().AddPath("/Elements/" + r_elem_name);

    if (mElementNames.size() == 1)
        WriteUniformElements(rElements);
    else
        WriteMixedElements(rElements);

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
    KRATOS_TRY;

    KRATOS_ERROR_IF(mpFile->GetTotalProcesses() != 1)
        << "Serial IO expects file access by a single process only." << std::endl;

    for (const std::string& r_elem_name : mElementNames)
    {
        KRATOS_ERROR_IF_NOT(KratosComponents<Element>::Has(r_elem_name))
            << "Element " << r_elem_name << " not registered." << std::endl;
    }

    for (const std::string& r_cond_name : mConditionNames)
    {
        KRATOS_ERROR_IF_NOT(KratosComponents<Condition>::Has(r_cond_name))
            << "Condition " << r_cond_name << " not registered." << std::endl;
    }

    KRATOS_CATCH("");
}

void HDF5ModelPartIO::WriteUniformElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    unsigned size1 = rElements.size();
    unsigned size2 = rElements.front().GetGeometry().size();
    HDF5File::Vector<int> elem_ids(size1);
    HDF5File::Vector<int> prop_ids(size1);
    HDF5File::Matrix<int> connectivities(size1, size2);

    const std::string& r_expected_print_info =
        mElementNamesToPrintInfos.at(mElementNames[0]);
    unsigned pos = 0;
    // Fill arrays and perform checks.
    for (const auto& r_elem : rElements)
    {
        // Check for matching print info string.
        std::stringstream elem_info;
        r_elem.PrintInfo(elem_info);
        const std::string& r_elem_info = elem_info.str();
        KRATOS_ERROR_IF(r_elem_info != r_expected_print_info)
            << "Element " << r_elem.Id() << " with print info string: \n"
            << r_elem_info << "\n is not " << mElementNames[0]
            << " with print info string: \n"
            << r_expected_print_info << std::endl;
        // Fill ids.
        elem_ids[pos] = r_elem.Id();
        prop_ids[pos] = r_elem.GetProperties().Id();
        // Check for matching node count.
        const Element::GeometryType& r_geom = r_elem.GetGeometry();
        KRATOS_ERROR_IF(r_geom.size() != size2)
            << "Geometry has non-uniform size." << std::endl;
        // Fill connectivities.
        for (unsigned i = 0; i < size2; ++i)
            connectivities(pos, i) = r_geom[i].Id();
        ++pos;
    }

    // Write to HDF5 file.
    GetFile().WriteDataSet("/Elements/" + mElementNames[0] + "/Id", elem_ids);
    GetFile().WriteDataSet("/Elements/" + mElementNames[0] + "/PropertyId", prop_ids);
    GetFile().WriteDataSet("/Elements/" + mElementNames[0] + "/Connectivity", connectivities);

    KRATOS_CATCH("");
}

void HDF5ModelPartIO::WriteMixedElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    std::unordered_map<std::string, unsigned> num_elems;
    std::unordered_map<std::string, unsigned> num_nodes;
    for (const std::string& r_elem_name : mElementNames)
    {
        num_elems[r_elem_name] = 0; // Initialize count.
        num_nodes[r_elem_name] = KratosComponents<Element>::Get(r_elem_name).GetGeometry().size();
    }

    // Count number of each element and perform checks.
    for (const auto& r_elem : rElements)
    {
        std::stringstream elem_info;
        r_elem.PrintInfo(elem_info);
        const std::string& r_elem_info = elem_info.str();
        const std::string& r_elem_name = mElementPrintInfosToNames.at(r_elem_info); // Throw here if print info not found in registry.
        // Check for matching node count.
        KRATOS_ERROR_IF(r_elem.GetGeometry().size() != num_nodes.at(r_elem_name))
            << "Non-uniform connectivity detected for element #"
            << r_elem.Id() << std::endl;
        // Increment element count.
        num_elems.at(r_elem_name) += 1;
    }
    
    // Initialized arrays.
    std::unordered_map<std::string, unsigned> pos;
    std::unordered_map<std::string, HDF5File::Vector<int>> elem_ids;
    std::unordered_map<std::string, HDF5File::Vector<int>> prop_ids;
    std::unordered_map<std::string, HDF5File::Matrix<int>> connectivities;
    for (const std::string& r_elem_name : mElementNames)
    {
        pos[r_elem_name] = 0; // Initialize position indices.
        elem_ids[r_elem_name].resize(num_elems[r_elem_name], false);
        prop_ids[r_elem_name].resize(num_elems[r_elem_name], false);
        connectivities[r_elem_name].resize(num_elems[r_elem_name], num_nodes[r_elem_name], false);
    }

    // Fill arrays. By now all element checks have been performed so no exceptions should occur.
    for (const auto& r_elem : rElements)
    {
        // Identity the element name from print info.
        std::stringstream elem_info;
        r_elem.PrintInfo(elem_info);
        const std::string& r_elem_info = elem_info.str();
        const std::string& r_elem_name = mElementPrintInfosToNames.at(r_elem_info);
        // Get the current position for this element.
        unsigned& r_pos = pos[r_elem_name];
        // Fill the element and property ids.
        elem_ids[r_elem_name][r_pos] = r_elem.Id();
        prop_ids[r_elem_name][r_pos] = r_elem.GetProperties().Id();
        // Fill the connectivities.
        HDF5File::Matrix<int>& r_connectivities = connectivities[r_elem_name];
        const Element::GeometryType& r_geom = r_elem.GetGeometry();
        for (unsigned i = 0; i < r_geom.size(); ++i)
            r_connectivities(r_pos, i) = r_geom[i].Id();
        ++r_pos;
    }

    // Write to HDF5 file.
    for (const std::string& r_elem_name : mElementNames)
    {
        GetFile().WriteDataSet("/Elements/" + r_elem_name + "/Id", elem_ids[r_elem_name]);
        GetFile().WriteDataSet("/Elements/" + r_elem_name + "/PropertyId", prop_ids[r_elem_name]);
        GetFile().WriteDataSet("/Elements/" + r_elem_name + "/Connectivity", connectivities[r_elem_name]);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
