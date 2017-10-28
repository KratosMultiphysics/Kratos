#include "hdf5_model_part_io.h"

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

    mElementPointers.resize(mElementNames.size());
    for (unsigned i = 0; i < mElementNames.size(); ++i)
    {
        const Element& r_elem = KratosComponents<Element>::Get(mElementNames[i]);
        mElementPointers[i] = &r_elem;
    }

    int num_threads = OpenMPUtils::GetNumThreads();
    mCurrentElement.resize(num_threads, 0);

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

    const Element& r_expected_element = *mElementPointers[0];
    unsigned size1 = rElements.size();
    unsigned size2 = r_expected_element.GetGeometry().size();
    HDF5File::Vector<int> elem_ids(size1);
    HDF5File::Vector<int> prop_ids(size1);
    HDF5File::Matrix<int> connectivities(size1, size2);

    // Fill arrays and perform checks.
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rElements.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        ElementsContainerType::const_iterator it = rElements.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            const Element& r_elem = *it;
            // Check element type.
            KRATOS_ERROR_IF(typeid(r_elem) != typeid(r_expected_element))
                << "Element #" << r_elem.Id() << " is not "
                << typeid(r_expected_element).name() << std::endl;
            // Fill ids.
            elem_ids[i] = r_elem.Id();
            prop_ids[i] = r_elem.GetProperties().Id();
            // Check for matching node count.
            const Element::GeometryType& r_geom = r_elem.GetGeometry();
            KRATOS_ERROR_IF(r_geom.size() != size2)
                << "Geometry has non-standard size for element #" << r_elem.Id()
                << std::endl;
            // Fill connectivities.
            for (unsigned k = 0; k < size2; ++k)
                connectivities(i, k) = r_geom[k].Id();
            ++it;
        }
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

    const unsigned num_elem_types = mElementNames.size();
    std::vector<unsigned> num_elems(num_elem_types);
    std::vector<unsigned> num_nodes(num_elem_types);
    for (unsigned i_type = 0; i_type < num_elem_types; ++i_type)
    {
        num_elems[i_type] = 0; // Initialize count.
        num_nodes[i_type] = mElementPointers[i_type]->GetGeometry().size();
    }

    // Count number of each element and perform checks.
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rElements.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        ElementsContainerType::const_iterator it = rElements.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            const Element& r_elem = *it;
            unsigned i_type = FindIndexOfMatchingReferenceElement(r_elem);
            // Check for matching node count.
            KRATOS_ERROR_IF(r_elem.GetGeometry().size() != num_nodes[i_type])
                << "Geometry has non-standard size for element #"
                << r_elem.Id() << std::endl;
            // Increment element count.
            unsigned& r_count = num_elems[i_type];
            #pragma omp atomic
            r_count += 1;
            ++it;
        }
    }

    // Initialized arrays.
    std::vector<unsigned> current_index(num_elem_types);
    std::vector<HDF5File::Vector<int>> elem_ids(num_elem_types);
    std::vector<HDF5File::Vector<int>> prop_ids(num_elem_types);
    std::vector<HDF5File::Matrix<int>> connectivities(num_elem_types);
    for (unsigned i_type = 0; i_type < num_elem_types; ++i_type)
    {
        current_index[i_type] = 0; // Initialize position indices.
        elem_ids[i_type].resize(num_elems[i_type], false);
        prop_ids[i_type].resize(num_elems[i_type], false);
        connectivities[i_type].resize(num_elems[i_type], num_nodes[i_type], false);
    }

    // Fill arrays. By now all element checks have been performed so no
    // exceptions should occur. Be careful with openmp here!!!
    for (const auto& r_elem : rElements)
    {
        unsigned i_type = FindIndexOfMatchingReferenceElement(r_elem);
        // Get the current position for this element.
        unsigned pos = current_index[i_type];
        // Fill the element and property ids.
        elem_ids[i_type][pos] = r_elem.Id();
        prop_ids[i_type][pos] = r_elem.GetProperties().Id();
        // Fill the connectivities.
        HDF5File::Matrix<int>& r_connectivities = connectivities[i_type];
        const Element::GeometryType& r_geom = r_elem.GetGeometry();
        for (unsigned k = 0; k < r_geom.size(); ++k)
            r_connectivities(pos, k) = r_geom[k].Id();
        current_index[i_type] += 1;
    }

    // Write to HDF5 file.
    for (unsigned i_type = 0; i_type < num_elem_types; ++i_type)
    {
        const std::string& r_elem_name = mElementNames[i_type];
        GetFile().WriteDataSet("/Elements/" + r_elem_name + "/Id", elem_ids[i_type]);
        GetFile().WriteDataSet("/Elements/" + r_elem_name + "/PropertyId", prop_ids[i_type]);
        GetFile().WriteDataSet("/Elements/" + r_elem_name + "/Connectivity", connectivities[i_type]);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
