#include "custom_io/hdf5_model_part_io.h"

#include "custom_utilities/hdf5_points_data.h"
#include "custom_utilities/hdf5_connectivities_data.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIO::ModelPartIO(Parameters& rParams, File::Pointer pFile)
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

    mConditionPointers.resize(mConditionNames.size());
    for (unsigned i = 0; i < mConditionNames.size(); ++i)
    {
        const Condition& r_cond = KratosComponents<Condition>::Get(mConditionNames[i]);
        mConditionPointers[i] = &r_cond;
    }

    int num_threads = OpenMPUtils::GetNumThreads();
    mCurrentFillPosition.resize(num_threads, 0);

    KRATOS_CATCH("");
}

bool ModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();

    const std::size_t num_nodes = ReadNodesNumber();
    File& r_file = GetFile();
    PointsData points;

    points.ReadData(r_file, "/Nodes/Local", 0, num_nodes);
    points.CreateNodes(rNodes);

    return true;
    KRATOS_CATCH("");
}

std::size_t ModelPartIO::ReadNodesNumber()
{
    const std::vector<unsigned> dims = GetFile().GetDataDimensions("/Nodes/Local/Ids");
    return dims[0];
}

void ModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    PointsData points;
    points.SetData(rNodes);
    File& r_file = GetFile();
    points.WriteData(r_file, "/Nodes/Local");
    
    KRATOS_CATCH("");
}

void ModelPartIO::ReadElements(NodesContainerType& rNodes,
                               PropertiesContainerType& rProperties,
                               ElementsContainerType& rElements)
{
    KRATOS_TRY;

    rElements.clear();

    File& r_file = GetFile();
    ConnectivitiesData connectivities;

    for (unsigned i = 0; i < mElementNames.size(); ++i)
    {
        const std::string elem_path = "/Elements/" + mElementNames[i];
        const unsigned num_elems = r_file.GetDataDimensions(elem_path + "/Ids")[0];
        connectivities.Clear();
        connectivities.ReadData(r_file, elem_path, 0, num_elems);
        const Element& r_elem = *mElementPointers[i];
        connectivities.CreateElements(r_elem, rNodes, rProperties, rElements);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteElements(ElementsContainerType const& rElements)
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

void ModelPartIO::ReadConditions(NodesContainerType& rNodes,
                                 PropertiesContainerType& rProperties,
                                 ConditionsContainerType& rConditions)
{
    KRATOS_TRY;

    rConditions.clear();

    File& r_file = GetFile();
    ConnectivitiesData connectivities;

    for (unsigned i = 0; i < mConditionNames.size(); ++i)
    {
        const std::string cond_path = "/Conditions/" + mConditionNames[i];
        const unsigned num_conds = r_file.GetDataDimensions(cond_path + "/Ids")[0];
        connectivities.Clear();
        connectivities.ReadData(r_file, cond_path, 0, num_conds);
        const Condition& r_cond = *mConditionPointers[i];
        connectivities.CreateConditions(r_cond, rNodes, rProperties, rConditions);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    // Write the headers.
    for (const std::string& r_cond_name : mConditionNames)
        GetFile().AddPath("/Conditions/" + r_cond_name);

    if (mConditionNames.size() == 1)
        WriteUniformConditions(rConditions);
    else
        WriteMixedConditions(rConditions);

    KRATOS_CATCH("");
}

void ModelPartIO::ReadInitialValues(ModelPart& rModelPart)
{
}

void ModelPartIO::ReadInitialValues(NodesContainerType& rNodes,
                                        ElementsContainerType& rElements,
                                        ConditionsContainerType& rConditions)
{
}

void ModelPartIO::ReadModelPart(ModelPart& rModelPart)
{
}

void ModelPartIO::WriteModelPart(ModelPart& rModelPart)
{
}

HDF5::File& ModelPartIO::GetFile() const
{
    return *mpFile;
}

void ModelPartIO::Check()
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

void ModelPartIO::WriteUniformElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    const Element& r_expected_element = *mElementPointers[0];
    const Element& r_elem = rElements.front();
    KRATOS_ERROR_IF(typeid(r_elem) != typeid(r_expected_element))
        << "Elements are not of type: " << typeid(r_expected_element).name()
        << std::endl;

    ConnectivitiesData connectivities;
    connectivities.SetData(rElements);
    File& r_file = GetFile();
    connectivities.WriteData(r_file, "/Elements/" + mElementNames[0]);

    KRATOS_CATCH("");
}

void ModelPartIO::WriteMixedElements(ElementsContainerType const& rElements)
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
    std::vector<HDF5::File::Vector<int>> elem_ids(num_elem_types);
    std::vector<HDF5::File::Vector<int>> prop_ids(num_elem_types);
    std::vector<HDF5::File::Matrix<int>> connectivities(num_elem_types);
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
        HDF5::File::Matrix<int>& r_connectivities = connectivities[i_type];
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

void ModelPartIO::WriteUniformConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    const Condition& r_expected_condition = *mConditionPointers[0];
    const Condition& r_cond = rConditions.front();
    KRATOS_ERROR_IF(typeid(r_cond) != typeid(r_expected_condition))
        << "Conditions are not of type: " << typeid(r_expected_condition).name()
        << std::endl;

    ConnectivitiesData connectivities;
    connectivities.SetData(rConditions);
    File& r_file = GetFile();
    connectivities.WriteData(r_file, "/Conditions/" + mConditionNames[0]);

    KRATOS_CATCH("");
}

void ModelPartIO::WriteMixedConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    const unsigned num_cond_types = mConditionNames.size();
    std::vector<unsigned> num_conds(num_cond_types);
    std::vector<unsigned> num_nodes(num_cond_types);
    for (unsigned i_type = 0; i_type < num_cond_types; ++i_type)
    {
        num_conds[i_type] = 0; // Initialize count.
        num_nodes[i_type] = mConditionPointers[i_type]->GetGeometry().size();
    }

    // Count number of each condition and perform checks.
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rConditions.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        ConditionsContainerType::const_iterator it = rConditions.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            const Condition& r_cond = *it;
            unsigned i_type = FindIndexOfMatchingReferenceCondition(r_cond);
            // Check for matching node count.
            KRATOS_ERROR_IF(r_cond.GetGeometry().size() != num_nodes[i_type])
                << "Geometry has non-standard size for condition #"
                << r_cond.Id() << std::endl;
            // Increment condition count.
            unsigned& r_count = num_conds[i_type];
            #pragma omp atomic
            r_count += 1;
            ++it;
        }
    }

    // Initialized arrays.
    std::vector<unsigned> current_index(num_cond_types);
    std::vector<HDF5::File::Vector<int>> cond_ids(num_cond_types);
    std::vector<HDF5::File::Vector<int>> prop_ids(num_cond_types);
    std::vector<HDF5::File::Matrix<int>> connectivities(num_cond_types);
    for (unsigned i_type = 0; i_type < num_cond_types; ++i_type)
    {
        current_index[i_type] = 0; // Initialize position indices.
        cond_ids[i_type].resize(num_conds[i_type], false);
        prop_ids[i_type].resize(num_conds[i_type], false);
        connectivities[i_type].resize(num_conds[i_type], num_nodes[i_type], false);
    }

    // Fill arrays. By now all condition checks have been performed so no
    // exceptions should occur. Be careful with openmp here!!!
    for (const auto& r_cond : rConditions)
    {
        unsigned i_type = FindIndexOfMatchingReferenceCondition(r_cond);
        // Get the current position for this condition.
        unsigned pos = current_index[i_type];
        // Fill the condition and property ids.
        cond_ids[i_type][pos] = r_cond.Id();
        prop_ids[i_type][pos] = r_cond.GetProperties().Id();
        // Fill the connectivities.
        HDF5::File::Matrix<int>& r_connectivities = connectivities[i_type];
        const Condition::GeometryType& r_geom = r_cond.GetGeometry();
        for (unsigned k = 0; k < r_geom.size(); ++k)
            r_connectivities(pos, k) = r_geom[k].Id();
        current_index[i_type] += 1;
    }

    // Write to HDF5 file.
    for (unsigned i_type = 0; i_type < num_cond_types; ++i_type)
    {
        const std::string& r_cond_name = mConditionNames[i_type];
        GetFile().WriteDataSet("/Conditions/" + r_cond_name + "/Id", cond_ids[i_type]);
        GetFile().WriteDataSet("/Conditions/" + r_cond_name + "/PropertyId", prop_ids[i_type]);
        GetFile().WriteDataSet("/Conditions/" + r_cond_name + "/Connectivity", connectivities[i_type]);
    }

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
