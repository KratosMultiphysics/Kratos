#include "custom_utilities/hdf5_connectivities_data.h"

#include "utilities/openmp_utils.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
void ConnectivitiesData::ReadData(File& rFile, std::string Path, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    Clear();
    rFile.ReadDataSet(Path + "/Ids", mIds, StartIndex, BlockSize);
    rFile.ReadDataSet(Path + "/PropertiesIds", mPropertiesIds, StartIndex, BlockSize);
    rFile.ReadDataSet(Path + "/Connectivities", mConnectivities, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void ConnectivitiesData::WriteData(File& rFile, std::string Path)
{
    KRATOS_TRY;
    rFile.WriteDataSet(Path + "/Ids", mIds);
    rFile.WriteDataSet(Path + "/PropertiesIds", mPropertiesIds);
    rFile.WriteDataSet(Path + "/Connectivities", mConnectivities);
    KRATOS_CATCH("");
}

void ConnectivitiesData::CreateEntities(ElementType const& rElementType,
                                        NodesContainerType& rNodes,
                                        PropertiesContainerType& rProperties,
                                        ElementsContainerType& rElements)
{
    KRATOS_TRY;
    const unsigned num_new_elems = mIds.size();
    const unsigned geometry_size = rElementType.GetGeometry().size();
    KRATOS_ERROR_IF(geometry_size != mConnectivities.size2())
        << "Non-matching geometry and connectivity sizes." << std::endl;
    rElements.reserve(rElements.size() + num_new_elems);
    ElementType::NodesArrayType nodes(geometry_size);

    for (unsigned i = 0; i < num_new_elems; ++i)
    {
        for (unsigned j = 0; j < geometry_size; ++j)
        {
            const int node_id = mConnectivities(i, j);
            nodes(j) = rNodes(node_id);
        }
        ElementType::Pointer p_elem =
            rElementType.Create(mIds[i], nodes, rProperties(mPropertiesIds[i]));
        rElements.push_back(p_elem);
    }
    KRATOS_CATCH("");
}

void ConnectivitiesData::CreateEntities(ConditionType const& rConditionType,
                                        NodesContainerType& rNodes,
                                        PropertiesContainerType& rProperties,
                                        ConditionsContainerType& rConditions)
{
    KRATOS_TRY;
    const unsigned num_new_conds = mIds.size();
    const unsigned geometry_size = rConditionType.GetGeometry().size();
    KRATOS_ERROR_IF(geometry_size != mConnectivities.size2())
        << "Non-matching geometry and connectivity sizes." << std::endl;
    rConditions.reserve(rConditions.size() + num_new_conds);
    ConditionType::NodesArrayType nodes(geometry_size);

    for (unsigned i = 0; i < num_new_conds; ++i)
    {
        for (unsigned j = 0; j < geometry_size; ++j)
        {
            const int node_id = mConnectivities(i, j);
            nodes(j) = rNodes(node_id);
        }
        Condition::Pointer p_cond =
            rConditionType.Create(mIds[i], nodes, rProperties(mPropertiesIds[i]));
        rConditions.push_back(p_cond);
    }
    KRATOS_CATCH("");
}

// Fill data from elements of a single element type.
void ConnectivitiesData::SetData(ConstElementsContainerType const& rElements)
{
    KRATOS_TRY;

    const unsigned num_elems = rElements.size();
    if (num_elems == 0)
    {
        Clear();
        return;
    }
    const ElementType& r_expected_element = *rElements.front();
    const unsigned geometry_size = r_expected_element.GetGeometry().size();
    mIds.resize(num_elems, false);
    mPropertiesIds.resize(num_elems, false);
    mConnectivities.resize(num_elems, geometry_size, false);

    // Fill arrays and perform checks.
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(num_elems, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            const ElementType& r_elem = *rElements[i];
            // Check element type. (this is done by the bins utility.)
            //KRATOS_ERROR_IF(typeid(r_elem) != typeid(r_expected_element))
            //    << "Element #" << r_elem.Id() << " is not "
            //    << typeid(r_expected_element).name() << std::endl;
            // Fill ids.
            mIds[i] = r_elem.Id();
            mPropertiesIds[i] = r_elem.GetProperties().Id();
            // Check for matching node count.
            const ElementType::GeometryType& r_geom = r_elem.GetGeometry();
            KRATOS_ERROR_IF(r_geom.size() != geometry_size)
                << "Geometry has non-standard size for element #" << r_elem.Id()
                << std::endl;
            // Fill connectivities.
            for (unsigned k = 0; k < geometry_size; ++k)
                mConnectivities(i, k) = r_geom[k].Id();
        }
    }

    KRATOS_CATCH("");
}

// Fill data from conditions of a single condition type.
void ConnectivitiesData::SetData(ConstConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    const unsigned num_conds = rConditions.size();
    if (num_conds == 0)
    {
        Clear();
        return;
    }
    const ConditionType& r_expected_condition = *rConditions.front();
    const unsigned geometry_size = r_expected_condition.GetGeometry().size();
    mIds.resize(num_conds, false);
    mPropertiesIds.resize(num_conds, false);
    mConnectivities.resize(num_conds, geometry_size, false);

    // Fill arrays and perform checks.
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(num_conds, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            const ConditionType& r_cond = *rConditions[i];
            // Check condition type. (this is done by the bins utility.)
            //KRATOS_ERROR_IF(typeid(r_cond) != typeid(r_expected_condition))
            //    << "Condition #" << r_cond.Id() << " is not "
            //    << typeid(r_expected_condition).name() << std::endl;
            // Fill ids.
            mIds[i] = r_cond.Id();
            mPropertiesIds[i] = r_cond.GetProperties().Id();
            // Check for matching node count.
            const ConditionType::GeometryType& r_geom = r_cond.GetGeometry();
            KRATOS_ERROR_IF(r_geom.size() != geometry_size)
                << "Geometry has non-standard size for condition #" << r_cond.Id()
                << std::endl;
            // Fill connectivities.
            for (unsigned k = 0; k < geometry_size; ++k)
                mConnectivities(i, k) = r_geom[k].Id();
        }
    }

    KRATOS_CATCH("");
}

void ConnectivitiesData::Clear()
{
    mIds.clear();
    mPropertiesIds.clear();
    mConnectivities.clear();
}
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
