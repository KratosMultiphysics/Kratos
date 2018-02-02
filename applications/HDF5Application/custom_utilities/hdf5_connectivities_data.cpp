#include "custom_utilities/hdf5_connectivities_data.h"

#include "includes/kratos_components.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/compare_element_and_condition_utility.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
void ConnectivitiesData::ReadData(File& rFile, const std::string& rPath, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    Clear();
    rFile.ReadAttribute(rPath, "Name", mName);
    if (KratosComponents<ElementType>::Has(mName) == false && KratosComponents<ConditionType>::Has(mName) == false)
        KRATOS_ERROR << "Invalid name \"" << mName << "\"!" << std::endl;
    rFile.ReadDataSet(rPath + "/Ids", mIds, StartIndex, BlockSize);
    rFile.ReadDataSet(rPath + "/PropertiesIds", mPropertiesIds, StartIndex, BlockSize);
    rFile.ReadDataSet(rPath + "/Connectivities", mConnectivities, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void ConnectivitiesData::WriteData(File& rFile, const std::string& rPath, WriteInfo& rInfo)
{
    KRATOS_TRY;
    rFile.WriteDataSet(rPath + "/Ids", mIds, rInfo);
    rFile.WriteDataSet(rPath + "/PropertiesIds", mPropertiesIds, rInfo);
    rFile.WriteDataSet(rPath + "/Connectivities", mConnectivities, rInfo);
    int ws_dim, dim, nnodes;
    if (KratosComponents<ElementType>::Has(mName))
    {
        const auto& r_elem = KratosComponents<ElementType>::Get(mName);
        ws_dim = r_elem.WorkingSpaceDimension();
        dim = r_elem.GetGeometry().Dimension();
        nnodes = r_elem.GetGeometry().size();
    }
    else
    {
        const auto& r_cond = KratosComponents<ConditionType>::Get(mName);
        ws_dim = r_cond.WorkingSpaceDimension();
        dim = r_cond.GetGeometry().Dimension();
        nnodes = r_cond.GetGeometry().size();
    }
    rFile.WriteAttribute(rPath, "Name", mName);
    rFile.WriteAttribute(rPath, "WorkingSpaceDimension", ws_dim);
    rFile.WriteAttribute(rPath, "Dimension", dim);
    rFile.WriteAttribute(rPath, "NumberOfNodes", nnodes);
    KRATOS_CATCH("");
}

void ConnectivitiesData::CreateEntities(NodesContainerType& rNodes,
                                        PropertiesContainerType& rProperties,
                                        ElementsContainerType& rElements)
{
    KRATOS_TRY;
    const ElementType& r_elem = KratosComponents<ElementType>::Get(mName);
    const unsigned num_new_elems = mIds.size();
    const unsigned geometry_size = r_elem.GetGeometry().size();
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
            r_elem.Create(mIds[i], nodes, rProperties(mPropertiesIds[i]));
        rElements.push_back(p_elem);
    }
    KRATOS_CATCH("");
}

void ConnectivitiesData::CreateEntities(NodesContainerType& rNodes,
                                        PropertiesContainerType& rProperties,
                                        ConditionsContainerType& rConditions)
{
    KRATOS_TRY;
    const ConditionType& r_cond = KratosComponents<ConditionType>::Get(mName);
    const unsigned num_new_conds = mIds.size();
    const unsigned geometry_size = r_cond.GetGeometry().size();
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
            r_cond.Create(mIds[i], nodes, rProperties(mPropertiesIds[i]));
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
    CompareElementAndConditionUtility::FindNameInRegistry(r_expected_element, mName);
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
    CompareElementAndConditionUtility::FindNameInRegistry(r_expected_condition, mName);
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
    mName = "";
    mIds.clear();
    mPropertiesIds.clear();
    mConnectivities.clear();
}
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
