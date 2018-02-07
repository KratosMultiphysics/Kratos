#include "custom_io/hdf5_connectivities_data.h"

#include "includes/kratos_components.h"
#include "utilities/openmp_utils.h"
#include "utilities/compare_elements_and_conditions_utility.h"

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
void ConnectivitiesData::SetData(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    const unsigned num_elems = rElements.size();
    if (num_elems == 0)
    {
        Clear();
        return;
    }
    mIds.resize(num_elems, false);
    mPropertiesIds.resize(num_elems, false);
    const unsigned geometry_size = rElements.front().GetGeometry().size();
    mConnectivities.resize(num_elems, geometry_size, false);
    CompareElementsAndConditionsUtility::FindNameInRegistry(rElements.front(), mName);

    // Fill arrays and perform checks.
    #pragma omp parallel for
    for (unsigned i = 0; i < num_elems; ++i)
    {
        ElementsContainerType::const_iterator it = rElements.begin() + i;
        // Check that the element and geometry types are the same.
        KRATOS_ERROR_IF_NOT(CompareElementsAndConditionsUtility::IsSame(*it, rElements.front()))
            << "Element #" << it->Id() << " is not the same as #" << rElements.front().Id() << '!' << std::endl;
        // Fill ids.
        mIds[i] = it->Id();
        mPropertiesIds[i] = it->GetProperties().Id();
        const ElementType::GeometryType& r_geom = it->GetGeometry();
        // Fill connectivities.
        for (unsigned k = 0; k < geometry_size; ++k)
            mConnectivities(i, k) = r_geom[k].Id();
    }

    KRATOS_CATCH("");
}

// Fill data from conditions of a single condition type.
void ConnectivitiesData::SetData(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    const unsigned num_conds = rConditions.size();
    if (num_conds == 0)
    {
        Clear();
        return;
    }
    mIds.resize(num_conds, false);
    mPropertiesIds.resize(num_conds, false);
    const unsigned geometry_size = rConditions.front().GetGeometry().size();
    mConnectivities.resize(num_conds, geometry_size, false);
    CompareElementsAndConditionsUtility::FindNameInRegistry(rConditions.front(), mName);

    // Fill arrays and perform checks.
    #pragma omp parallel for
    for (unsigned i = 0; i < num_conds; ++i)
    {
        ConditionsContainerType::const_iterator it = rConditions.begin() + i;
        // Check that the condition and geometry types are the same.
        KRATOS_ERROR_IF_NOT(CompareElementsAndConditionsUtility::IsSame(*it, rConditions.front()))
            << "Condition #" << it->Id() << " is not the same as #" << rConditions.front().Id() << '!' << std::endl;
        // Fill ids.
        mIds[i] = it->Id();
        mPropertiesIds[i] = it->GetProperties().Id();
        const ConditionType::GeometryType& r_geom = it->GetGeometry();
        // Fill connectivities.
        for (unsigned k = 0; k < geometry_size; ++k)
            mConnectivities(i, k) = r_geom[k].Id();
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
