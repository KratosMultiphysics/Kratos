#include "custom_io/hdf5_connectivities_data.h"

#include "includes/kratos_components.h"
#include "utilities/openmp_utils.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "custom_io/hdf5_file.h"

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
    rFile.ReadDataSet(rPath + "/Ids", mIds, StartIndex, BlockSize);
    rFile.ReadDataSet(rPath + "/PropertiesIds", mPropertiesIds, StartIndex, BlockSize);
    rFile.ReadDataSet(rPath + "/Connectivities", mConnectivities, StartIndex, BlockSize);
    if (KratosComponents<ElementType>::Has(mName))
    {
        const ElementType& r_elem = KratosComponents<ElementType>::Get(mName);
        KRATOS_ERROR_IF(r_elem.GetGeometry().size() != mConnectivities.size2())
            << "Element \"" << mName << "\""
            << " has a different connectivities size than the file data." << std::endl;
    }
    else if (KratosComponents<ConditionType>::Has(mName))
    {
        const ConditionType& r_cond = KratosComponents<ConditionType>::Get(mName);
        KRATOS_ERROR_IF(r_cond.GetGeometry().size() != mConnectivities.size2())
            << "Condition \"" << mName << "\""
            << " has a different connectivities size than the file data." << std::endl;
    }
    else
    {
        KRATOS_ERROR << "Invalid name \"" << mName << "\"!" << std::endl;
    }
    KRATOS_CATCH("");
}

void ConnectivitiesData::WriteData(File& rFile, const std::string& rPath, WriteInfo& rInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR_IF(mName == "") << "Cannot write an empty data." << std::endl;
    rFile.WriteDataSet(rPath + "/Ids", mIds, rInfo);
    rFile.WriteDataSet(rPath + "/PropertiesIds", mPropertiesIds, rInfo);
    rFile.WriteDataSet(rPath + "/Connectivities", mConnectivities, rInfo);
    int ws_dim, dim, num_nodes;
    if (KratosComponents<ElementType>::Has(mName))
    {
        const auto& r_geom = KratosComponents<ElementType>::Get(mName).GetGeometry();
        ws_dim = r_geom.WorkingSpaceDimension();
        dim = r_geom.Dimension();
        num_nodes = r_geom.size();
    }
    else
    {
        const auto& r_geom = KratosComponents<ConditionType>::Get(mName).GetGeometry();
        ws_dim = r_geom.WorkingSpaceDimension();
        dim = r_geom.Dimension();
        num_nodes = r_geom.size();
    }
    rFile.WriteAttribute(rPath, "Name", mName);
    rFile.WriteAttribute(rPath, "WorkingSpaceDimension", ws_dim);
    rFile.WriteAttribute(rPath, "Dimension", dim);
    rFile.WriteAttribute(rPath, "NumberOfNodes", num_nodes);
    KRATOS_CATCH("");
}

void ConnectivitiesData::CreateEntities(NodesContainerType& rNodes,
                                        PropertiesContainerType& rProperties,
                                        ElementsContainerType& rElements) const
{
    KRATOS_TRY;
    if (mName == "")
        return; // Do not append any new elements.
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
                                        ConditionsContainerType& rConditions) const
{
    KRATOS_TRY;
    if (mName == "")
        return; // Do not append any new conditions.
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

void ConnectivitiesData::SetData(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rElements.empty())
        << "Must pass non-empty container to this function." << std::endl;

    std::string name;
    CompareElementsAndConditionsUtility::GetRegisteredName(rElements.front(), name);
    SetData(name, rElements);

    KRATOS_CATCH("");
}

void ConnectivitiesData::SetData(const std::string& rName, ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(KratosComponents<ElementType>::Has(rName))
        << "Name \"" << rName << "\" is not registered as an element." << std::endl;

    if (rElements.empty())
    {
        Clear();
        // Correctly set the name and array sizes. Otherwise, MPI collective write
        // can fail due to inconsistent data across processes.
        mName = rName;
        const Element& r_elem = KratosComponents<ElementType>::Get(mName);
        mConnectivities.resize(0, r_elem.GetGeometry().size());
        return;
    }

    CompareElementsAndConditionsUtility::GetRegisteredName(rElements.front(), mName);
    KRATOS_ERROR_IF(mName != rName) << "Name \"" << rName << "\" is not the same as \"" << mName << "\"." << std::endl;
    const int num_elems = rElements.size();
    mIds.resize(num_elems, false);
    mPropertiesIds.resize(num_elems, false);
    const unsigned geometry_size = rElements.front().GetGeometry().size();
    mConnectivities.resize(num_elems, geometry_size, false);

    // Fill arrays and perform checks.
    #pragma omp parallel for
    for (int i = 0; i < num_elems; ++i)
    {
        ElementsContainerType::const_iterator it = rElements.begin() + i;
        // Check that the element and geometry types are the same.
        KRATOS_ERROR_IF_NOT(GeometricalObject::IsSame(*it, rElements.front()))
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

    KRATOS_ERROR_IF(rConditions.empty())
        << "Must pass non-empty container to this function." << std::endl;

    std::string name;
    CompareElementsAndConditionsUtility::GetRegisteredName(rConditions.front(), name);
    SetData(name, rConditions);

    KRATOS_CATCH("");
}

void ConnectivitiesData::SetData(const std::string& rName, ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(KratosComponents<ConditionType>::Has(rName))
        << "Name \"" << rName << "\" is not registered as a condition." << std::endl;

    if (rConditions.empty())
    {
        Clear();
        // Correctly set the name and array sizes. Otherwise, MPI collective write
        // can fail due to inconsistent data across processes.
        mName = rName;
        const Condition& r_cond = KratosComponents<ConditionType>::Get(mName);
        mConnectivities.resize(0, r_cond.GetGeometry().size());
        return;
    }

    CompareElementsAndConditionsUtility::GetRegisteredName(rConditions.front(), mName);
    KRATOS_ERROR_IF(mName != rName) << "Name \"" << rName << "\" is not the same as \"" << mName << "\"." << std::endl;
    const int num_conds = rConditions.size();
    mIds.resize(num_conds, false);
    mPropertiesIds.resize(num_conds, false);
    const unsigned geometry_size = rConditions.front().GetGeometry().size();
    mConnectivities.resize(num_conds, geometry_size, false);

    // Fill arrays and perform checks.
    #pragma omp parallel for
    for (int i = 0; i < num_conds; ++i)
    {
        ConditionsContainerType::const_iterator it = rConditions.begin() + i;
        // Check that the condition and geometry types are the same.
        KRATOS_ERROR_IF_NOT(GeometricalObject::IsSame(*it, rConditions.front()))
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
    mIds.resize(0);
    mPropertiesIds.resize(0);
    mConnectivities.resize(0,0);
}
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
