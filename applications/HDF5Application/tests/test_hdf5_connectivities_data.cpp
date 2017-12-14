//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/node.h"
#include "includes/kratos_components.h"

// Application includes
#include "custom_io/hdf5_file_serial.h"
#include "custom_utilities/hdf5_connectivities_data.h"
#include "custom_utilities/hdf5_pointer_bins_utility.h"

namespace Kratos
{
namespace Testing
{
void CreateTestMesh(HDF5::NodesContainerType& rNodes,
                    HDF5::PropertiesContainerType& rProperties,
                    HDF5::ElementsContainerType& rElements,
                    HDF5::Vector<int>& rElementIds,
                    HDF5::Vector<int>& rPropertiesIds,
                    HDF5::Matrix<int>& rConnectivities)
{
    const unsigned num_elems = 2;
    const unsigned num_nodes = num_elems + 3;

    rNodes.clear();
    rProperties.clear();
    rElements.clear();
    rElementIds.resize(num_elems, false);
    rPropertiesIds.resize(num_elems, false);
    rConnectivities.resize(num_elems, 3, false);

    const HDF5::ElementType& Element2D3N = KratosComponents<HDF5::ElementType>::Get("Element2D3N");
    // Create nodes.
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        auto p_node = boost::make_shared<HDF5::NodeType>(
            i + 1, 1.0, 2.0, 3.0);
        rNodes.push_back(p_node);
    }
    // Create elements.
    Element::NodesArrayType geom_nodes(3);
    for (unsigned i = 0; i < num_elems; ++i)
    {
        rElementIds[i] = i + 1;
        rPropertiesIds[i] = 1;
        for (unsigned j = 0; j < 3; ++j)
        {
            rConnectivities(i, j) = i + j + 1;
            geom_nodes(j) = rNodes(i + j + 1);
        }
        auto p_elem = Element2D3N.Create(i + 1, geom_nodes, rProperties(1));
        rElements.push_back(p_elem);
    }
}

void CreateTestMesh(HDF5::NodesContainerType& rNodes,
                    HDF5::PropertiesContainerType& rProperties,
                    HDF5::ConditionsContainerType& rConditions,
                    HDF5::Vector<int>& rConditionIds,
                    HDF5::Vector<int>& rPropertiesIds,
                    HDF5::Matrix<int>& rConnectivities)
{
    const unsigned num_conds = 2;
    const unsigned num_nodes = num_conds + 3;

    rNodes.clear();
    rProperties.clear();
    rConditions.clear();
    rConditionIds.resize(num_conds, false);
    rPropertiesIds.resize(num_conds, false);
    rConnectivities.resize(num_conds, 3, false);

    const HDF5::ConditionType& SurfaceCondition3D3N = KratosComponents<HDF5::ConditionType>::Get("SurfaceCondition3D3N");
    // Create nodes.
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        auto p_node = boost::make_shared<HDF5::NodeType>(
            i + 1, 1.0, 2.0, 3.0);
        rNodes.push_back(p_node);
    }
    // Create conditions.
    Condition::NodesArrayType geom_nodes(3);
    for (unsigned i = 0; i < num_conds; ++i)
    {
        rConditionIds[i] = i + 1;
        rPropertiesIds[i] = 1;
        for (unsigned j = 0; j < 3; ++j)
        {
            rConnectivities(i, j) = i + j + 1;
            geom_nodes(j) = rNodes(i + j + 1);
        }
        auto p_cond = SurfaceCondition3D3N.Create(i + 1, geom_nodes, rProperties(1));
        rConditions.push_back(p_cond);
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_ReadData, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
        {
            "file_name": "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);

    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ElementsContainerType elements;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, elements, ids, pids, connectivities);

    test_file.WriteDataSet("/Elements/Ids", ids);
    test_file.WriteDataSet("/Elements/PropertiesIds", pids);
    test_file.WriteDataSet("/Elements/Connectivities", connectivities);

    HDF5::Internals::ConnectivitiesData data;
    data.ReadData(test_file, "/Elements", 0, ids.size());
    KRATOS_CHECK(data.size() == ids.size());
    KRATOS_CHECK(data.GetConnectivities().size2() == connectivities.size2());
    for (unsigned i = 0; i < data.size(); ++i)
    {
        KRATOS_CHECK(data.GetIds()[i] == ids[i]);
        KRATOS_CHECK(data.GetPropertiesIds()[i] == pids[i]);
        for (unsigned j = 0; j < connectivities.size2(); ++j)
            KRATOS_CHECK(data.GetConnectivities()(i,j) == connectivities(i,j));
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_CreateElements, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);

    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ElementsContainerType elements;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, elements, ids, pids, connectivities);

    test_file.WriteDataSet("/Elements/Ids", ids);
    test_file.WriteDataSet("/Elements/PropertiesIds", pids);
    test_file.WriteDataSet("/Elements/Connectivities", connectivities);

   const HDF5::ElementType& Element2D3N = KratosComponents<HDF5::ElementType>::Get("Element2D3N");

    HDF5::Internals::ConnectivitiesData data;
    data.ReadData(test_file, "/Elements", 0, ids.size());
    HDF5::ElementsContainerType new_elements;
    data.CreateEntities(Element2D3N, nodes, properties, new_elements);

    KRATOS_CHECK(new_elements.size() == elements.size());
    for (Element& r_new_elem : new_elements)
    {
        Element& r_elem = elements[r_new_elem.Id()];
        KRATOS_CHECK(r_new_elem.GetProperties().Id() == r_elem.GetProperties().Id());
        KRATOS_CHECK(r_new_elem.GetGeometry().size() == r_elem.GetGeometry().size());
        for (unsigned j = 0; j < r_elem.GetGeometry().size(); ++j)
            KRATOS_CHECK(r_new_elem.GetGeometry()[j].Id() == r_elem.GetGeometry()[j].Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_CreateConditions, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);

    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ConditionsContainerType conditions;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, conditions, ids, pids, connectivities);

    test_file.WriteDataSet("/Conditions/Ids", ids);
    test_file.WriteDataSet("/Conditions/PropertiesIds", pids);
    test_file.WriteDataSet("/Conditions/Connectivities", connectivities);

   const HDF5::ConditionType& SurfaceCondition3D3N = KratosComponents<HDF5::ConditionType>::Get("SurfaceCondition3D3N");

    HDF5::Internals::ConnectivitiesData data;
    data.ReadData(test_file, "/Conditions", 0, ids.size());
    HDF5::ConditionsContainerType new_conditions;
    data.CreateEntities(SurfaceCondition3D3N, nodes, properties, new_conditions);

    KRATOS_CHECK(new_conditions.size() == conditions.size());
    for (Condition& r_new_cond : new_conditions)
    {
        Condition& r_cond = conditions[r_new_cond.Id()];
        KRATOS_CHECK(r_new_cond.GetProperties().Id() == r_cond.GetProperties().Id());
        KRATOS_CHECK(r_new_cond.GetGeometry().size() == r_cond.GetGeometry().size());
        for (unsigned j = 0; j < r_cond.GetGeometry().size(); ++j)
            KRATOS_CHECK(r_new_cond.GetGeometry()[j].Id() == r_cond.GetGeometry()[j].Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_SetData1, KratosHDF5TestSuite)
{
    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ElementsContainerType elements;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, elements, ids, pids, connectivities);

    HDF5::Internals::ConnectivitiesData data;
    std::vector<const HDF5::ElementType*> bin_keys{&elements.front()};
    HDF5::Internals::PointerBinsUtility<HDF5::ElementType> elem_bins(bin_keys);
    elem_bins.CreateBins(elements);
    HDF5::ConstElementsContainerType& element_ptrs = elem_bins.GetBin(bin_keys.front());
    data.SetData(element_ptrs);
    KRATOS_CHECK(data.size() == element_ptrs.size());
    for (unsigned i = 0; i < data.size(); ++i)
    {
        Element const& r_elem = *element_ptrs[i];
        KRATOS_CHECK(r_elem.Id() == static_cast<unsigned>(data.GetIds()[i]));
        KRATOS_CHECK(r_elem.GetProperties().Id() == static_cast<unsigned>(data.GetPropertiesIds()[i]));
        KRATOS_CHECK(r_elem.GetGeometry().size() == static_cast<unsigned>(data.GetConnectivities().size2()));
        for (unsigned j = 0; j < data.GetConnectivities().size2(); ++j)
            KRATOS_CHECK(r_elem.GetGeometry()[j].Id() == static_cast<unsigned>(data.GetConnectivities()(i, j)));
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_SetData2, KratosHDF5TestSuite)
{
    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ConditionsContainerType conditions;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, conditions, ids, pids, connectivities);

    HDF5::Internals::ConnectivitiesData data;
    std::vector<const HDF5::ConditionType*> bin_keys{&conditions.front()};
    HDF5::Internals::PointerBinsUtility<HDF5::ConditionType> cond_bins(bin_keys);
    cond_bins.CreateBins(conditions);
    HDF5::ConstConditionsContainerType& condition_ptrs = cond_bins.GetBin(bin_keys.front());
    data.SetData(condition_ptrs);
    KRATOS_CHECK(data.size() == conditions.size());
    for (unsigned i = 0; i < data.size(); ++i)
    {
        Condition const& r_cond = *condition_ptrs[i];
        KRATOS_CHECK(r_cond.Id() == static_cast<unsigned>(data.GetIds()[i]));
        KRATOS_CHECK(r_cond.GetProperties().Id() == static_cast<unsigned>(data.GetPropertiesIds()[i]));
        KRATOS_CHECK(r_cond.GetGeometry().size() == static_cast<unsigned>(data.GetConnectivities().size2()));
        for (unsigned j = 0; j < data.GetConnectivities().size2(); ++j)
            KRATOS_CHECK(r_cond.GetGeometry()[j].Id() == static_cast<unsigned>(data.GetConnectivities()(i, j)));
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_WriteData, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
        {
            "file_name": "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);

    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ElementsContainerType elements;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, elements, ids, pids, connectivities);

    HDF5::Internals::ConnectivitiesData data;
    std::vector<const HDF5::ElementType*> bin_keys{&elements.front()};
    HDF5::Internals::PointerBinsUtility<HDF5::ElementType> elem_bins(bin_keys);
    elem_bins.CreateBins(elements);
    HDF5::ConstElementsContainerType& element_ptrs = elem_bins.GetBin(bin_keys.front());
    data.SetData(element_ptrs);
    data.WriteData(test_file, "/Elements");

    HDF5::Vector<int> new_ids, new_pids;
    HDF5::Matrix<int> new_connectivities;
    test_file.ReadDataSet("/Elements/Ids", new_ids, 0, ids.size());
    test_file.ReadDataSet("/Elements/PropertiesIds", new_pids, 0, pids.size());
    test_file.ReadDataSet("/Elements/Connectivities", new_connectivities, 0, connectivities.size1());
    KRATOS_CHECK(new_connectivities.size2() == connectivities.size2());
    for (unsigned i = 0; i < ids.size(); ++i)
    {
        KRATOS_CHECK(new_ids[i] == ids[i]);
        KRATOS_CHECK(new_pids[i] == pids[i]);
        for (unsigned j = 0; j < connectivities.size2(); ++j)
            KRATOS_CHECK(new_connectivities(i, j) == connectivities(i, j));
    }
}

} // namespace Testing
} // namespace Kratos.
