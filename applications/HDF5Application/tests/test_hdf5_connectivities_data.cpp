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
#include "utilities/compare_elements_and_conditions_utility.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_io/hdf5_file_serial.h"
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_ConnectivitiesData_CreateEntities1, KratosHDF5TestSuite)
{
    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ElementsContainerType elements;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, elements, ids, pids, connectivities);
    Parameters test_params(R"(
        {
            "file_name": "test.h5",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::WriteInfo info;
    test_file.WriteDataSet("/Elements/Ids", ids, info);
    test_file.WriteDataSet("/Elements/PropertiesIds", pids, info);
    test_file.WriteDataSet("/Elements/Connectivities", connectivities, info);
    std::string name;
    CompareElementsAndConditionsUtility::GetRegisteredName(elements.front(), name);
    test_file.WriteAttribute("/Elements", "Name", name);
    HDF5::Internals::ConnectivitiesData data;
    data.ReadData(test_file, "/Elements", info.StartIndex, info.BlockSize);
    HDF5::ElementsContainerType new_elements;
    data.CreateEntities(nodes, properties, new_elements);
    KRATOS_CHECK(new_elements.size() == elements.size());
    for (Element& r_new_elem : new_elements)
    {
        Element& r_elem = elements[r_new_elem.Id()];
        KRATOS_CHECK(r_new_elem.GetProperties().Id() == r_elem.GetProperties().Id());
        KRATOS_CHECK(GeometricalObject::IsSame(r_new_elem, r_elem));
        for (std::size_t j = 0; j < r_new_elem.GetGeometry().size(); ++j)
            KRATOS_CHECK(r_new_elem.GetGeometry()[j].Id() == r_elem.GetGeometry()[j].Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ConnectivitiesData_CreateEntities2, KratosHDF5TestSuite)
{
    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ConditionsContainerType conditions;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, conditions, ids, pids, connectivities);
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::WriteInfo info;
    test_file.WriteDataSet("/Conditions/Ids", ids, info);
    test_file.WriteDataSet("/Conditions/PropertiesIds", pids, info);
    test_file.WriteDataSet("/Conditions/Connectivities", connectivities, info);
    std::string name;
    CompareElementsAndConditionsUtility::GetRegisteredName(conditions.front(), name);
    test_file.WriteAttribute("/Conditions", "Name", name);
    HDF5::Internals::ConnectivitiesData data;
    data.ReadData(test_file, "/Conditions", info.StartIndex, info.BlockSize);
    HDF5::ConditionsContainerType new_conditions;
    data.CreateEntities(nodes, properties, new_conditions);
    KRATOS_CHECK(new_conditions.size() == conditions.size());
    for (Condition& r_new_cond : new_conditions)
    {
        Condition& r_cond = conditions[r_new_cond.Id()];
        KRATOS_CHECK(r_new_cond.GetProperties().Id() == r_cond.GetProperties().Id());
        KRATOS_CHECK(GeometricalObject::IsSame(r_new_cond, r_cond));
        for (std::size_t j = 0; j < r_cond.GetGeometry().size(); ++j)
            KRATOS_CHECK(r_new_cond.GetGeometry()[j].Id() == r_cond.GetGeometry()[j].Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_WriteData1, KratosHDF5TestSuite)
{
    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ElementsContainerType elements;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, elements, ids, pids, connectivities);
    HDF5::Internals::ConnectivitiesData data;
    data.SetData(FactorElements(elements).front());
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::WriteInfo info;
    data.WriteData(test_file, "/Elements", info);
    data.Clear();
    data.ReadData(test_file, "/Elements", info.StartIndex, info.BlockSize);
    HDF5::ElementsContainerType new_elements;
    data.CreateEntities(nodes, properties, new_elements);
    KRATOS_CHECK(new_elements.size() == elements.size());
    for (Element& r_new_elem : new_elements)
    {
        Element& r_elem = elements[r_new_elem.Id()];
        KRATOS_CHECK(r_new_elem.GetProperties().Id() == r_elem.GetProperties().Id());
        KRATOS_CHECK(GeometricalObject::IsSame(r_new_elem, r_elem));
        for (std::size_t j = 0; j < r_elem.GetGeometry().size(); ++j)
            KRATOS_CHECK(r_new_elem.GetGeometry()[j].Id() == r_elem.GetGeometry()[j].Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5ConnectivitiesData_WriteData2, KratosHDF5TestSuite)
{
    HDF5::NodesContainerType nodes;
    HDF5::PropertiesContainerType properties;
    HDF5::ConditionsContainerType conditions;
    HDF5::Vector<int> ids, pids;
    HDF5::Matrix<int> connectivities;
    CreateTestMesh(nodes, properties, conditions, ids, pids, connectivities);
    HDF5::Internals::ConnectivitiesData data;
    data.SetData(FactorConditions(conditions).front());
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::WriteInfo info;
    data.WriteData(test_file, "/Conditions", info);
    data.Clear();
    data.ReadData(test_file, "/Conditions", info.StartIndex, info.BlockSize);
    HDF5::ConditionsContainerType new_conditions;
    data.CreateEntities(nodes, properties, new_conditions);
    KRATOS_CHECK(new_conditions.size() == conditions.size());
    for (Condition& r_new_cond : new_conditions)
    {
        Condition& r_cond = conditions[r_new_cond.Id()];
        KRATOS_CHECK(r_new_cond.GetProperties().Id() == r_cond.GetProperties().Id());
        KRATOS_CHECK(GeometricalObject::IsSame(r_new_cond, r_cond));
        for (std::size_t j = 0; j < r_cond.GetGeometry().size(); ++j)
            KRATOS_CHECK(r_new_cond.GetGeometry()[j].Id() == r_cond.GetGeometry()[j].Id());
    }
}

} // namespace Testing
} // namespace Kratos.
