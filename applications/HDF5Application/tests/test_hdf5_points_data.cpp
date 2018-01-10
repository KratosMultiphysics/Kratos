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
#include "custom_io/hdf5_file_serial.h"
#include "custom_utilities/hdf5_points_data.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_ReadData, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::Vector<int> ids(3);
    HDF5::Vector<array_1d<double, 3>> coords(3);
    for (unsigned i = 0; i < 3; ++i)
    {
        ids[i] = i + 1;
        coords[i] = array_1d<double, 3>(3, 1.2345);
    }
    test_file.WriteDataSet("/Nodes/Ids", ids);
    test_file.WriteDataSet("/Nodes/Coordinates", coords);
    HDF5::Internals::PointsData data;
    data.ReadData(test_file, "/Nodes", 0, 3);
    for (unsigned i = 0; i < 3; ++i)
    {
        KRATOS_CHECK(data.GetIds()[i] == ids[i]);
        for (unsigned j = 0; j < 3; ++j)
            KRATOS_CHECK(data.GetCoordinates()[i][j] == coords[i][j]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_CreateNodes, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::Vector<int> ids(3);
    HDF5::Vector<array_1d<double, 3>> coords(3);
    for (unsigned i = 0; i < 3; ++i)
    {
        ids[i] = i + 1;
        coords[i] = array_1d<double, 3>(3, 1.2345);
    }
    test_file.WriteDataSet("/Nodes/Ids", ids);
    test_file.WriteDataSet("/Nodes/Coordinates", coords);
    HDF5::Internals::PointsData data;
    data.ReadData(test_file, "/Nodes", 0, 3);
    HDF5::NodesContainerType nodes;
    data.CreateNodes(nodes);
    KRATOS_CHECK(nodes.size() == 3);
    for (unsigned i = 0; i < 3; ++i)
    {
        HDF5::NodeType& r_node = nodes[ids[i]];
        for (unsigned j = 0; j < 3; ++j)
            KRATOS_CHECK(r_node.Coordinates()[j] == coords[i][j]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_SetData, KratosHDF5TestSuite)
{
    HDF5::NodesContainerType nodes;
    for (unsigned i = 0; i < 3; ++i)
    {
        auto p_node = boost::make_shared<HDF5::NodeType>(
            i + 1, 1.2345, 1.2345, 1.2345);
        nodes.push_back(p_node);
    }
    HDF5::Internals::PointsData data;
    data.SetData(nodes);
    KRATOS_CHECK(data.size() == 3);
    for (int i = 0; i < 3; ++i)
    {
        KRATOS_CHECK(data.GetIds()[i] == i + 1);
        for (unsigned j = 0; j < 3; ++j)
            KRATOS_CHECK(data.GetCoordinates()[i][j] == 1.2345);
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_WriteData, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::NodesContainerType nodes;
    for (unsigned i = 0; i < 3; ++i)
    {
        auto p_node = boost::make_shared<HDF5::NodeType>(
            i + 1, 1.2345, 1.2345, 1.2345);
        nodes.push_back(p_node);
    }
    HDF5::Internals::PointsData data;
    data.SetData(nodes);
    data.WriteData(test_file, "/Nodes");
    HDF5::Vector<int> ids(3);
    HDF5::Vector<array_1d<double, 3>> coords(3);
    test_file.ReadDataSet("/Nodes/Ids", ids, 0, 3);
    test_file.ReadDataSet("/Nodes/Coordinates", coords, 0, 3);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            KRATOS_CHECK(coords[i][j] == nodes[ids[i]].Coordinates()[j]);
}

} // namespace Testing
} // namespace Kratos.
