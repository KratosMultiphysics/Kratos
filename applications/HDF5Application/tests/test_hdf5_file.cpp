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
#include <cmath>
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(HDF5Test1, KratosHDF5TestSuite)
{
    Parameters test_params(R"(
                {
                    "file_name" : "test.h5",
                    "file_access_mode": "exclusive",
                    "file_driver": "core"
                })");

    HDF5File test_file(test_params);
    // Check IsPath().
    KRATOS_CHECK(HDF5Utils::IsPath("") == false);
    KRATOS_CHECK(HDF5Utils::IsPath("/") == false);
    KRATOS_CHECK(HDF5Utils::IsPath("/foo") == true);
    KRATOS_CHECK(HDF5Utils::IsPath("/foo/") == false);
    KRATOS_CHECK(HDF5Utils::IsPath("/foo/bar") == true);
    // Check HasPath().
    KRATOS_CHECK(test_file.HasPath("/foo/bar") == false);
    test_file.CreateGroup("/foo");
    test_file.CreateGroup("/foo/bar");
    KRATOS_CHECK(test_file.HasPath("/foo/bar") == true);
    // Check IsGroup().
    KRATOS_CHECK(test_file.IsGroup("/foo") == true);
    KRATOS_CHECK(test_file.IsGroup("/foo/bar") == true);
    KRATOS_CHECK(test_file.IsGroup("/abcdef") == false);
    // Check IsDataSet(), WriteDataSet().
    std::vector<int> write_data_1{1, 2, 3, 4, 5};
    test_file.WriteDataSet("/foo/data1", write_data_1);
    KRATOS_CHECK(test_file.IsDataSet("/foo/bar") == false);
    KRATOS_CHECK(test_file.IsDataSet("/foo/data1") == true);
    array_1d<double, 3> x1(3, 2.718282);
    array_1d<double, 3> x2(3, 3.141592);
    std::vector<array_1d<double, 3>> write_data_2{x1, x2};
    test_file.WriteDataSet("/foo/bar/data2", write_data_2);
    KRATOS_CHECK(test_file.IsDataSet("/foo/bar/data2") == true);
    // Check GetDataDimensions().
    std::vector<unsigned> dims1 = test_file.GetDataDimensions("/foo/data1");
    KRATOS_CHECK(dims1.size() == 1);
    KRATOS_CHECK(dims1[0] == 5);
    std::vector<unsigned> dims2 = test_file.GetDataDimensions("/foo/bar/data2");
    KRATOS_CHECK(dims2.size() == 2);
    KRATOS_CHECK(dims2[0] == 2);
    KRATOS_CHECK(dims2[1] == 3);
    // Check ReadDataSet() reading full data set.
    std::vector<int> read_data_1;
    std::vector<array_1d<double, 3>> read_data_2;
    test_file.ReadDataSet("/foo/data1", read_data_1, 5);
    test_file.ReadDataSet("/foo/bar/data2", read_data_2, 2);
    for (unsigned i = 0; i < write_data_1.size(); ++i)
        KRATOS_CHECK(write_data_1[i] == read_data_1[i]);
    for (unsigned i = 0; i < write_data_2.size(); ++i)
        for (unsigned j = 0; j < write_data_2[i].size(); ++j)
            KRATOS_CHECK(write_data_2[i][j] == read_data_2[i][j]);
    // Check ReadDataSet() reading partial data set.
    read_data_2.resize(0);
    test_file.ReadDataSet("/foo/bar/data2", read_data_2, 1);
    KRATOS_CHECK(read_data_2.size() == 1);
    for (unsigned j = 0; j < write_data_2[0].size(); ++j)
        KRATOS_CHECK(write_data_2[0][j] == read_data_2[0][j]);
    // Check HasIntDataType().
    KRATOS_CHECK(test_file.HasIntDataType("/foo/data1") == true);
    KRATOS_CHECK(test_file.HasFloatDataType("/foo/data1") == false);
    KRATOS_CHECK(test_file.HasIntDataType("/foo/bar/data2") == false);
    KRATOS_CHECK(test_file.HasFloatDataType("/foo/bar/data2") == true);
    // Check GetFileName().
    KRATOS_CHECK(test_file.GetFileName() == "test.h5");
}
} // namespace Testing
} // namespace Kratos.
