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
#include <string>

// External includes

// Project includes
#include "testing/testing.h"
#include "custom_io/hdf5_file_serial.h"

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

    HDF5FileSerial test_file(test_params);
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
    HDF5File::Vector<int> write_data_vec_int(5);
    for (unsigned i = 0; i < 5; ++i)
        write_data_vec_int[i] = i + 1;
    test_file.WriteDataSet("/foo/data_vec_int", write_data_vec_int);
    KRATOS_CHECK(test_file.IsDataSet("/foo/bar") == false);
    KRATOS_CHECK(test_file.IsDataSet("/foo/data_vec_int") == true);
    array_1d<double, 3> x1(3, 2.718282);
    array_1d<double, 3> x2(3, 3.141592);
    HDF5File::Vector<array_1d<double, 3>> write_data_array1d(2);
    write_data_array1d[0] = x1;
    write_data_array1d[1] = x2;
    test_file.WriteDataSet("/foo/bar/data_array1d", write_data_array1d);
    KRATOS_CHECK(test_file.IsDataSet("/foo/bar/data_array1d") == true);
    HDF5File::Matrix<int> write_data_mat_int(2, 2);
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j)
            write_data_mat_int(i, j) = 10 * i + j;
    test_file.WriteDataSet("/foo/data_mat_int", write_data_mat_int);

    // Check GetDataDimensions().
    std::vector<unsigned> dims_vec_int = test_file.GetDataDimensions("/foo/data_vec_int");
    KRATOS_CHECK(dims_vec_int.size() == 1);
    KRATOS_CHECK(dims_vec_int[0] == 5);
    std::vector<unsigned> dims_array1d = test_file.GetDataDimensions("/foo/bar/data_array1d");
    KRATOS_CHECK(dims_array1d.size() == 2);
    KRATOS_CHECK(dims_array1d[0] == 2);
    KRATOS_CHECK(dims_array1d[1] == 3);
    std::vector<unsigned> dims_mat_int = test_file.GetDataDimensions("/foo/data_mat_int");
    KRATOS_CHECK(dims_mat_int.size() == 2);
    KRATOS_CHECK(dims_mat_int[0] == 2);
    KRATOS_CHECK(dims_mat_int[1] == 2);

    // Check ReadDataSet() reading full data set.
    HDF5File::Vector<int> read_data_vec_int;
    HDF5File::Vector<array_1d<double, 3>> read_data_array1d;
    HDF5File::Matrix<int> read_data_mat_int;
    test_file.ReadDataSet("/foo/data_vec_int", read_data_vec_int, 0, 5);
    test_file.ReadDataSet("/foo/bar/data_array1d", read_data_array1d, 0, 2);
    test_file.ReadDataSet("/foo/data_mat_int", read_data_mat_int, 0, 2);
    for (unsigned i = 0; i < write_data_vec_int.size(); ++i)
        KRATOS_CHECK(write_data_vec_int[i] == read_data_vec_int[i]);
    for (unsigned i = 0; i < write_data_array1d.size(); ++i)
        for (unsigned j = 0; j < write_data_array1d[i].size(); ++j)
            KRATOS_CHECK(write_data_array1d[i][j] == read_data_array1d[i][j]);
    for (unsigned i = 0; i < write_data_mat_int.size1(); ++i)
        for (unsigned j = 0; j < write_data_mat_int.size2(); ++j)
            KRATOS_CHECK(write_data_mat_int(i,j) == read_data_mat_int(i,j));

    // Check ReadDataSet() reading partial data set.
    read_data_array1d.resize(0);
    test_file.ReadDataSet("/foo/bar/data_array1d", read_data_array1d, 0, 1);
    KRATOS_CHECK(read_data_array1d.size() == 1);
    for (unsigned j = 0; j < write_data_array1d[0].size(); ++j)
        KRATOS_CHECK(write_data_array1d[0][j] == read_data_array1d[0][j]);

    // Check HasIntDataType().
    KRATOS_CHECK(test_file.HasIntDataType("/foo/data_vec_int") == true);
    KRATOS_CHECK(test_file.HasFloatDataType("/foo/data_vec_int") == false);
    KRATOS_CHECK(test_file.HasIntDataType("/foo/bar/data_array1d") == false);
    KRATOS_CHECK(test_file.HasFloatDataType("/foo/bar/data_array1d") == true);
    KRATOS_CHECK(test_file.HasIntDataType("/foo/data_mat_int") == true);

    // Check GetFileName().
    KRATOS_CHECK(test_file.GetFileName() == "test.h5");

    // Check scalar attribute.
    test_file.WriteAttribute("/foo", "DENSITY", 1.2);
    test_file.WriteAttribute("/foo", "INERTIA", write_data_vec_int);
    test_file.WriteAttribute("/foo", "MATERIAL", write_data_mat_int);
    std::vector<std::string> attr_names;
    test_file.GetAttributeNames("/foo", attr_names);
    KRATOS_CHECK(attr_names.size() == 3);
    KRATOS_CHECK(attr_names[0] == "DENSITY");
    KRATOS_CHECK(attr_names[1] == "INERTIA");
    KRATOS_CHECK(attr_names[2] == "MATERIAL");
}
} // namespace Testing
} // namespace Kratos.
