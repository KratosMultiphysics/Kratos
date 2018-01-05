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
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_file_serial.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5Utils_IsPath, KratosHDF5TestSuite)
{
    KRATOS_CHECK(HDF5::Internals::IsPath("") == false);
    KRATOS_CHECK(HDF5::Internals::IsPath("/") == false);
    KRATOS_CHECK(HDF5::Internals::IsPath("/foo") == true);
    KRATOS_CHECK(HDF5::Internals::IsPath("/foo/") == false);
    KRATOS_CHECK(HDF5::Internals::IsPath("/foo/bar") == true);
    KRATOS_CHECK(HDF5::Internals::IsPath("/foo//bar") == false);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5Utils_Split, KratosHDF5TestSuite)
{
    std::vector<std::string> result = HDF5::Internals::Split("/foo//bar", '/');
    KRATOS_CHECK(result.size() == 2);
    KRATOS_CHECK(result[0] == "foo");
    KRATOS_CHECK(result[1] == "bar");
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_HDF5File1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        Parameters test_params(R"(
            {
                "file_name" : "test.h5",
                "file_access_mode": "read_only",
                "file_driver": "sec2"
            })");
        HDF5::File test_file(test_params);
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
        , "Invalid HDF5 file: test.h5");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_HDF5File2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        Parameters test_params(R"(
            {
                "file_name" : "test.h5",
                "file_access_mode": "bad_access_mode",
                "file_driver": "core"
            })");
        HDF5::File test_file(test_params);
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
        , "Invalid \"file_access_mode\": bad_access_mode");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_HDF5File3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        Parameters test_params(R"(
            {
                "file_name" : "test.h5",
                "file_access_mode": "exclusive",
                "file_driver": "bad_file_driver"
            })");
        HDF5::File test_file(test_params);
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
        , "Unsupported \"file_driver\": bad_file_driver");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_HasPath, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.HasPath("invalid_path");
        , "Invalid path: \"invalid_path\"");
    KRATOS_CHECK(test_file.HasPath("/foo/bar") == false);
    test_file.CreateGroup("/foo");
    test_file.CreateGroup("/foo/bar");
    KRATOS_CHECK(test_file.HasPath("/foo/bar") == true);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_IsGroup, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);
    test_file.AddPath("/foo/bar");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.IsGroup("invalid_path");
        , "Invalid path: \"invalid_path\"");
    KRATOS_CHECK(test_file.IsGroup("/foo/bar") == true);
    KRATOS_CHECK(test_file.IsGroup("/abcdef") == false);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_IsDataSet, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<int> data(5, 0);
    test_file.WriteDataSet("/foo/data", data);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.IsGroup("invalid_path");
        , "Invalid path: \"invalid_path\"");
    KRATOS_CHECK(test_file.IsDataSet("/foo") == false);
    KRATOS_CHECK(test_file.IsDataSet("/foo/data") == true);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_HasAttribute, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.HasAttribute("/invalid/path", "DENSITY");
        , "H5Aexists_by_name failed");
    test_file.AddPath("/foo");
    KRATOS_CHECK(test_file.HasAttribute("/foo", "DENSITY") == false);
    test_file.WriteAttribute("/foo", "DENSITY", 1.2);
    KRATOS_CHECK(test_file.HasAttribute("/foo", "DENSITY") == true);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_GetAttributeNames, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);
    std::vector<std::string> names;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.GetAttributeNames("/invalid/path", names);
        , "H5Oopen failed");
    test_file.AddPath("/foo");
    test_file.GetAttributeNames("/foo", names);
    KRATOS_CHECK(names.size() == 0);
    test_file.WriteAttribute("/foo", "DENSITY", 1.2);
    test_file.WriteAttribute("/foo", "VISCOSITY", 1e-5);
    test_file.GetAttributeNames("/foo", names);
    KRATOS_CHECK(names.size() == 2);
    KRATOS_CHECK(names[0] == "DENSITY");
    KRATOS_CHECK(names[1] == "VISCOSITY");
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_CreateGroup, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.CreateGroup("/");
        , "H5Gcreate failed");
    test_file.CreateGroup("/foo");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.CreateGroup("/foo");
        , "H5Gcreate failed");
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_AddPath, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.AddPath("invalid_path");
        , "Invalid path: invalid_path");
    HDF5::File::Vector<double> data(5, 1.23);
    test_file.WriteDataSet("/a_data_set", data);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.AddPath("/a_data_set");
        , "Path exists and is not a group: /a_data_set");
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_GetDataDimensions, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<double> data1(5);
    HDF5::File::Vector<array_1d<double,3>> data2(3);
    HDF5::File::Matrix<double> data3(2,2);
    test_file.WriteDataSet("/data1", data1);
    test_file.WriteDataSet("/data2", data2);
    test_file.WriteDataSet("/data3", data3);
    
    std::vector<unsigned> dims;
    dims = test_file.GetDataDimensions("/data1");
    KRATOS_CHECK(dims.size() == 1);
    KRATOS_CHECK(dims[0] == 5);
    dims = test_file.GetDataDimensions("/data2");
    KRATOS_CHECK(dims.size() == 2);
    KRATOS_CHECK(dims[0] == 3);
    KRATOS_CHECK(dims[1] == 3);
    dims = test_file.GetDataDimensions("/data3");
    KRATOS_CHECK(dims.size() == 2);
    KRATOS_CHECK(dims[0] == 2);
    KRATOS_CHECK(dims[1] == 2);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_HasIntDataType, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<int> data1(3);
    HDF5::File::Vector<double> data2(3);
    test_file.WriteDataSet("/data1", data1);
    test_file.WriteDataSet("/data2", data2);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.HasIntDataType("/invalid/path");
        , "H5Dopen failed");
    KRATOS_CHECK(test_file.HasIntDataType("/data1") == true);
    KRATOS_CHECK(test_file.HasIntDataType("/data2") == false);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_HasFloatDataType, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<int> data1(3);
    HDF5::File::Vector<double> data2(3);
    test_file.WriteDataSet("/data1", data1);
    test_file.WriteDataSet("/data2", data2);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.HasIntDataType("/invalid/path");
        , "H5Dopen failed");
    KRATOS_CHECK(test_file.HasFloatDataType("/data1") == false);
    KRATOS_CHECK(test_file.HasFloatDataType("/data2") == true);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_GetFileName, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);
    KRATOS_CHECK(test_file.GetFileName() == "test.h5");
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadDataSet1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<double> data_in, data_out(3, 0.0);
    test_file.WriteDataSet("/data", data_out);
    HDF5::File::Matrix<double> bad_data_container;
    HDF5::File::Vector<int> bad_data_type;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("invalid_path", data_in, 0, 3);
        , "Invalid path: \"invalid_path\"");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/not/a/dataset", data_in, 0, 3);
        , "Path is not a data set: /not/a/dataset");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", bad_data_container, 0, 3);
        , "Invalid data set dimension.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", data_in, 10, 3);
        , "StartIndex (10) + BlockSize (3) > size of data set (3).");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", bad_data_type, 0, 3);
        , "Data type is not int: /data")
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadDataSet2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<array_1d<double,3>> data_in, data_out(3);
    test_file.WriteDataSet("/data", data_out);
    HDF5::File::Vector<double> bad_data_container;
    HDF5::File::Vector<array_1d<int,3>> bad_data_type;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("invalid_path", data_in, 0, 3);
        , "Invalid path: \"invalid_path\"");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/not/a/dataset", data_in, 0, 3);
        , "Path is not a data set: /not/a/dataset");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", bad_data_container, 0, 3);
        , "Invalid data set dimension.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", data_in, 10, 3);
        , "StartIndex (10) + BlockSize (3) > size of data set (3).");
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadDataSet3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Matrix<double> data_in, data_out(3,3);
    test_file.WriteDataSet("/data", data_out);
    HDF5::File::Vector<double> bad_data_container;
    HDF5::File::Matrix<int> bad_data_type;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("invalid_path", data_in, 0, 3);
        , "Invalid path: \"invalid_path\"");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/not/a/dataset", data_in, 0, 3);
        , "Path is not a data set: /not/a/dataset");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", bad_data_container, 0, 3);
        , "Invalid data set dimension.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", data_in, 10, 3);
        , "StartIndex (10) + BlockSize (3) > size of data set (3).");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadDataSet("/data", bad_data_type, 0, 3);
        , "Data type is not int: /data")
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadDataSet4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<int> data_in, data_out(3);
    for (int i = 0; i < 3; ++i)
        data_out[i] = i;
    test_file.WriteDataSet("/data", data_out);
    test_file.ReadDataSet("/data", data_in, 0, 3);
    KRATOS_CHECK(data_in.size() == 3);
    for (int i = 0; i < 3; ++i)
        KRATOS_CHECK(data_in[i] == data_out[i]);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadDataSet5, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<double> data_in, data_out(3);
    for (int i = 0; i < 3; ++i)
        data_out[i] = i;
    test_file.WriteDataSet("/data", data_out);
    test_file.ReadDataSet("/data", data_in, 0, 3);
    KRATOS_CHECK(data_in.size() == 3);
    for (int i = 0; i < 3; ++i)
        KRATOS_CHECK(data_in[i] == data_out[i]);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadDataSet6, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Vector<array_1d<double,3>> data_in, data_out(3);
    for (int i = 0; i < 3; ++i)
        data_out[i] = array_1d<double, 3>(3, 2.718282);
    test_file.WriteDataSet("/data", data_out);
    test_file.ReadDataSet("/data", data_in, 0, 3);
    KRATOS_CHECK(data_in.size() == 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            KRATOS_CHECK(data_in[i][j] == data_out[i][j]);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadDataSet7, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::FileSerial test_file(test_params);
    HDF5::File::Matrix<double> data_in, data_out(3,2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            data_out(i,j) = 2.718282;
    test_file.WriteDataSet("/data", data_out);
    test_file.ReadDataSet("/data", data_in, 0, 3);
    KRATOS_CHECK(data_in.size1() == 3);
    KRATOS_CHECK(data_in.size2() == 2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            KRATOS_CHECK(data_in(i,j) == data_out(i,j));
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    const double property_out = 1.2;
    double property_in;
    
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.WriteAttribute("/foo", "PROPERTY", property_out);
        , "H5Acreate_by_name failed.")
    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY", property_out);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.WriteAttribute("/foo", "PROPERTY", property_out);
        , "H5Acreate_by_name failed.")
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadAttribute("/foo", "PROPERTYY", property_in);
        , "H5Aopen_by_name failed.")
    test_file.ReadAttribute("/foo", "PROPERTY", property_in);
    KRATOS_CHECK(property_in == property_out);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    const double property_out = 1.2;
    int property_bad_type;
    
    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY", property_out);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadAttribute("/foo", "PROPERTY", property_bad_type);
        , "Memory and file data types are different.")
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    double property_in;
    HDF5::File::Vector<double> property_bad_container(2, 1.2);
    
    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY_VEC", property_bad_container);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
       test_file.ReadAttribute("/foo", "PROPERTY_VEC", property_in);
       , "Attribute \"PROPERTY_VEC\" is not scalar.")
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    HDF5::File::Vector<double> property_out(2);
    property_out[0] = property_out[1] = 1.2345;
    HDF5::File::Vector<double> property_in;

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.WriteAttribute("/foo", "PROPERTY", property_out);
        , "H5Acreate_by_name failed.")
    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY", property_out);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.WriteAttribute("/foo", "PROPERTY", property_out);
        , "H5Acreate_by_name failed.")
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadAttribute("/foo", "PROPERTYY", property_in);
        , "H5Aopen_by_name failed.")
    test_file.ReadAttribute("/foo", "PROPERTY", property_in);
    KRATOS_CHECK(property_in.size() == property_out.size());
    KRATOS_CHECK(property_in[0] == property_out[0]);
    KRATOS_CHECK(property_in[1] == property_out[1]);
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute5, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    HDF5::File::Vector<double> property_out(2);
    property_out[0] = property_out[1] = 1.2345;
    HDF5::File::Vector<int> property_bad_type;

    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY", property_out);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadAttribute("/foo", "PROPERTY", property_bad_type);
        , "Memory and file data types are different.")
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute6, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    HDF5::File::Vector<double> property_in;
    double property_bad_container = 1.2345;

    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY_SCAL", property_bad_container);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
       test_file.ReadAttribute("/foo", "PROPERTY_SCAL", property_in);
       , "Attribute \"PROPERTY_SCAL\" is not vector.")
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute7, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    HDF5::File::Matrix<double> property_out(2,2);
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j)
            property_out(i,j) = 1.2345;
    HDF5::File::Matrix<double> property_in;
    
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.WriteAttribute("/foo", "PROPERTY", property_out);
        , "H5Acreate_by_name failed.")
    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY", property_out);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.WriteAttribute("/foo", "PROPERTY", property_out);
        , "H5Acreate_by_name failed.")
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadAttribute("/foo", "PROPERTYY", property_in);
        , "H5Aopen_by_name failed.")
    test_file.ReadAttribute("/foo", "PROPERTY", property_in);
    KRATOS_CHECK(property_in.size1() == property_out.size1());
    KRATOS_CHECK(property_in.size2() == property_out.size2());
    for (unsigned i = 0; i < property_out.size1(); ++i)
        for (unsigned j = 0; j < property_out.size2(); ++j)
            KRATOS_CHECK(property_in(i,j) == property_out(i,j));
    KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute8, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    HDF5::File::Matrix<double> property_out(2,2);
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j)
            property_out(i,j) = 1.2345;
    HDF5::File::Matrix<int> property_bad_type;

    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY", property_out);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadAttribute("/foo", "PROPERTY", property_bad_type);
        , "Memory and file data types are different.")
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5File_ReadAttribute9, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    HDF5::File test_file(test_params);

    HDF5::File::Matrix<double> property_in;
    double property_bad_container = 1.2345;
    
    test_file.CreateGroup("/foo");
    test_file.WriteAttribute("/foo", "PROPERTY_SCAL", property_bad_container);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        test_file.ReadAttribute("/foo", "PROPERTY_SCAL", property_in);
        , "Attribute \"PROPERTY_SCAL\" is not matrix.")
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

} // namespace Testing
} // namespace Kratos.
