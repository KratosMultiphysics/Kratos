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
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_parameters.h"
#include "includes/parallel_environment.h"

// Application includes
#include "tests/test_utils.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_IsPath1, KratosHDF5TestSuite)
{
    KRATOS_EXPECT_TRUE(HDF5::Internals::IsPath("") == false);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_IsPath2, KratosHDF5TestSuite)
{
    KRATOS_EXPECT_TRUE(HDF5::Internals::IsPath("/") == false);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_IsPath3, KratosHDF5TestSuite)
{
    KRATOS_EXPECT_TRUE(HDF5::Internals::IsPath("/foo") == true);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_IsPath4, KratosHDF5TestSuite)
{
    KRATOS_EXPECT_TRUE(HDF5::Internals::IsPath("/foo/") == false);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_IsPath5, KratosHDF5TestSuite)
{
    KRATOS_EXPECT_TRUE(HDF5::Internals::IsPath("/foo/bar") == true);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_IsPath6, KratosHDF5TestSuite)
{
    KRATOS_EXPECT_TRUE(HDF5::Internals::IsPath("/foo//bar") == false);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_File1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        HDF5::File test_file(Testing::GetDefaultDataCommunicator(), test_params);
        , "Invalid file name: PLEASE_SPECIFY_HDF5_FILENAME");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_File2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "abcd"
        })");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        HDF5::File test_file(Testing::GetDefaultDataCommunicator(), test_params);
        , "Unsupported \"file_driver\": abcd");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_File3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "abcd",
            "file_driver" : "core"
        })");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        HDF5::File test_file(Testing::GetDefaultDataCommunicator(), test_params);
        , "Invalid \"file_access_mode\": abcd");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_File4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "read_only"
        })");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        HDF5::File test_file(Testing::GetDefaultDataCommunicator(), test_params);
        , "Invalid HDF5 file: test.h5");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasPath1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GetTestFile().HasPath("invalid_path");
        , "Invalid path: \"invalid_path\"");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasPath2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        KRATOS_EXPECT_TRUE(test_file.HasPath("/foo/bar") == false);
        test_file.AddPath("/foo/bar");
        KRATOS_EXPECT_TRUE(test_file.HasPath("/foo/bar") == true);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsGroup1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GetTestFile().IsGroup("invalid_path");
                                     , "Invalid path: \"invalid_path\". Path should start with \"/\" and should only have characters A-Z, a-z, 0-9, \"/\", and \"_\".");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsGroup2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo/bar");
        KRATOS_EXPECT_TRUE(test_file.IsGroup("/foo/bar") == true);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsGroup3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        KRATOS_EXPECT_TRUE(test_file.IsGroup("/abcdef") == false);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsGroup4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::Vector<int> data = TestVector<int>();
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/foo", data, info);
        KRATOS_EXPECT_TRUE(test_file.IsGroup("/foo") == false);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsDataSet1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GetTestFile().IsDataSet("invalid_path");
                                     , "Invalid path: \"invalid_path\"");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsDataSet2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<int> data = TestVector<int>();
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/foo/data", data, info);
        KRATOS_EXPECT_TRUE(test_file.IsDataSet("/foo/data") == true);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsDataSet3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        KRATOS_EXPECT_TRUE(test_file.IsDataSet("/foo") == false);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_IsDataSet4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo/bar");
        KRATOS_EXPECT_TRUE(test_file.IsDataSet("/foo/bar") == false);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasAttribute1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GetTestFile().HasAttribute("/abcd", "DENSITY");
        , "Error occured in H5Aexists_by_name [ hdf5 error code = -1 ].");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasAttribute2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        KRATOS_EXPECT_TRUE(test_file.HasAttribute("/foo", "DENSITY") == false);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasAttribute3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        test_file.WriteAttribute("/foo", "DENSITY", 1.2);
        KRATOS_EXPECT_TRUE(test_file.HasAttribute("/foo", "DENSITY") == true);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetAttributeNames1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GetTestFile().GetAttributeNames("/invalid/path");, "Error occured in H5Oopen [ hdf5 error code = -1 ].");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetAttributeNames2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        KRATOS_EXPECT_TRUE(test_file.GetAttributeNames("/foo").size() == 0);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetAttributeNames3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        test_file.WriteAttribute("/foo", "DENSITY", 1.2);
        test_file.WriteAttribute("/foo", "VISCOSITY", 1e-5);
        auto names = test_file.GetAttributeNames("/foo");
        KRATOS_EXPECT_TRUE(names.size() == 2);
        KRATOS_EXPECT_TRUE(names[0] == "DENSITY");
        KRATOS_EXPECT_TRUE(names[1] == "VISCOSITY");
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_CreateGroup1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GetTestFile().CreateGroup("/");
                                     , "Error occured in H5Gcreate [ hdf5 error code = -1 ]");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_CreateGroup2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.CreateGroup("/foo");
        KRATOS_EXPECT_TRUE(test_file.IsGroup("/foo"));
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_CreateGroup3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        H5_stderr_muter muter;
        auto test_file = GetTestFile();
        test_file.CreateGroup("/foo");
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(test_file.CreateGroup("/foo");
                                         , "Error occured in H5Gcreate [ hdf5 error code = -1 ]");
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetLinkNames1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GetTestFile().GetLinkNames("/foo");
                                     , "Error occured in H5Gopen [ hdf5 error code = -1 ].");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetLinkNames2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        auto names = test_file.GetLinkNames("/foo");
        KRATOS_EXPECT_TRUE(names.size() == 0);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetLinkNames3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        test_file.AddPath("/foo/group");
        HDF5::File::Vector<double> data = TestVector<double>();
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/foo/data", data, info);
        HDF5::File& r_test_file = test_file;
        auto names = r_test_file.GetLinkNames("/foo");
        KRATOS_EXPECT_TRUE(names.size() == 2);
        KRATOS_EXPECT_TRUE(names[0] == "data");
        KRATOS_EXPECT_TRUE(names[1] == "group");
        KRATOS_EXPECT_TRUE(r_test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetGroupNames1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        test_file.AddPath("/foo/group");
        HDF5::File::Vector<double> data = TestVector<double>();
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/foo/data", data, info);
        HDF5::File& r_test_file = test_file;
        auto names = r_test_file.GetGroupNames("/foo");
        KRATOS_EXPECT_TRUE(names.size() == 1);
        KRATOS_EXPECT_TRUE(names[0] == "group");
        KRATOS_EXPECT_TRUE(r_test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_AddPath1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GetTestFile().AddPath("invalid_path");
                                     , "Invalid path: \"invalid_path\". Path should start with \"/\" and should only have characters A-Z, a-z, 0-9, \"/\", and \"_\".");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_AddPath2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo/bar");
        KRATOS_EXPECT_TRUE(test_file.IsGroup("/foo/bar"));
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_AddPath3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo/bar");
        test_file.AddPath("/foo/bar"); // Redundant call must not throw.
        KRATOS_EXPECT_TRUE(test_file.IsGroup("/foo/bar"));
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_AddPath4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        HDF5::File::Vector<double> data = TestVector<double>();
        HDF5::WriteInfo info;
        auto test_file = GetTestSerialFile();
        test_file.WriteDataSet("/data", data, info);
        HDF5::File& r_test_file = test_file;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            r_test_file.AddPath("/data");
            , "Path exists and is not a group: /data");
        KRATOS_EXPECT_TRUE(r_test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_WriteAttribute1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        const int int_attr = 1;
        const double double_attr = 1.0;
        const std::string string_attr = "value";
        const HDF5::Vector<int> vec_int_attr = TestVector<int>();
        const HDF5::Vector<double> vec_double_attr = TestVector<double>();
        const HDF5::Matrix<int> mat_int_attr = TestMatrix<int>();
        const HDF5::Matrix<double> mat_double_attr = TestMatrix<double>();
        const array_1d<double, 3> array_1d_attr = array_1d<double, 3>(3, 1.234);
        test_file.WriteAttribute("/foo", "INT_ATTRIBUTE", int_attr);
        test_file.WriteAttribute("/foo", "DOUBLE_ATTRIBUTE", double_attr);
        test_file.WriteAttribute("/foo", "STRING_ATTRIBUTE", string_attr);
        test_file.WriteAttribute("/foo", "VEC_INT_ATTRIBUTE", vec_int_attr);
        test_file.WriteAttribute("/foo", "VEC_DOUBLE_ATTRIBUTE", vec_double_attr);
        test_file.WriteAttribute("/foo", "MAT_INT_ATTRIBUTE", mat_int_attr);
        test_file.WriteAttribute("/foo", "MAT_DOUBLE_ATTRIBUTE", mat_double_attr);
        test_file.WriteAttribute("/foo", "ARRAY1D_ATTRIBUTE", array_1d_attr);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_WriteAttribute2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        H5_stderr_muter muter;
        auto test_file = GetTestFile();
        const int attr = 1;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.WriteAttribute("/foo", "ATTRIBUTE", attr);
            , "Error occured in H5Aexists_by_name [ hdf5 error code = -1 ].");
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_WriteAttribute3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        H5_stderr_muter muter;
        auto test_file = GetTestFile();
        const double double_attr = 1.0;
        test_file.CreateGroup("/foo");
        test_file.WriteAttribute("/foo", "ATTRIBUTE", double_attr);
        const HDF5::Matrix<int> mat_int_attr = TestMatrix<int>();
        test_file.WriteAttribute("/foo", "ATTRIBUTE", mat_int_attr);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_DeleteAttribute1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        H5_stderr_muter muter;
        auto test_file = GetTestFile();
        const double attr = 1.0;
        test_file.CreateGroup("/foo");
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr);
        test_file.DeleteAttribute("/foo", "ATTRIBUTE");
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.DeleteAttribute("/foo", "ATTRIBUTE");
            , "Error occured in H5Adelete_by_name [ hdf5 error code = -1 ].");
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetDataDimensions1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<double> data1 = TestVector<double>(5);
        HDF5::File::Vector<array_1d<double,3>> data2 = TestVector<array_1d<double,3>>(3);
        HDF5::File::Matrix<double> data3 = TestMatrix<double>(2, 2);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data1", data1, info);
        test_file.WriteDataSet("/data2", data2, info);
        test_file.WriteDataSet("/data3", data3, info);
        HDF5::File& r_test_file = test_file;
        std::vector<unsigned> dims;
        dims = r_test_file.GetDataDimensions("/data1");
        KRATOS_EXPECT_TRUE(dims.size() == 1);
        KRATOS_EXPECT_TRUE(dims[0] == 5);
        dims = r_test_file.GetDataDimensions("/data2");
        KRATOS_EXPECT_TRUE(dims.size() == 2);
        KRATOS_EXPECT_TRUE(dims[0] == 3);
        KRATOS_EXPECT_TRUE(dims[1] == 3);
        dims = r_test_file.GetDataDimensions("/data3");
        KRATOS_EXPECT_TRUE(dims.size() == 2);
        KRATOS_EXPECT_TRUE(dims[0] == 2);
        KRATOS_EXPECT_TRUE(dims[1] == 2);
        KRATOS_EXPECT_TRUE(r_test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetDataDimensions2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        H5_stderr_muter muter;
        auto test_file = GetTestFile();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(test_file.GetDataDimensions("/data");
                                         , "Error occured in H5Dopen [ hdf5 error code = -1 ].");
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasIntDataType1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<int> data1 = TestVector<int>(3);
        HDF5::File::Vector<double> data2 = TestVector<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data1", data1, info);
        test_file.WriteDataSet("/data2", data2, info);
        HDF5::File& r_test_file = test_file;
        KRATOS_EXPECT_TRUE(r_test_file.HasIntDataType("/data1") == true);
        KRATOS_EXPECT_TRUE(r_test_file.HasIntDataType("/data2") == false);
        KRATOS_EXPECT_TRUE(r_test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasIntDataType2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GetTestFile().HasIntDataType("/invalid/path");
                                     , "Error occured in H5Dopen [ hdf5 error code = -1 ].");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasFloatDataType1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<int> data1 = TestVector<int>(3);
        HDF5::File::Vector<double> data2 = TestVector<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data1", data1, info);
        test_file.WriteDataSet("/data2", data2, info);
        HDF5::File& r_test_file = test_file;
        KRATOS_EXPECT_TRUE(r_test_file.HasFloatDataType("/data1") == false);
        KRATOS_EXPECT_TRUE(r_test_file.HasFloatDataType("/data2") == true);
        KRATOS_EXPECT_TRUE(r_test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_HasFloatDataType2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    H5_stderr_muter muter;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GetTestSerialFile().HasFloatDataType("/invalid/path");
                                     , "Error occured in H5Dopen [ hdf5 error code = -1 ].");
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_GetFileName1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        KRATOS_EXPECT_TRUE(test_file.GetFileName() == "test.h5");
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        const double double_attr_out = 1.0;
        test_file.WriteAttribute("/foo", "DOUBLE_ATTRIBUTE", double_attr_out);
        double double_attr_in;
        test_file.ReadAttribute("/foo", "DOUBLE_ATTRIBUTE", double_attr_in);
        KRATOS_EXPECT_TRUE(double_attr_in == double_attr_out);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        const std::string string_attr_out = "value";
        test_file.WriteAttribute("/foo", "STRING_ATTRIBUTE", string_attr_out);
        std::string string_attr_in;
        test_file.ReadAttribute("/foo", "STRING_ATTRIBUTE", string_attr_in);
        KRATOS_EXPECT_TRUE(string_attr_in == string_attr_out);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        const HDF5::Vector<double> vec_double_attr_out = TestVector<double>();
        test_file.WriteAttribute("/foo", "VEC_DOUBLE_ATTRIBUTE", vec_double_attr_out);
        HDF5::Vector<double> vec_double_attr_in;
        test_file.ReadAttribute("/foo", "VEC_DOUBLE_ATTRIBUTE", vec_double_attr_in);
        KRATOS_EXPECT_TRUE(vec_double_attr_in.size() == vec_double_attr_out.size());
        for (std::size_t i = 0; i < vec_double_attr_in.size(); ++i)
            KRATOS_EXPECT_TRUE(vec_double_attr_in(i) == vec_double_attr_out(i));
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}


KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        const HDF5::Matrix<double> mat_double_attr_out = TestMatrix<double>();
        test_file.WriteAttribute("/foo", "MAT_DOUBLE_ATTRIBUTE", mat_double_attr_out);
        HDF5::Matrix<double> mat_double_attr_in;
        test_file.ReadAttribute("/foo", "MAT_DOUBLE_ATTRIBUTE", mat_double_attr_in);
        KRATOS_EXPECT_TRUE(mat_double_attr_in.size1() == mat_double_attr_out.size1());
        KRATOS_EXPECT_TRUE(mat_double_attr_in.size2() == mat_double_attr_out.size2());
        for (std::size_t i = 0; i < mat_double_attr_in.size1(); ++i)
            for (std::size_t j = 0; j < mat_double_attr_in.size2(); ++j)
                KRATOS_EXPECT_TRUE(mat_double_attr_in(i, j) == mat_double_attr_out(i, j));
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute5, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.AddPath("/foo");
        const array_1d<double, 3> array_1d_attr_out = array_1d<double, 3>(3, 1.234);
        test_file.WriteAttribute("/foo", "ARRAY1D_ATTRIBUTE", array_1d_attr_out);
        array_1d<double, 3> array_1d_attr_in;
        test_file.ReadAttribute("/foo", "ARRAY1D_ATTRIBUTE", array_1d_attr_in);
        KRATOS_EXPECT_TRUE(array_1d_attr_in.size() == array_1d_attr_out.size());
        for (std::size_t i = 0; i < array_1d_attr_in.size(); ++i)
            KRATOS_EXPECT_TRUE(array_1d_attr_in[i] == array_1d_attr_out[i]);
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute6, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        H5_stderr_muter muter;
        auto test_file = GetTestFile();
        double attr;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr);
            , "Error occured in H5Aopen_by_name [ hdf5 error code = -1 ].")
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute7, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        H5_stderr_muter muter;
        auto test_file = GetTestFile();
        double attr;
        test_file.CreateGroup("/foo");
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr);
            , "Error occured in H5Aopen_by_name [ hdf5 error code = -1 ].")
        KRATOS_EXPECT_TRUE(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute8, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        HDF5::File::Vector<double> attr_out = TestVector<double>();
        test_file.CreateGroup("/foo");
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr_out);
        double attr_in;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr_in);
            , "Attribute \"ATTRIBUTE\" has dimension mismatch [ memory dimension = 0, file dimension = 1 ].")
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute9, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.CreateGroup("/foo");
        double attr_out{};
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr_out);
        HDF5::File::Vector<double> attr_in;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr_in);
            , "Attribute \"ATTRIBUTE\" has dimension mismatch [ memory dimension = 1, file dimension = 0 ].")
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute10, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.CreateGroup("/foo");
        double attr_out{};
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr_out);
        HDF5::File::Matrix<double> attr_in;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr_in);
            , "Attribute \"ATTRIBUTE\" has dimension mismatch [ memory dimension = 2, file dimension = 0 ]")
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute11, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        test_file.CreateGroup("/foo");
        double attr_out{};
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr_out);
        std::string attr_in;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr_in);
            , "Memory and file data types are different.")
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute12, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        double attr_out{};
        test_file.CreateGroup("/foo");
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr_out);
        int attr_in{};
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr_in);
            , "Memory and file data types are different.")
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute13, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        HDF5::File::Vector<double> attr_out = TestVector<double>();
        test_file.CreateGroup("/foo");
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr_out);
        HDF5::File::Vector<int> attr_in;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr_in);
            , "Memory and file data types are different.")
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadAttribute14, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestFile();
        HDF5::File::Matrix<double> attr_out = TestMatrix<double>();
        test_file.CreateGroup("/foo");
        test_file.WriteAttribute("/foo", "ATTRIBUTE", attr_out);
        HDF5::File::Matrix<int> attr_in;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            test_file.ReadAttribute("/foo", "ATTRIBUTE", attr_in);
            , "Memory and file data types are different.")
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_File_ReadWriteParameterAttribute1, KratosHDF5TestSuite)
{
    KRATOS_TRY
    {
        Parameters attr_out = Parameters(R"(
        {
            "int"      : 2,
            "double"   : 3.0,
            "string"   : "this is a string",
            "vector"   : [1.0, 2.0, 3.0, 4.0, 5.0],
            "matrix"   : [[1.0, 2.0], [2.0, 3.0]],
            "int_array": [1, 2, 3, 4]
        })" );

        auto test_file = GetTestFile();
        test_file.CreateGroup("/foo");
        test_file.WriteAttribute("/foo", attr_out);

        const auto attr_in = test_file.ReadAttribute("/foo");

        KRATOS_EXPECT_TRUE(attr_in["int"].IsInt() && attr_in["int"].GetInt() == 2);
        KRATOS_EXPECT_TRUE(attr_in["double"].IsDouble() && attr_in["double"].GetDouble() == 3.0);
        KRATOS_EXPECT_TRUE(attr_in["string"].IsString() && attr_in["string"].GetString() == "this is a string");
        array_1d<double, 5> vec_values{1, 2, 3, 4, 5};
        KRATOS_EXPECT_TRUE(attr_in["vector"].IsVector());
        KRATOS_EXPECT_VECTOR_EQ(attr_in["vector"].GetVector(), vec_values);
        Matrix mat_values(2, 2);
        mat_values.data()[0] = 1;
        mat_values.data()[1] = 2;
        mat_values.data()[2] = 2;
        mat_values.data()[3] = 3;
        KRATOS_EXPECT_TRUE(attr_in["matrix"].IsMatrix());
        KRATOS_EXPECT_MATRIX_EQ(attr_in["matrix"].GetMatrix(), mat_values);
        KRATOS_EXPECT_TRUE(attr_in["int_array"].IsArray());
        int local_index = 1;
        for (const auto& v : attr_in["int_array"]) {
            KRATOS_EXPECT_TRUE(v.IsInt() && v.GetInt() == (local_index++));
        }
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

} // namespace Testing
} // namespace Kratos.
