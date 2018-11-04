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

// Application includes
#include "tests/test_utils.h"
#include "custom_io/hdf5_file_serial.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet1, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<int> data_out = TestVector<int>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Vector<int> data_in;
        test_file.ReadDataSet("/data", data_in, 0, data_out.size());
        KRATOS_CHECK(data_in.size() == data_out.size());
        for (std::size_t i = 0; i < data_out.size(); ++i)
            KRATOS_CHECK(data_in[i] == data_out[i]);
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet2, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<double> data_out = TestVector<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Vector<double> data_in;
        test_file.ReadDataSet("/data", data_in, 0, data_out.size());
        KRATOS_CHECK(data_in.size() == data_out.size());
        for (std::size_t i = 0; i < data_out.size(); ++i)
            KRATOS_CHECK(data_in[i] == data_out[i]);
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}


KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet3, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<array_1d<double,3>> data_out = TestVector<array_1d<double,3>>();
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Vector<array_1d<double,3>> data_in;
        test_file.ReadDataSet("/data", data_in, 0, data_out.size());
        KRATOS_CHECK(data_in.size() == data_out.size());
        for (std::size_t i = 0; i < data_out.size(); ++i)
            for (std::size_t j = 0; j < 3; ++j)
                KRATOS_CHECK(data_in(i)(j) == data_out(i)(j));
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet4, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Matrix<double> data_out = TestMatrix<double>(3, 2);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Matrix<double> data_in;
        test_file.ReadDataSet("/data", data_in, 0, data_out.size1());
        KRATOS_CHECK(data_in.size1() == data_out.size1());
        KRATOS_CHECK(data_in.size2() == data_out.size2());
        for (std::size_t i = 0; i < data_out.size1(); ++i)
            for (std::size_t j = 0; j < data_out.size2(); ++j)
                KRATOS_CHECK(data_in(i,j) == data_out(i,j));
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet5, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("invalid_path", data_in, 0, 1);
            , "Invalid path: \"invalid_path\"");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet6, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/not/a/dataset", data_in, 0, 1);
            , "Path is not a data set: /not/a/dataset");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet7, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<double> data_out = TestVector<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Matrix<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/data", data_in, 0, data_out.size());
            , "Invalid data set dimension.");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet8, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<double> data_out = TestVector<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Vector<int> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/data", data_in, 0, data_out.size());
            , "Data type is not int: /data");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet9, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Vector<double> data_out = TestVector<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Vector<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/data", data_in, 10, 3);
            , "StartIndex (10) + BlockSize (3) > size of data set (3).");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet10, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Matrix<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("invalid_path", data_in, 0, 1);
            , "Invalid path: \"invalid_path\"");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet11, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Matrix<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/not/a/dataset", data_in, 0, 1);
            , "Path is not a data set: /not/a/dataset");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet12, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Matrix<double> data_out = TestMatrix<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Vector<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/data", data_in, 0, data_out.size1());
            , "Invalid data set dimension.");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet13, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Matrix<double> data_out = TestMatrix<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Matrix<int> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/data", data_in, 0, data_out.size1());
            , "Data type is not int: /data");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_FileSerial_ReadDataSet14, KratosHDF5TestSuite)
{
    KRATOS_TRY;
    {
        auto test_file = GetTestSerialFile();
        HDF5::File::Matrix<double> data_out = TestMatrix<double>(3);
        HDF5::WriteInfo info;
        test_file.WriteDataSet("/data", data_out, info);
        HDF5::File::Matrix<double> data_in;
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            test_file.ReadDataSet("/data", data_in, 10, 3);
            , "StartIndex (10) + BlockSize (3) > size of data set (3).");
        KRATOS_CHECK(test_file.GetOpenObjectsCount() == 1); // Check for leaks.
    }
    H5close(); // Clean HDF5 for next unit test.
    KRATOS_CATCH_WITH_BLOCK("", H5close(););
}

} // namespace Testing
} // namespace Kratos.
