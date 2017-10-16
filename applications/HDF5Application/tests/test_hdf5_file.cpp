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
                    "file_driver": "sec2"
                })");
            
            HDF5File test_file(test_params);
            // Check IsPath().
            KRATOS_CHECK(test_file.IsPath("") == false);
            KRATOS_CHECK(test_file.IsPath("/") == false);
            KRATOS_CHECK(test_file.IsPath("/foo") == true);
            KRATOS_CHECK(test_file.IsPath("/foo/") == false);
            KRATOS_CHECK(test_file.IsPath("/foo/bar") == true);
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
            std::vector<int> wd1{1,2,3,4,5};
            test_file.WriteDataSet("/foo/d1", wd1);
            KRATOS_CHECK(test_file.IsDataSet("/foo/bar") == false);
            KRATOS_CHECK(test_file.IsDataSet("/foo/d1") == true);
            array_1d<double, 3> x1(3,2.718282);
            array_1d<double, 3> x2(3,3.141592);
            std::vector<array_1d<double,3>> wd2{x1, x2};
            test_file.WriteDataSet("/foo/bar/d2", wd2);
            KRATOS_CHECK(test_file.IsDataSet("/foo/bar/d2") == true);
            // Check GetDataDimensions(), ReadDataSet().
            std::vector<unsigned> dims1 = test_file.GetDataDimensions("/foo/d1");
            KRATOS_CHECK(dims1.size() == 1);
            KRATOS_CHECK(dims1[0] == 5);
            std::vector<unsigned> dims2 = test_file.GetDataDimensions("/foo/bar/d2");
            KRATOS_CHECK(dims2.size() == 2);
            KRATOS_CHECK(dims2[0] == 2);
            KRATOS_CHECK(dims2[1] == 3);
            std::vector<int> rd1;
            std::vector<array_1d<double,3>> rd2;
            test_file.ReadDataSet("/foo/d1", rd1, 5);
            test_file.ReadDataSet("/foo/bar/d2", rd2, 2);
            //for (unsigned i=0; i < wd1.size(); ++i)
            //    KRATOS_CHECK(wd1[i] == rd1[i]);
            //for (unsigned i=0; i < wd2.size(); ++i)
            //    for (unsigned j=0; j < wd2[i].size(); ++j)
            //        KRATOS_CHECK(wd2[i][j] == rd2[i][j]);
        }
    } // namespace Testing
}  // namespace Kratos.

