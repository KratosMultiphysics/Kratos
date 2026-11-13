//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#ifdef KRATOS_FUTURE_MDSPAN

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/mdspan.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MdspanConstruction, KratosCoreFastSuite)
{
    // 1D mdspan
    Kratos::Future::extents<int, 5> extents1d;
    std::array<int, 5> data1d = {1, 2, 3, 4, 5};
    Kratos::Future::mdspan<int, Kratos::Future::extents<int, 5>> mds1d(data1d.data(), extents1d);
    KRATOS_EXPECT_EQ(mds1d.rank(), 1);
    KRATOS_EXPECT_EQ(mds1d.extent(0), 5);

    // 2D mdspan with double
    Kratos::Future::extents<int, 2, 3> extents2d_double;
    std::array<double, 6> data2d_double = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 2, 3>> mds2d_double(data2d_double.data(), extents2d_double);
    KRATOS_EXPECT_EQ(mds2d_double.rank(), 2);
    KRATOS_EXPECT_EQ(mds2d_double.extent(0), 2);
    KRATOS_EXPECT_EQ(mds2d_double.extent(1), 3);

    // 2D mdspan with float
    Kratos::Future::extents<int, 3, 2> extents2d_float;
    std::array<float, 6> data2d_float = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    Kratos::Future::mdspan<float, Kratos::Future::extents<int, 3, 2>> mds2d_float(data2d_float.data(), extents2d_float);
    KRATOS_EXPECT_EQ(mds2d_float.rank(), 2);
    KRATOS_EXPECT_EQ(mds2d_float.extent(0), 3);
    KRATOS_EXPECT_EQ(mds2d_float.extent(1), 2);

    // 3D mdspan
    Kratos::Future::extents<int, 2, 3, 4> extents3d;
    std::array<double, 24> data3d; // Initialize if needed, not strictly necessary for this construction test part
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 2, 3, 4>> mds3d(data3d.data(), extents3d);
    KRATOS_EXPECT_EQ(mds3d.rank(), 3);
    KRATOS_EXPECT_EQ(mds3d.extent(0), 2);
    KRATOS_EXPECT_EQ(mds3d.extent(1), 3);
    KRATOS_EXPECT_EQ(mds3d.extent(2), 4);
}

KRATOS_TEST_CASE_IN_SUITE(MdspanDataAccess, KratosCoreFastSuite)
{
    // 1D mdspan - int
    Kratos::Future::extents<int, 5> extents1d_int;
    std::array<int, 5> data1d_int = {1, 2, 3, 4, 5};
    Kratos::Future::mdspan<int, Kratos::Future::extents<int, 5>> mds1d_int(data1d_int.data(), extents1d_int);
    KRATOS_EXPECT_EQ(mds1d_int(0), 1);
    KRATOS_EXPECT_EQ(mds1d_int(4), 5);
    mds1d_int(2) = 10;
    KRATOS_EXPECT_EQ(mds1d_int(2), 10);

    // 2D mdspan - double
    Kratos::Future::extents<int, 2, 3> extents2d_double;
    std::array<double, 6> data2d_double = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 2, 3>> mds2d_double(data2d_double.data(), extents2d_double);
    KRATOS_EXPECT_DOUBLE_EQ(mds2d_double(0, 0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds2d_double(1, 2), 6.0);
    mds2d_double(0, 1) = 7.5;
    KRATOS_EXPECT_DOUBLE_EQ(mds2d_double(0, 1), 7.5);

    // 3D mdspan - double
    Kratos::Future::extents<int, 2, 2, 2> extents3d_double;
    std::array<double, 8> data3d_double = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 2, 2, 2>> mds3d_double(data3d_double.data(), extents3d_double);
    KRATOS_EXPECT_DOUBLE_EQ(mds3d_double(0, 0, 0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds3d_double(1, 1, 1), 8.0);
    mds3d_double(0, 1, 0) = 9.5;
    KRATOS_EXPECT_DOUBLE_EQ(mds3d_double(0, 1, 0), 9.5);
}

KRATOS_TEST_CASE_IN_SUITE(MdspanLayoutLeftAndRightMapping, KratosCoreFastSuite)
{
    // Layout Left
    Kratos::Future::extents<int, 2, 3> extents_left;
    std::array<double, 6> data_left = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // Column-major
    Kratos::Future::layout_left::mapping<Kratos::Future::extents<int, 2, 3>> mapping_left(extents_left);
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 2, 3>, Kratos::Future::layout_left> mds_left(data_left.data(), mapping_left);

    KRATOS_EXPECT_DOUBLE_EQ(mds_left(0, 0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_left(1, 0), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_left(0, 1), 3.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_left(1, 1), 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_left(0, 2), 5.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_left(1, 2), 6.0);

    // Layout Right
    Kratos::Future::extents<int, 2, 3> extents_right;
    std::array<double, 6> data_right = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // Row-major
    Kratos::Future::layout_right::mapping<Kratos::Future::extents<int, 2, 3>> mapping_right(extents_right);
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 2, 3>, Kratos::Future::layout_right> mds_right(data_right.data(), mapping_right);

    KRATOS_EXPECT_DOUBLE_EQ(mds_right(0, 0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_right(0, 1), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_right(0, 2), 3.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_right(1, 0), 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_right(1, 1), 5.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds_right(1, 2), 6.0);
}

KRATOS_TEST_CASE_IN_SUITE(MdspanLayoutStride, KratosCoreFastSuite)
{
    Kratos::Future::extents<int, 2, 3> extents;
    std::array<double, 12> data = { // Data array is larger to accommodate custom strides
        1.0, 0.0, 2.0, 0.0, 3.0, 0.0, // Row 0 with stride
        4.0, 0.0, 5.0, 0.0, 6.0, 0.0  // Row 1 with stride
    };
    std::array<int, 2> strides = {6, 2}; // Stride for first dim, stride for second dim
    Kratos::Future::layout_stride::mapping<Kratos::Future::extents<int, 2, 3>> mapping(extents, strides);
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 2, 3>, Kratos::Future::layout_stride> mds(data.data(), mapping);

    KRATOS_EXPECT_DOUBLE_EQ(mds(0, 0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds(0, 1), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds(0, 2), 3.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds(1, 0), 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds(1, 1), 5.0);
    KRATOS_EXPECT_DOUBLE_EQ(mds(1, 2), 6.0);

    // Test with non-contiguous elements in the original data due to stride
    KRATOS_EXPECT_DOUBLE_EQ(data[mapping(0,0)], 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(data[mapping(0,1)], 2.0); // data[2]
    KRATOS_EXPECT_DOUBLE_EQ(data[mapping(0,2)], 3.0); // data[4]
    KRATOS_EXPECT_DOUBLE_EQ(data[mapping(1,0)], 4.0); // data[6]
    KRATOS_EXPECT_DOUBLE_EQ(data[mapping(1,1)], 5.0); // data[8]
    KRATOS_EXPECT_DOUBLE_EQ(data[mapping(1,2)], 6.0); // data[10]
}

KRATOS_TEST_CASE_IN_SUITE(MdspanViewSlice, KratosCoreFastSuite)
{
    Kratos::Future::extents<int, 3, 4> extents;
    std::array<double, 12> data = {
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0
    };
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 3, 4>> mds(data.data(), extents);

    // Create a slice: view of the second row (index 1)
    // auto slice = Kratos::Future::mdspan_view(mds, 1, Kratos::Future::full_extent); // Alias template deduction only available with ‘-std=c++20’ or ‘-std=gnu++20’
    // In C++17, we explicitly state the resulting type instead of using 'auto'.
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 4>> slice = Kratos::Future::submdspan(mds, 1, Kratos::Future::full_extent);

    KRATOS_EXPECT_EQ(slice.rank(), 1);
    KRATOS_EXPECT_EQ(slice.extent(0), 4);
    KRATOS_EXPECT_DOUBLE_EQ(slice(0), 5.0);
    KRATOS_EXPECT_DOUBLE_EQ(slice(1), 6.0);
    KRATOS_EXPECT_DOUBLE_EQ(slice(2), 7.0);
    KRATOS_EXPECT_DOUBLE_EQ(slice(3), 8.0);

    // Modify data through slice
    slice(1) = 60.0;
    KRATOS_EXPECT_DOUBLE_EQ(mds(1,1), 60.0);

    // Create a slice: view of the second column (index 1)
    // auto col_slice = Kratos::Future::mdspan_view(mds, Kratos::Future::full_extent, 1); // Alias template deduction only available with ‘-std=c++20’ or ‘-std=gnu++20’
    // In C++17, we explicitly state the resulting type instead of using 'auto'.
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 3>> col_slice = Kratos::Future::submdspan(mds, Kratos::Future::full_extent, 1);
    KRATOS_EXPECT_EQ(col_slice.rank(), 1);
    KRATOS_EXPECT_EQ(col_slice.extent(0), 3);
    KRATOS_EXPECT_DOUBLE_EQ(col_slice(0), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(col_slice(1), 60.0); // Modified value
    KRATOS_EXPECT_DOUBLE_EQ(col_slice(2), 10.0);
}

KRATOS_TEST_CASE_IN_SUITE(MdspanViewSubmdspan, KratosCoreFastSuite)
{
    Kratos::Future::extents<int, 3, 4> extents;
    std::array<double, 12> data = {
        1.0,  2.0,  3.0,  4.0,
        5.0,  6.0,  7.0,  8.0,
        9.0, 10.0, 11.0, 12.0
    };
    Kratos::Future::mdspan<double, Kratos::Future::extents<int, 3, 4>> mds(data.data(), extents);

    // Create a submdspan: e.g., a 2x2 sub-array starting from mds(1,1)
    // This would be rows 1-2 and columns 1-2 of the original mdspan
    auto sub = Kratos::Future::submdspan(mds, std::pair<int, int>{1, 3}, std::pair<int, int>{1, 3});

    KRATOS_EXPECT_EQ(sub.rank(), 2);
    KRATOS_EXPECT_EQ(sub.extent(0), 2); // 3-1
    KRATOS_EXPECT_EQ(sub.extent(1), 2); // 3-1

    KRATOS_EXPECT_DOUBLE_EQ(sub(0,0), 6.0); // Corresponds to mds(1,1)
    KRATOS_EXPECT_DOUBLE_EQ(sub(0,1), 7.0); // Corresponds to mds(1,2)
    KRATOS_EXPECT_DOUBLE_EQ(sub(1,0), 10.0); // Corresponds to mds(2,1)
    KRATOS_EXPECT_DOUBLE_EQ(sub(1,1), 11.0); // Corresponds to mds(2,2)

    // Modify data through submdspan
    sub(0,0) = 600.0;
    KRATOS_EXPECT_DOUBLE_EQ(mds(1,1), 600.0);

    // Another submdspan example: first row, first two elements
    auto sub2 = submdspan(mds, 0, std::pair<int,int>{0,2});
    KRATOS_EXPECT_EQ(sub2.rank(), 1);
    KRATOS_EXPECT_EQ(sub2.extent(0), 2);
    KRATOS_EXPECT_DOUBLE_EQ(sub2(0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(sub2(1), 2.0);
}

KRATOS_TEST_CASE_IN_SUITE(Mdspan1DViewAlias, KratosCoreFastSuite)
{
    // Test with int
    std::array<int, 5> data_int = {1, 2, 3, 4, 5};
    Kratos::Future::mdspan_1d_view<int> view_int(data_int.data(), data_int.size());

    KRATOS_EXPECT_EQ(view_int.rank(), 1);
    KRATOS_EXPECT_EQ(view_int.extent(0), 5);
    KRATOS_EXPECT_EQ(view_int(0), 1);
    KRATOS_EXPECT_EQ(view_int(2), 3);
    KRATOS_EXPECT_EQ(view_int(4), 5);

    view_int(1) = 20;
    KRATOS_EXPECT_EQ(data_int[1], 20);
    KRATOS_EXPECT_EQ(view_int(1), 20);

    // Test with double
    std::array<double, 3> data_double = {1.1, 2.2, 3.3};
    Kratos::Future::mdspan_1d_view<double> view_double(data_double.data(), data_double.size());

    KRATOS_EXPECT_EQ(view_double.rank(), 1);
    KRATOS_EXPECT_EQ(view_double.extent(0), 3);
    KRATOS_EXPECT_DOUBLE_EQ(view_double(0), 1.1);
    KRATOS_EXPECT_DOUBLE_EQ(view_double(1), 2.2);

    view_double(2) = 33.3;
    KRATOS_EXPECT_DOUBLE_EQ(data_double[2], 33.3);
    KRATOS_EXPECT_DOUBLE_EQ(view_double(2), 33.3);
}

KRATOS_TEST_CASE_IN_SUITE(Mdspan2DViewAlias, KratosCoreFastSuite)
{
    // Test with int, layout_right (default)
    std::array<int, 6> data_int_lr = {1, 2, 3, 4, 5, 6}; // 2x3 row-major
    Kratos::Future::mdspan_2d_view<int> view_int_lr(data_int_lr.data(), 2, 3);

    KRATOS_EXPECT_EQ(view_int_lr.rank(), 2);
    KRATOS_EXPECT_EQ(view_int_lr.extent(0), 2);
    KRATOS_EXPECT_EQ(view_int_lr.extent(1), 3);
    KRATOS_EXPECT_EQ(view_int_lr(0,0), 1);
    KRATOS_EXPECT_EQ(view_int_lr(0,2), 3);
    KRATOS_EXPECT_EQ(view_int_lr(1,0), 4);
    KRATOS_EXPECT_EQ(view_int_lr(1,2), 6);
    view_int_lr(0,1) = 20;
    KRATOS_EXPECT_EQ(data_int_lr[1], 20);
    KRATOS_EXPECT_EQ(view_int_lr(0,1), 20);

    // Test with double, layout_left
    std::array<double, 6> data_double_ll = {1.1, 4.4, 2.2, 5.5, 3.3, 6.6}; // 2x3 column-major
    Kratos::Future::mdspan_2d_view<double> view_double_ll(data_double_ll.data(), 2, 3);

    KRATOS_EXPECT_EQ(view_double_ll.rank(), 2);
    KRATOS_EXPECT_EQ(view_double_ll.extent(0), 2);
    KRATOS_EXPECT_EQ(view_double_ll.extent(1), 3);
    KRATOS_EXPECT_DOUBLE_EQ(view_double_ll(0,0), 1.1);
    KRATOS_EXPECT_DOUBLE_EQ(view_double_ll(1,0), 5.5);
    KRATOS_EXPECT_DOUBLE_EQ(view_double_ll(0,2), 2.2);
    KRATOS_EXPECT_DOUBLE_EQ(view_double_ll(1,2), 6.6);
    view_double_ll(0,1) = 22.2;
    KRATOS_EXPECT_DOUBLE_EQ(view_double_ll(0,1), 22.2);
}

} // namespace Kratos::Testing

#endif // KRATOS_FUTURE_MDSPAN