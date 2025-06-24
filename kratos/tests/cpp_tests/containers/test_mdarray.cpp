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
#include <vector>
#include <array>

// External includes

// Project includes
#include "containers/mdarray.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MdarrayConstruction, KratosCoreFastSuite)
{
    // Default construction for mdarray_dynamic (if extents are all dynamic, needs runtime sizes)
    Kratos::Future::mdarray_dynamic<int, 2> arr_default_dyn; // Default constructor for mdarray itself is tricky without extents
                                                            // The mdarray_dynamic alias implies dextents, which IS default constructible.
                                                            // However, the mdarray class itself expects extents.
                                                            // Let's test construction with explicit extents.

    // Construction with Kratos::Future::extents
    Kratos::Future::mdarray<int, Kratos::Future::extents<std::size_t, 3, 4>> arr_static_ext;
    KRATOS_EXPECT_EQ(arr_static_ext.rank(), 2);
    KRATOS_EXPECT_EQ(arr_static_ext.static_extent(0), 3);
    KRATOS_EXPECT_EQ(arr_static_ext.static_extent(1), 4);
    KRATOS_EXPECT_EQ(arr_static_ext.extent(0), 3);
    KRATOS_EXPECT_EQ(arr_static_ext.extent(1), 4);
    KRATOS_EXPECT_EQ(arr_static_ext.container().size(), 3 * 4);

    // Construction with Kratos::Future::dextents (dynamic extents)
    Kratos::Future::mdarray<double, Kratos::Future::dextents<std::size_t, 2>> arr_dyn_ext(2, 5);
    KRATOS_EXPECT_EQ(arr_dyn_ext.rank(), 2);
    KRATOS_EXPECT_EQ(arr_dyn_ext.extent(0), 2);
    KRATOS_EXPECT_EQ(arr_dyn_ext.extent(1), 5);
    KRATOS_EXPECT_EQ(arr_dyn_ext.container().size(), 2 * 5);

    // Construction of mdarray_1d and mdarray_2d aliases
    Kratos::Future::mdarray_1d<float> arr_1d_alias(6);
    KRATOS_EXPECT_EQ(arr_1d_alias.rank(), 1);
    KRATOS_EXPECT_EQ(arr_1d_alias.extent(0), 6);
    KRATOS_EXPECT_EQ(arr_1d_alias.container().size(), 6);

    Kratos::Future::mdarray_2d<int> arr_2d_alias(3, 2);
    KRATOS_EXPECT_EQ(arr_2d_alias.rank(), 2);
    KRATOS_EXPECT_EQ(arr_2d_alias.extent(0), 3);
    KRATOS_EXPECT_EQ(arr_2d_alias.extent(1), 2);
    KRATOS_EXPECT_EQ(arr_2d_alias.container().size(), 3 * 2);

    // Copy construction
    Kratos::Future::mdarray_2d<int> arr_2d_copy_src(2, 2);
    arr_2d_copy_src(0,0) = 1; arr_2d_copy_src(0,1) = 2;
    arr_2d_copy_src(1,0) = 3; arr_2d_copy_src(1,1) = 4;

    Kratos::Future::mdarray_2d<int> arr_2d_copy_dest(arr_2d_copy_src);
    KRATOS_EXPECT_EQ(arr_2d_copy_dest.extent(0), 2);
    KRATOS_EXPECT_EQ(arr_2d_copy_dest.extent(1), 2);
    KRATOS_EXPECT_EQ(arr_2d_copy_dest(0,0), 1);
    KRATOS_EXPECT_EQ(arr_2d_copy_dest(1,1), 4);
    KRATOS_EXPECT_NE(arr_2d_copy_dest.data(), arr_2d_copy_src.data()); // Should be a deep copy

    // Move construction
    Kratos::Future::mdarray_2d<int> arr_2d_move_dest(std::move(arr_2d_copy_src));
    KRATOS_EXPECT_EQ(arr_2d_move_dest.extent(0), 2);
    KRATOS_EXPECT_EQ(arr_2d_move_dest.extent(1), 2);
    KRATOS_EXPECT_EQ(arr_2d_move_dest(0,0), 1);
    // arr_2d_copy_src is now in a valid but unspecified state. Accessing its data might be problematic.
    // KRATOS_EXPECT_EQ(arr_2d_copy_src.container().size(), 0); // A common outcome for moved-from std::vector

    // Construction with a specific container (std::array)
    // Note: mdarray with std::array container requires extents to be fully static.
    using StaticExtents = Kratos::Future::extents<std::size_t, 2, 3>;
    using ArrayContainer = std::array<double, 2*3>;
    Kratos::Future::mdarray<double, StaticExtents, Kratos::Future::layout_right, ArrayContainer> arr_std_array;
    KRATOS_EXPECT_EQ(arr_std_array.rank(), 2);
    KRATOS_EXPECT_EQ(arr_std_array.extent(0), 2);
    KRATOS_EXPECT_EQ(arr_std_array.extent(1), 3);
    arr_std_array(0,0) = 1.0; // Should compile and work
}

KRATOS_TEST_CASE_IN_SUITE(MdarrayDataAccess, KratosCoreFastSuite)
{
    // 1D mdarray - int
    Kratos::Future::mdarray_1d<int> arr1d(5);
    for(int i=0; i<5; ++i) arr1d(i) = i + 1;

    KRATOS_EXPECT_EQ(arr1d(0), 1);
    KRATOS_EXPECT_EQ(arr1d(4), 5);
    arr1d(2) = 10;
    KRATOS_EXPECT_EQ(arr1d(2), 10);

    // // Test operator[] for 1D. Only available for mdarray_1d and MDSPAN_USE_BRACKET_OPERATOR
    // KRATOS_EXPECT_EQ(arr1d[0], 1);
    // arr1d[3] = 13;
    // KRATOS_EXPECT_EQ(arr1d[3], 13);
    // KRATOS_EXPECT_EQ(arr1d(3), 13)

    // 2D mdarray - double
    Kratos::Future::mdarray_2d<double> arr2d(2, 3);
    arr2d(0,0) = 1.0; arr2d(0,1) = 2.0; arr2d(0,2) = 3.0;
    arr2d(1,0) = 4.0; arr2d(1,1) = 5.0; arr2d(1,2) = 6.0;

    KRATOS_EXPECT_DOUBLE_EQ(arr2d(0, 0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(arr2d(1, 2), 6.0);
    arr2d(0, 1) = 7.5;
    KRATOS_EXPECT_DOUBLE_EQ(arr2d(0, 1), 7.5);

    // 3D mdarray - double
    Kratos::Future::mdarray_dynamic<double, 3> arr3d(2, 2, 2);
    double val = 1.0;
    for(std::size_t i=0; i<2; ++i)
        for(std::size_t j=0; j<2; ++j)
            for(std::size_t k=0; k<2; ++k)
                arr3d(i,j,k) = val++;

    KRATOS_EXPECT_DOUBLE_EQ(arr3d(0,0,0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(arr3d(1,1,1), 8.0);
    arr3d(0,1,0) = 9.5;
    KRATOS_EXPECT_DOUBLE_EQ(arr3d(0,1,0), 9.5);
}

KRATOS_TEST_CASE_IN_SUITE(MdarrayLayout, KratosCoreFastSuite)
{
    // Layout Right (default)
    Kratos::Future::mdarray_2d<int, Kratos::Future::layout_right> arr_lr(2, 3);
    int count_lr = 0;
    for(std::size_t i=0; i<2; ++i)
        for(std::size_t j=0; j<3; ++j)
            arr_lr(i,j) = count_lr++;
    // Expected data order: 0 1 2 / 3 4 5
    KRATOS_EXPECT_EQ(arr_lr.data()[0], 0); KRATOS_EXPECT_EQ(arr_lr.data()[1], 1); KRATOS_EXPECT_EQ(arr_lr.data()[2], 2);
    KRATOS_EXPECT_EQ(arr_lr.data()[3], 3); KRATOS_EXPECT_EQ(arr_lr.data()[4], 4); KRATOS_EXPECT_EQ(arr_lr.data()[5], 5);

    KRATOS_EXPECT_EQ(arr_lr(0,0), 0); KRATOS_EXPECT_EQ(arr_lr(0,1), 1); KRATOS_EXPECT_EQ(arr_lr(0,2), 2);
    KRATOS_EXPECT_EQ(arr_lr(1,0), 3); KRATOS_EXPECT_EQ(arr_lr(1,1), 4); KRATOS_EXPECT_EQ(arr_lr(1,2), 5);

    // Layout Left
    Kratos::Future::mdarray_2d<int, Kratos::Future::layout_left> arr_ll(2, 3);
    int count_ll = 0;
    for(std::size_t j=0; j<3; ++j) // Column major filling for clarity
        for(std::size_t i=0; i<2; ++i)
            arr_ll(i,j) = count_ll++;
    // Expected data order: 0 1 / 2 3 / 4 5
    KRATOS_EXPECT_EQ(arr_ll.data()[0], 0); KRATOS_EXPECT_EQ(arr_ll.data()[1], 1);
    KRATOS_EXPECT_EQ(arr_ll.data()[2], 2); KRATOS_EXPECT_EQ(arr_ll.data()[3], 3);
    KRATOS_EXPECT_EQ(arr_ll.data()[4], 4); KRATOS_EXPECT_EQ(arr_ll.data()[5], 5);

    KRATOS_EXPECT_EQ(arr_ll(0,0), 0); KRATOS_EXPECT_EQ(arr_ll(1,0), 1);
    KRATOS_EXPECT_EQ(arr_ll(0,1), 2); KRATOS_EXPECT_EQ(arr_ll(1,1), 3);
    KRATOS_EXPECT_EQ(arr_ll(0,2), 4); KRATOS_EXPECT_EQ(arr_ll(1,2), 5);
}

KRATOS_TEST_CASE_IN_SUITE(MdarrayContainerInteraction, KratosCoreFastSuite)
{
    Kratos::Future::mdarray_2d<double> arr(2,3);
    KRATOS_EXPECT_EQ(arr.container().size(), 2*3);
    KRATOS_EXPECT_EQ(arr.size(), 2*3); // mdarray::size() should be total elements

    arr(0,0) = 1.1;
    arr(1,2) = 6.6;

    KRATOS_EXPECT_DOUBLE_EQ(arr.data()[0], 1.1); // Assuming layout_right
    KRATOS_EXPECT_DOUBLE_EQ(arr.data()[5], 6.6); // Assuming layout_right

    // Check data() pointer consistency
    double* ptr = arr.data();
    ptr[1] = 2.2;
    KRATOS_EXPECT_DOUBLE_EQ(arr(0,1), 2.2);

    // Check container() reference consistency
    std::vector<double>& vec_ref = arr.container();
    vec_ref[2] = 3.3;
    KRATOS_EXPECT_DOUBLE_EQ(arr(0,2), 3.3);
}

KRATOS_TEST_CASE_IN_SUITE(MdarrayToMdspanConversion, KratosCoreFastSuite)
{
    Kratos::Future::mdarray_2d<int> arr(2,3);
    arr(0,0) = 1; arr(0,1) = 2; arr(0,2) = 3;
    arr(1,0) = 4; arr(1,1) = 5; arr(1,2) = 6;

    // Implicit conversion
    Kratos::Future::mdspan_2d_view<int> mds_view = arr;
    KRATOS_EXPECT_EQ(mds_view.rank(), 2);
    KRATOS_EXPECT_EQ(mds_view.extent(0), 2);
    KRATOS_EXPECT_EQ(mds_view.extent(1), 3);
    KRATOS_EXPECT_EQ(mds_view(0,0), 1);
    KRATOS_EXPECT_EQ(mds_view(1,2), 6);
    KRATOS_EXPECT_EQ(mds_view.data_handle(), arr.data());

    // Modify through mdspan view
    mds_view(0,1) = 22;
    KRATOS_EXPECT_EQ(arr(0,1), 22);

    // to_mdspan()
    auto mds_view_explicit = arr.to_mdspan();
    KRATOS_EXPECT_EQ(mds_view_explicit(0,0), 1);
    mds_view_explicit(1,1) = 55;
    KRATOS_EXPECT_EQ(arr(1,1), 55);

    // Const conversion
    const Kratos::Future::mdarray_2d<int>& const_arr = arr;
    Kratos::Future::mdspan_2d_view<const int> const_mds_view = const_arr;
    KRATOS_EXPECT_EQ(const_mds_view(1,1), 55);
    // const_mds_view(0,0) = 11; // This should not compile if uncommented

    auto const_mds_view_explicit = const_arr.to_mdspan();
    KRATOS_EXPECT_EQ(const_mds_view_explicit(1,1), 55);
}


KRATOS_TEST_CASE_IN_SUITE(MdarrayStaticExtentsWithStdArray, KratosCoreFastSuite)
{
    // mdarray with std::array container
    using StaticExtents = Kratos::Future::extents<std::size_t, 2, 2>;
    using ArrayContainer = std::array<int, 2*2>;
    Kratos::Future::mdarray<int, StaticExtents, Kratos::Future::layout_right, ArrayContainer> arr_std_array;

    arr_std_array(0,0) = 10;
    arr_std_array(0,1) = 20;
    arr_std_array(1,0) = 30;
    arr_std_array(1,1) = 40;

    KRATOS_EXPECT_EQ(arr_std_array(0,0), 10);
    KRATOS_EXPECT_EQ(arr_std_array(1,1), 40);
    KRATOS_EXPECT_EQ(arr_std_array.container()[0], 10); // Assuming layout_right

    // Test conversion to mdspan
    Kratos::Future::mdspan<int, StaticExtents> mds_view = arr_std_array;
    KRATOS_EXPECT_EQ(mds_view(1,0), 30);
}

} // namespace Kratos::Testing

#endif // KRATOS_FUTURE_MDSPAN
