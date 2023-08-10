//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/data_type_traits.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsInt, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<int>;

    static_assert(std::is_same_v<type_trait::PrimitiveDataType, int>);
    static_assert(type_trait::HasContiguousPrimitiveData);
    static_assert(!type_trait::HasDynamicMemoryAllocation);

    int test = 1;

    KRATOS_CHECK_EQUAL(type_trait::Size(test), 1U);
    KRATOS_CHECK_EQUAL(type_trait::Shape(test), std::vector<unsigned int>{});
    KRATOS_CHECK_IS_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{}));

    #ifdef KRATOS_DEBUG
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            type_trait::Reshape(test, std::vector<unsigned int>{1}),
            "Invalid shape given for primitive data type [ Expected shape = [], provided shape = [1] ]."
        );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsArray1dDouble, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<array_1d<double, 5>>;

    static_assert(std::is_same_v<type_trait::ContainerType, array_1d<double, 5>>);
    static_assert(std::is_same_v<type_trait::ValueType, double>);
    static_assert(std::is_same_v<type_trait::PrimitiveDataType, double>);
    static_assert(type_trait::HasContiguousPrimitiveData);
    static_assert(!type_trait::HasDynamicMemoryAllocation);

    array_1d<double, 5> test{};

    KRATOS_CHECK_EQUAL(type_trait::Size(test), 5);
    KRATOS_CHECK_EQUAL(type_trait::Shape(test), std::vector<unsigned int>{5});
    KRATOS_CHECK_IS_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{5}));

    #ifdef KRATOS_DEBUG
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            type_trait::Reshape(test, std::vector<unsigned int>{}),
            "Invalid shape given for array_1d data type [ Expected shape = [5], provided shape = [] ]."
        );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsVectorDouble, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<Vector>;

    static_assert(std::is_same_v<type_trait::ContainerType, Vector>);
    static_assert(std::is_same_v<type_trait::ValueType, double>);
    static_assert(std::is_same_v<type_trait::PrimitiveDataType, double>);
    static_assert(type_trait::HasContiguousPrimitiveData);
    static_assert(type_trait::HasDynamicMemoryAllocation);

    Vector test(7, 1);

    KRATOS_CHECK_EQUAL(type_trait::Size(test), 7);
    KRATOS_CHECK_EQUAL(type_trait::Shape(test), std::vector<unsigned int>{7});
    KRATOS_CHECK_IS_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{7}));
    KRATOS_CHECK(type_trait::Reshape(test, std::vector<unsigned int>{9}));
    KRATOS_CHECK_EQUAL(test.size(), 9);

    #ifdef KRATOS_DEBUG
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            type_trait::Reshape(test, std::vector<unsigned int>{}),
            "Invalid shape given for DenseVector data type."
        );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsMatrixDouble, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<Matrix>;

    static_assert(std::is_same_v<type_trait::ContainerType, Matrix>);
    static_assert(std::is_same_v<type_trait::ValueType, double>);
    static_assert(std::is_same_v<type_trait::PrimitiveDataType, double>);
    static_assert(type_trait::HasContiguousPrimitiveData);
    static_assert(type_trait::HasDynamicMemoryAllocation);

    Matrix test(3, 4, 1);

    KRATOS_CHECK_EQUAL(type_trait::Size(test), 12);

    std::vector<unsigned int> shape{3, 4};
    KRATOS_CHECK_EQUAL(type_trait::Shape(test), shape);
    KRATOS_CHECK_IS_FALSE(type_trait::Reshape(test, shape));

    shape[0] = 4;
    shape[1] = 5;
    KRATOS_CHECK(type_trait::Reshape(test, shape));
    KRATOS_CHECK_EQUAL(test.size1(), 4);
    KRATOS_CHECK_EQUAL(test.size2(), 5);

    #ifdef KRATOS_DEBUG
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            type_trait::Reshape(test, std::vector<unsigned int>{}),
            "Invalid shape given for DenseMatrix data type."
        );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsString, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<std::string>;

    static_assert(std::is_same_v<type_trait::ContainerType, std::string>);
    static_assert(std::is_same_v<type_trait::ValueType, char>);
    static_assert(std::is_same_v<type_trait::PrimitiveDataType, char>);
    static_assert(type_trait::HasContiguousPrimitiveData);
    static_assert(type_trait::HasDynamicMemoryAllocation);

    std::string test = "test";

    KRATOS_CHECK_EQUAL(type_trait::Size(test), 4);
    KRATOS_CHECK_EQUAL(type_trait::Shape(test), std::vector<unsigned int>{4});
    KRATOS_CHECK_IS_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{4}));
    KRATOS_CHECK(type_trait::Reshape(test, std::vector<unsigned int>{6}));
    KRATOS_CHECK_EQUAL(test.size(), 6);

    #ifdef KRATOS_DEBUG
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            type_trait::Reshape(test, std::vector<unsigned int>{}),
            "Invalid shape given for std::string data type."
        );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsStdVectorInt, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<std::vector<int>>;

    static_assert(std::is_same_v<type_trait::ContainerType, std::vector<int>>);
    static_assert(std::is_same_v<type_trait::ValueType, int>);
    static_assert(std::is_same_v<type_trait::PrimitiveDataType, int>);
    static_assert(type_trait::HasContiguousPrimitiveData);
    static_assert(type_trait::HasDynamicMemoryAllocation);

    std::vector<int> test(7, 1);

    KRATOS_CHECK_EQUAL(type_trait::Size(test), 7);
    KRATOS_CHECK_EQUAL(type_trait::Shape(test), std::vector<unsigned int>{7});
    KRATOS_CHECK_IS_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{7}));
    KRATOS_CHECK(type_trait::Reshape(test, std::vector<unsigned int>{9}));
    KRATOS_CHECK_EQUAL(test.size(), 9);

    #ifdef KRATOS_DEBUG
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            type_trait::Reshape(test, std::vector<unsigned int>{}),
            "Invalid shape given for std::vector data type."
        );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsArray1dStaticNested, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<array_1d<array_1d<array_1d<int, 10>, 3>, 5>>;

    static_assert(std::is_same_v<type_trait::ContainerType, array_1d<array_1d<array_1d<int, 10>, 3>, 5>>);
    static_assert(std::is_same_v<type_trait::ValueType, array_1d<array_1d<int, 10>, 3>>);
    static_assert(std::is_same_v<type_trait::PrimitiveDataType, int>);
    static_assert(!type_trait::HasContiguousPrimitiveData);
    static_assert(!type_trait::HasDynamicMemoryAllocation);

    array_1d<array_1d<array_1d<int, 10>, 3>, 5> test{};

    KRATOS_CHECK_EQUAL(type_trait::Size(test), 150);

    std::vector<unsigned int> shape{5, 3, 10};
    KRATOS_CHECK_EQUAL(type_trait::Shape(test), shape);
    KRATOS_CHECK_IS_FALSE(type_trait::Reshape(test, shape));

    #ifdef KRATOS_DEBUG
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            type_trait::Reshape(test, std::vector<unsigned int>{}),
            "Invalid shape given for array_1d data type [ Expected shape = [5], provided shape = [] ]."
        );
    #endif
}

} // namespace Kratos::Testing