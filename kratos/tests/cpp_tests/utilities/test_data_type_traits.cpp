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

    static_assert(std::is_same_v<type_trait::PrimitiveType, int>);
    static_assert(type_trait::Dimension == 0);
    static_assert(type_trait::IsContiguous);
    static_assert(!type_trait::IsDynamic);
    static_assert(!type_trait::IsDimensionDynamic<0>());

    int test = 1;

    KRATOS_EXPECT_EQ(type_trait::Size(test), 1U);
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), std::vector<unsigned int>{});
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{}));
    KRATOS_EXPECT_EQ(type_trait::GetContiguousData(test), &test);

    int dummy = -1;
    type_trait::CopyToContiguousData(&dummy, test);
    KRATOS_EXPECT_EQ(dummy, test);

    dummy = -2;
    type_trait::CopyFromContiguousData(dummy, &test);
    KRATOS_EXPECT_EQ(dummy, test);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        type_trait::Reshape(test, std::vector<unsigned int>{1}),
        "Invalid shape/dimension given for primitive data type [ Expected shape = [], provided shape = [1] ]."
    );
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsArray1dDouble, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<array_1d<double, 5>>;

    static_assert(std::is_same_v<type_trait::ContainerType, array_1d<double, 5>>);
    static_assert(std::is_same_v<type_trait::ValueType, double>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, double>);
    static_assert(type_trait::Dimension == 1);
    static_assert(type_trait::IsContiguous);
    static_assert(!type_trait::IsDynamic);
    static_assert(!type_trait::IsDimensionDynamic<0>());

    array_1d<double, 5> test{1, 2, 3, 4, 5};

    KRATOS_EXPECT_EQ(type_trait::Size(test), 5);
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), std::vector<unsigned int>{5});
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{5}));
    KRATOS_EXPECT_EQ(type_trait::GetContiguousData(test), &test[0]);

    array_1d<double, 5> dummy;
    type_trait::CopyToContiguousData(dummy.data().data(), test);
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    dummy = array_1d<double, 5>(5, -2);
    type_trait::CopyFromContiguousData(dummy, test.data().data());
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        type_trait::Reshape(test, std::vector<unsigned int>{}),
        "Invalid shape/dimension given for array_1d data type [ Expected shape = [5], provided shape = [] ]."
    );
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsVectorDouble, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<Vector>;

    static_assert(std::is_same_v<type_trait::ContainerType, Vector>);
    static_assert(std::is_same_v<type_trait::ValueType, double>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, double>);
    static_assert(type_trait::Dimension == 1);
    static_assert(type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::IsDimensionDynamic<0>());

    Vector test(7, 1);

    KRATOS_EXPECT_EQ(type_trait::Size(test), 7);
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), std::vector<unsigned int>{7});
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{7}));
    KRATOS_EXPECT_TRUE(type_trait::Reshape(test, std::vector<unsigned int>{9}));
    KRATOS_EXPECT_EQ(test.size(), 9);
    KRATOS_EXPECT_EQ(type_trait::GetContiguousData(test), &test[0]);

    test = Vector(9, 2);

    Vector dummy(9, -1);
    type_trait::CopyToContiguousData(dummy.data().begin(), test);
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    dummy = Vector(9, -2);
    type_trait::CopyFromContiguousData(dummy, test.data().begin());
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        type_trait::Reshape(test, std::vector<unsigned int>{}),
        "Invalid shape/dimension given for DenseVector data type [ Expected = [9], provided = [] ]"
    );
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsMatrixDouble, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<Matrix>;

    static_assert(std::is_same_v<type_trait::ContainerType, Matrix>);
    static_assert(std::is_same_v<type_trait::ValueType, double>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, double>);
    static_assert(type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 2);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(type_trait::IsDimensionDynamic<1>());

    Matrix test(3, 4, 1);

    KRATOS_EXPECT_EQ(type_trait::Size(test), 12);

    std::vector<unsigned int> shape{3, 4};
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), shape);
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, shape));

    shape[0] = 4;
    shape[1] = 5;
    KRATOS_EXPECT_TRUE(type_trait::Reshape(test, shape));
    KRATOS_EXPECT_EQ(test.size1(), 4);
    KRATOS_EXPECT_EQ(test.size2(), 5);
    KRATOS_EXPECT_EQ(type_trait::GetContiguousData(test), &test(0, 0));

    test = Matrix(4, 5, 2);

    Matrix dummy(4, 5);
    type_trait::CopyToContiguousData(dummy.data().begin(), test);
    KRATOS_EXPECT_MATRIX_EQ(dummy, test);

    dummy = Matrix(4, 5, -2);
    type_trait::CopyFromContiguousData(dummy, test.data().begin());
    KRATOS_EXPECT_MATRIX_EQ(dummy, test);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        type_trait::Reshape(test, std::vector<unsigned int>{}),
        "Invalid shape/dimension given for DenseMatrix data type [ Expected = [4, 5], provided = [] ]."
    );
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsString, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<std::string>;

    static_assert(std::is_same_v<type_trait::ContainerType, std::string>);
    static_assert(std::is_same_v<type_trait::ValueType, char>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, char>);
    static_assert(type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 1);
    static_assert(type_trait::IsDimensionDynamic<0>());

    std::string test = "test";

    KRATOS_EXPECT_EQ(type_trait::Size(test), 4);
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), std::vector<unsigned int>{4});
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{4}));
    KRATOS_EXPECT_TRUE(type_trait::Reshape(test, std::vector<unsigned int>{6}));
    KRATOS_EXPECT_EQ(test.size(), 6);
    KRATOS_EXPECT_EQ(type_trait::GetContiguousData(test), &test[0]);

    test = "test01";

    std::string dummy = "000000";
    type_trait::CopyToContiguousData(dummy.data(), test);
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    dummy = "000000";
    type_trait::CopyFromContiguousData(dummy, test.data());
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        type_trait::Reshape(test, std::vector<unsigned int>{}),
        "Invalid shape/dimension given for std::string data type [ Expected = [6], provided = [] ]."
    );
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsStdVectorInt, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<std::vector<int>>;

    static_assert(std::is_same_v<type_trait::ContainerType, std::vector<int>>);
    static_assert(std::is_same_v<type_trait::ValueType, int>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, int>);
    static_assert(type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 1);
    static_assert(type_trait::IsDimensionDynamic<0>());

    std::vector<int> test(7, 1);

    KRATOS_EXPECT_EQ(type_trait::Size(test), 7);
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), std::vector<unsigned int>{7});
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, std::vector<unsigned int>{7}));
    KRATOS_EXPECT_TRUE(type_trait::Reshape(test, std::vector<unsigned int>{9}));
    KRATOS_EXPECT_EQ(test.size(), 9);
    KRATOS_EXPECT_EQ(type_trait::GetContiguousData(test), &test[0]);

    test = std::vector<int>(9, 2);

    std::vector<int> dummy(9);
    type_trait::CopyToContiguousData(dummy.data(), test);
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    dummy = std::vector<int>(9, -2);
    type_trait::CopyFromContiguousData(dummy, test.data());
    KRATOS_EXPECT_VECTOR_EQ(dummy, test);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        type_trait::Reshape(test, std::vector<unsigned int>{}),
        "Invalid shape/dimension given for std::vector data type [ Expected = [9], provided = [] ]."
    );
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsArray1dNested, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<array_1d<array_1d<array_1d<int, 10>, 3>, 5>>;

    static_assert(std::is_same_v<type_trait::ContainerType, array_1d<array_1d<array_1d<int, 10>, 3>, 5>>);
    static_assert(std::is_same_v<type_trait::ValueType, array_1d<array_1d<int, 10>, 3>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, int>);
    static_assert(type_trait::IsContiguous);
    static_assert(!type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 3);
    static_assert(!type_trait::IsDimensionDynamic<0>());
    static_assert(!type_trait::IsDimensionDynamic<1>());
    static_assert(!type_trait::IsDimensionDynamic<2>());

    array_1d<array_1d<array_1d<int, 10>, 3>, 5> test{}, result{};
    std::vector<int> ref_values(150);
    for (unsigned int i = 0; i < 5; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 10; ++k) {
                test[i][j][k] = (i+1)*(j+1)*(k+1);
                ref_values[i * 30 + j * 10 + k] = test[i][j][k];
            }
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(test), 150);

    std::vector<unsigned int> shape{5, 3, 10};
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), shape);
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, shape));

    std::vector<int> values(150, -1);
    type_trait::CopyToContiguousData(values.data(), test);
    KRATOS_EXPECT_EQ(values, ref_values);

    type_trait::CopyFromContiguousData(result, values.data());
    KRATOS_EXPECT_EQ(test, result);
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsVectorNested, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<DenseVector<DenseVector<Vector>>>;

    static_assert(std::is_same_v<type_trait::ContainerType, DenseVector<DenseVector<Vector>>>);
    static_assert(std::is_same_v<type_trait::ValueType, DenseVector<Vector>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, double>);
    static_assert(!type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 3);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(type_trait::IsDimensionDynamic<1>());
    static_assert(type_trait::IsDimensionDynamic<2>());

    DenseVector<DenseVector<Vector>> test(2), result(2);
    std::vector<double> ref_values(24);
    for (unsigned int i = 0; i < 2; ++i) {
        test[i] = DenseVector<Vector>(3);
        result[i] = DenseVector<Vector>(3);
        for (unsigned int j = 0; j < 3; ++j) {
            test[i][j] = Vector(4);
            result[i][j] = Vector(4, -1);
            for (unsigned int k = 0; k < 4; ++k) {
                test[i][j][k] = (i+1)*(j+1)*(k+1);
                ref_values[i * 12 + j * 4 + k] = test[i][j][k];
            }
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(test), 24);

    std::vector<unsigned int> shape{2, 3, 4};
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), shape);
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, shape));

    auto temp = test;
    shape[0] = 5;
    shape[1] = 6;
    shape[2] = 7;
    KRATOS_EXPECT_TRUE(type_trait::Reshape(temp, shape));
    KRATOS_EXPECT_EQ(temp.size(), 5);
    for (unsigned int i = 0; i < 5; ++i) {
        KRATOS_EXPECT_EQ(temp[i].size(), 6);
        for (unsigned int j = 0; j < 6; ++j) {
            KRATOS_EXPECT_EQ(temp[i][j].size(), 7);
        }
    }

    std::vector<double> values(24, -1);
    type_trait::CopyToContiguousData(values.data(), test);
    KRATOS_EXPECT_EQ(values, ref_values);

    type_trait::CopyFromContiguousData(result, values.data());
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            KRATOS_EXPECT_VECTOR_EQ(result[i][j], test[i][j]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsMatrixNested, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<DenseMatrix<DenseMatrix<Matrix>>>;

    static_assert(std::is_same_v<type_trait::ContainerType, DenseMatrix<DenseMatrix<Matrix>>>);
    static_assert(std::is_same_v<type_trait::ValueType, DenseMatrix<Matrix>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, double>);
    static_assert(!type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 6);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(type_trait::IsDimensionDynamic<1>());
    static_assert(type_trait::IsDimensionDynamic<2>());
    static_assert(type_trait::IsDimensionDynamic<3>());
    static_assert(type_trait::IsDimensionDynamic<4>());
    static_assert(type_trait::IsDimensionDynamic<5>());

    DenseMatrix<DenseMatrix<Matrix>> test(2, 3), result(2, 3);
    std::vector<double> ref_values(5040);
    for (unsigned int i = 0; i < 6; ++i) {
        test.data()[i] = DenseMatrix<Matrix>(4, 5);
        result.data()[i] = DenseMatrix<Matrix>(4, 5);
        for (unsigned int j = 0; j < 20; ++j) {
            test.data()[i].data()[j] = Matrix(6, 7);
            result.data()[i].data()[j] = Matrix(6, 7, -1);
            for (unsigned int k = 0; k < 42; ++k) {
                test.data()[i].data()[j].data()[k] = (i+1)*(j+1)*(k+1);
                ref_values[i * 840 + j * 42 + k] = (i+1)*(j+1)*(k+1);
            }
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(test), 2*3*4*5*6*7);

    std::vector<unsigned int> shape{2, 3, 4, 5, 6, 7};
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), shape);
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, shape));

    auto temp = test;
    shape[0] = 3;
    shape[1] = 4;
    shape[2] = 5;
    shape[3] = 6;
    shape[4] = 7;
    shape[5] = 8;
    KRATOS_EXPECT_TRUE(type_trait::Reshape(temp, shape));
    KRATOS_EXPECT_EQ(temp.size1(), 3);
    KRATOS_EXPECT_EQ(temp.size2(), 4);
    for (unsigned int i = 0; i < 12; ++i) {
        KRATOS_EXPECT_EQ(temp.data()[i].size1(), 5);
        KRATOS_EXPECT_EQ(temp.data()[i].size2(), 6);
        for (unsigned int j = 0; j < 30; ++j) {
            KRATOS_EXPECT_EQ(temp.data()[i].data()[j].size1(), 7);
            KRATOS_EXPECT_EQ(temp.data()[i].data()[j].size2(), 8);
        }
    }

    std::vector<double> values(5040, -1);
    type_trait::CopyToContiguousData(values.data(), test);
    KRATOS_EXPECT_EQ(values, ref_values);

    type_trait::CopyFromContiguousData(result, values.data());
    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 20; ++j) {
            KRATOS_EXPECT_MATRIX_EQ(result.data()[i].data()[j], test.data()[i].data()[j]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsStdVectorNested, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<std::vector<std::vector<int>>>;

    static_assert(std::is_same_v<type_trait::ContainerType, std::vector<std::vector<int>>>);
    static_assert(std::is_same_v<type_trait::ValueType, std::vector<int>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, int>);
    static_assert(!type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 2);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(type_trait::IsDimensionDynamic<1>());

    std::vector<std::vector<int>> test(2), result(2);
    std::vector<int> ref_values(6);
    for (unsigned int i = 0; i < 2; ++i) {
        test[i] = std::vector<int>(3);
        result[i] = std::vector<int>(3, -1);
        for (unsigned int j = 0; j < 3; ++j) {
            test[i][j] = (i+1)*(j+1);
            ref_values[i * 3 + j] = test[i][j];
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(test), 6);

    std::vector<unsigned int> shape{2, 3};
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), shape);
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, shape));

    auto temp = test;
    shape[0] = 3;
    shape[1] = 4;
    KRATOS_EXPECT_TRUE(type_trait::Reshape(temp, shape));
    KRATOS_EXPECT_EQ(temp.size(), 3);
    for (unsigned int i = 0; i < 3; ++i) {
        KRATOS_EXPECT_EQ(temp[i].size(), 4);
    }

    std::vector<int> values(6, -1);
    type_trait::CopyToContiguousData(values.data(), test);
    KRATOS_EXPECT_VECTOR_EQ(values, ref_values);

    type_trait::CopyFromContiguousData(result, values.data());
    for (unsigned int i = 0; i < 2; ++i) {
        KRATOS_EXPECT_VECTOR_EQ(result[i], test[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsMixedNested1, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<std::vector<DenseVector<DenseMatrix<array_1d<std::string, 6>>>>>;

    static_assert(std::is_same_v<type_trait::ContainerType, std::vector<DenseVector<DenseMatrix<array_1d<std::string, 6>>>>>);
    static_assert(std::is_same_v<type_trait::ValueType, DenseVector<DenseMatrix<array_1d<std::string, 6>>>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, char>);
    static_assert(!type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 6);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(type_trait::IsDimensionDynamic<1>());
    static_assert(type_trait::IsDimensionDynamic<2>());
    static_assert(type_trait::IsDimensionDynamic<3>());
    static_assert(!type_trait::IsDimensionDynamic<4>());
    static_assert(type_trait::IsDimensionDynamic<5>());

    type_trait::ContainerType test(3), result(3);
    std::vector<char> ref_values(2160);
    for (unsigned int i = 0; i < 3; ++i) {
        // loop of std::vector
        test[i] = DenseVector<DenseMatrix<array_1d<std::string, 6>>>(4);
        result[i] = DenseVector<DenseMatrix<array_1d<std::string, 6>>>(4);
        for (unsigned int j = 0; j < 4; ++j) {
            // loop of DenseVector
            test[i][j] = DenseMatrix<array_1d<std::string, 6>>(5, 6);
            result[i][j] = DenseMatrix<array_1d<std::string, 6>>(5, 6);
            for (unsigned int k = 0; k < 30; ++k) {
                test[i][j].data()[k] = array_1d<std::string, 6>(6, std::to_string(k % 10));
                result[i][j].data()[k] = array_1d<std::string ,6>(6, "a");
                std::fill(ref_values.begin() + i * 720 + j * 180 + k * 6, ref_values.begin() + i * 720 + j * 180 + k * 6 + 6, test[i][j].data()[k][0][0]);
            }
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(test), 2160);

    std::vector<unsigned int> shape{3, 4, 5, 6, 6, 1};
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), shape);
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, shape));

    auto temp = test;
    shape[0] = 4;
    shape[1] = 5;
    shape[2] = 6;
    shape[3] = 7;
    shape[4] = 6;
    shape[5] = 9;
    KRATOS_EXPECT_TRUE(type_trait::Reshape(temp, shape));
    KRATOS_EXPECT_EQ(temp.size(), 4);
    for (unsigned int i = 0; i < 4; ++i) {
        KRATOS_EXPECT_EQ(temp[i].size(), 5);
        for (unsigned int j = 0; j < 5; ++j) {
            KRATOS_EXPECT_EQ(temp[i][j].size1(), 6);
            KRATOS_EXPECT_EQ(temp[i][j].size2(), 7);
            for (unsigned int k = 0; k < 42; ++k) {
                KRATOS_EXPECT_EQ(temp[i][j].data()[k].size(), 6);
                for (unsigned int l = 0; l < 6; ++l) {
                    KRATOS_EXPECT_EQ(temp[i][j].data()[k][l].size(), 9);
                }
            }
        }
    }

    std::vector<char> values(2160, 0);
    type_trait::CopyToContiguousData(values.data(), test);
    KRATOS_EXPECT_EQ(values, ref_values);

    type_trait::CopyFromContiguousData(result, values.data());
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            for (unsigned int k = 0; k < 30; ++k) {
                for (unsigned int l = 0; l < 6; ++l) {
                    KRATOS_EXPECT_EQ(result[i][j].data()[k][l], test[i][j].data()[k][l]);
                }
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsMixedNested2, KratosCoreFastSuite)
{
    using type_trait = DataTypeTraits<std::vector<DenseVector<DenseMatrix<array_1d<double, 6>>>>>;

    static_assert(std::is_same_v<type_trait::ContainerType, std::vector<DenseVector<DenseMatrix<array_1d<double, 6>>>>>);
    static_assert(std::is_same_v<type_trait::ValueType, DenseVector<DenseMatrix<array_1d<double, 6>>>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, double>);
    static_assert(!type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 5);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(type_trait::IsDimensionDynamic<1>());
    static_assert(type_trait::IsDimensionDynamic<2>());
    static_assert(type_trait::IsDimensionDynamic<3>());
    static_assert(!type_trait::IsDimensionDynamic<4>());

    type_trait::ContainerType test(3), result(3);
    std::vector<double> ref_values(2160);
    for (unsigned int i = 0; i < 3; ++i) {
        // loop of std::vector
        test[i] = DenseVector<DenseMatrix<array_1d<double, 6>>>(4);
        result[i] = DenseVector<DenseMatrix<array_1d<double, 6>>>(4);
        for (unsigned int j = 0; j < 4; ++j) {
            // loop of DenseVector
            test[i][j] = DenseMatrix<array_1d<double, 6>>(5, 6);
            result[i][j] = DenseMatrix<array_1d<double, 6>>(5, 6);
            for (unsigned int k = 0; k < 30; ++k) {
                test[i][j].data()[k] = array_1d<double, 6>(6);
                result[i][j].data()[k] = array_1d<double ,6>(6, -1);
                for (unsigned int l = 0; l < 6; ++l) {
                    test[i][j].data()[k][l] = (i+1) * (j+1) * (k+1) * (l+1);
                    ref_values[i * 720 + j * 180 + k * 6 + l] = test[i][j].data()[k][l];
                }
            }
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(test), 2160);

    std::vector<unsigned int> shape{3, 4, 5, 6, 6};
    KRATOS_EXPECT_EQ(type_trait::Shape(test).size(), type_trait::Dimension);
    KRATOS_EXPECT_EQ(type_trait::Shape(test), shape);
    KRATOS_EXPECT_FALSE(type_trait::Reshape(test, shape));

    auto temp = test;
    shape[0] = 4;
    shape[1] = 5;
    shape[2] = 6;
    shape[3] = 7;
    shape[4] = 6;
    KRATOS_EXPECT_TRUE(type_trait::Reshape(temp, shape));
    KRATOS_EXPECT_EQ(temp.size(), 4);
    for (unsigned int i = 0; i < 4; ++i) {
        KRATOS_EXPECT_EQ(temp[i].size(), 5);
        for (unsigned int j = 0; j < 5; ++j) {
            KRATOS_EXPECT_EQ(temp[i][j].size1(), 6);
            KRATOS_EXPECT_EQ(temp[i][j].size2(), 7);
            for (unsigned int k = 0; k < 42; ++k) {
                KRATOS_EXPECT_EQ(temp[i][j].data()[k].size(), 6);
            }
        }
    }

    std::vector<double> values(2160, 0);
    type_trait::CopyToContiguousData(values.data(), test);
    KRATOS_EXPECT_EQ(values, ref_values);

    type_trait::CopyFromContiguousData(result, values.data());
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            for (unsigned int k = 0; k < 30; ++k) {
                for (unsigned int l = 0; l < 6; ++l) {
                    KRATOS_EXPECT_EQ(result[i][j].data()[k][l], test[i][j].data()[k][l]);
                }
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsNestedArray1d, KratosCoreFastSuite)
{
    using data_static_type = array_1d<array_1d<array_1d<int, 4>, 5>, 6>;

    using static_type_trait = DataTypeTraits<data_static_type>;

    static_assert(std::is_same_v<static_type_trait::ContainerType, data_static_type>);
    static_assert(std::is_same_v<static_type_trait::ValueType, array_1d<array_1d<int, 4>, 5>>);
    static_assert(std::is_same_v<static_type_trait::PrimitiveType, int>);
    static_assert(static_type_trait::IsContiguous);
    static_assert(!static_type_trait::IsDynamic);
    static_assert(static_type_trait::Dimension == 3);
    static_assert(!static_type_trait::IsDimensionDynamic<0>());
    static_assert(!static_type_trait::IsDimensionDynamic<1>());
    static_assert(!static_type_trait::IsDimensionDynamic<2>());

    data_static_type static_test;
    std::vector<int> ref_values(120, -1);
    unsigned int local_index = 0;
    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 5; ++j) {
            for (unsigned int k = 0; k < 4; ++k) {
                static_test[i][j][k] = (i+1)*(j+1)*(k+1);
                ref_values[local_index++] = static_test[i][j][k];
            }
        }
    }

    KRATOS_EXPECT_EQ(static_type_trait::Size(static_test), 120);
    const auto const_static_test = static_test;
    int* p_static_start = static_type_trait::GetContiguousData(static_test);
    int const* p_const_static_start =  static_type_trait::GetContiguousData(const_static_test);
    for (unsigned int i = 0; i < 120; ++i) {
        KRATOS_EXPECT_EQ(ref_values[i], *(p_static_start++));
        KRATOS_EXPECT_EQ(ref_values[i], *(p_const_static_start++));
    }

    // now check for non static data type combinations
    static_assert(!DataTypeTraits<array_1d<array_1d<Vector, 3>, 4>>::IsContiguous);
    static_assert(DataTypeTraits<array_1d<array_1d<Vector, 3>, 4>>::IsDynamic);
    static_assert(std::is_same_v<DataTypeTraits<array_1d<array_1d<Vector, 3>, 4>>::PrimitiveType, double>);
    static_assert(!DataTypeTraits<array_1d<array_1d<Vector, 3>, 4>>::IsDimensionDynamic<0>());
    static_assert(!DataTypeTraits<array_1d<array_1d<Vector, 3>, 4>>::IsDimensionDynamic<1>());
    static_assert(DataTypeTraits<array_1d<array_1d<Vector, 3>, 4>>::IsDimensionDynamic<2>());

    static_assert(!DataTypeTraits<array_1d<array_1d<DenseVector<array_1d<int, 2>>, 3>, 4>>::IsContiguous);
    static_assert(DataTypeTraits<array_1d<array_1d<DenseVector<array_1d<int, 2>>, 3>, 4>>::IsDynamic);
    static_assert(std::is_same_v<DataTypeTraits<array_1d<array_1d<DenseVector<array_1d<int, 2>>, 3>, 4>>::PrimitiveType, int>);
    static_assert(!DataTypeTraits<array_1d<array_1d<DenseVector<array_1d<int, 2>>, 3>, 4>>::IsDimensionDynamic<0>());
    static_assert(!DataTypeTraits<array_1d<array_1d<DenseVector<array_1d<int, 2>>, 3>, 4>>::IsDimensionDynamic<1>());
    static_assert(DataTypeTraits<array_1d<array_1d<DenseVector<array_1d<int, 2>>, 3>, 4>>::IsDimensionDynamic<2>());
    static_assert(!DataTypeTraits<array_1d<array_1d<DenseVector<array_1d<int, 2>>, 3>, 4>>::IsDimensionDynamic<3>());
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsNestedDenseVector, KratosCoreFastSuite)
{
    using data_type = DenseVector<array_1d<array_1d<int, 4>, 5>>;

    using type_trait = DataTypeTraits<data_type>;

    static_assert(std::is_same_v<type_trait::ContainerType, data_type>);
    static_assert(std::is_same_v<type_trait::ValueType, array_1d<array_1d<int, 4>, 5>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, int>);
    static_assert(type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 3);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(!type_trait::IsDimensionDynamic<1>());
    static_assert(!type_trait::IsDimensionDynamic<2>());

    data_type static_test(6);
    std::vector<int> ref_values(120, -1);
    unsigned int local_index = 0;
    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 5; ++j) {
            for (unsigned int k = 0; k < 4; ++k) {
                static_test[i][j][k] = (i+1)*(j+1)*(k+1);
                ref_values[local_index++] = static_test[i][j][k];
            }
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(static_test), 120);
    const auto const_static_test = static_test;
    int* p_static_start = type_trait::GetContiguousData(static_test);
    int const* p_const_static_start =  type_trait::GetContiguousData(const_static_test);
    for (unsigned int i = 0; i < 120; ++i) {
        KRATOS_EXPECT_EQ(ref_values[i], *(p_static_start++));
        KRATOS_EXPECT_EQ(ref_values[i], *(p_const_static_start++));
    }

    // now check for non static data type combinations
    static_assert(!DataTypeTraits<DenseVector<array_1d<Vector, 3>>>::IsContiguous);
    static_assert(DataTypeTraits<DenseVector<array_1d<Vector, 3>>>::IsDynamic);
    static_assert(std::is_same_v<DataTypeTraits<DenseVector<array_1d<Vector, 3>>>::PrimitiveType, double>);
    static_assert(DataTypeTraits<DenseVector<array_1d<Vector, 3>>>::IsDimensionDynamic<0>());
    static_assert(!DataTypeTraits<DenseVector<array_1d<Vector, 3>>>::IsDimensionDynamic<1>());
    static_assert(DataTypeTraits<DenseVector<array_1d<Vector, 3>>>::IsDimensionDynamic<2>());

    static_assert(!DataTypeTraits<DenseVector<array_1d<DenseVector<array_1d<int, 2>>, 3>>>::IsContiguous);
    static_assert(DataTypeTraits<DenseVector<array_1d<DenseVector<array_1d<int, 2>>, 3>>>::IsDynamic);
    static_assert(std::is_same_v<DataTypeTraits<DenseVector<array_1d<DenseVector<array_1d<int, 2>>, 3>>>::PrimitiveType, int>);
    static_assert(DataTypeTraits<DenseVector<array_1d<DenseVector<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<0>());
    static_assert(!DataTypeTraits<DenseVector<array_1d<DenseVector<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<1>());
    static_assert(DataTypeTraits<DenseVector<array_1d<DenseVector<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<2>());
    static_assert(!DataTypeTraits<DenseVector<array_1d<DenseVector<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<3>());
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsNestedDenseMatrix, KratosCoreFastSuite)
{
    using data_type = DenseMatrix<array_1d<array_1d<int, 4>, 5>>;

    using type_trait = DataTypeTraits<data_type>;

    static_assert(std::is_same_v<type_trait::ContainerType, data_type>);
    static_assert(std::is_same_v<type_trait::ValueType, array_1d<array_1d<int, 4>, 5>>);
    static_assert(std::is_same_v<type_trait::PrimitiveType, int>);
    static_assert(type_trait::IsContiguous);
    static_assert(type_trait::IsDynamic);
    static_assert(type_trait::Dimension == 4);
    static_assert(type_trait::IsDimensionDynamic<0>());
    static_assert(type_trait::IsDimensionDynamic<1>());
    static_assert(!type_trait::IsDimensionDynamic<2>());
    static_assert(!type_trait::IsDimensionDynamic<3>());

    data_type static_test(2, 3);
    std::vector<int> ref_values(120, -1);
    unsigned int local_index = 0;
    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 5; ++j) {
            for (unsigned int k = 0; k < 4; ++k) {
                static_test.data()[i][j][k] = (i+1)*(j+1)*(k+1);
                ref_values[local_index++] = static_test.data()[i][j][k];
            }
        }
    }

    KRATOS_EXPECT_EQ(type_trait::Size(static_test), 120);
    const auto const_static_test = static_test;
    int* p_static_start = type_trait::GetContiguousData(static_test);
    int const* p_const_static_start =  type_trait::GetContiguousData(const_static_test);
    for (unsigned int i = 0; i < 120; ++i) {
        KRATOS_EXPECT_EQ(ref_values[i], *(p_static_start++));
        KRATOS_EXPECT_EQ(ref_values[i], *(p_const_static_start++));
    }

    // now check for non static data type combinations
    static_assert(!DataTypeTraits<DenseMatrix<array_1d<Vector, 3>>>::IsContiguous);
    static_assert(DataTypeTraits<DenseMatrix<array_1d<Vector, 3>>>::IsDynamic);
    static_assert(std::is_same_v<DataTypeTraits<DenseMatrix<array_1d<Vector, 3>>>::PrimitiveType, double>);
    static_assert(DataTypeTraits<DenseMatrix<array_1d<Vector, 3>>>::IsDimensionDynamic<0>());
    static_assert(DataTypeTraits<DenseMatrix<array_1d<Vector, 3>>>::IsDimensionDynamic<1>());
    static_assert(!DataTypeTraits<DenseMatrix<array_1d<Vector, 3>>>::IsDimensionDynamic<2>());
    static_assert(DataTypeTraits<DenseMatrix<array_1d<Vector, 3>>>::IsDimensionDynamic<3>());

    static_assert(!DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsContiguous);
    static_assert(DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsDynamic);
    static_assert(std::is_same_v<DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::PrimitiveType, int>);
    static_assert(DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<0>());
    static_assert(DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<1>());
    static_assert(!DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<2>());
    static_assert(DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<3>());
    static_assert(DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<4>());
    static_assert(!DataTypeTraits<DenseMatrix<array_1d<DenseMatrix<array_1d<int, 2>>, 3>>>::IsDimensionDynamic<5>());
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsNestedZeroSize, KratosCoreFastSuite)
{
    using data_type_1 = std::vector<DenseVector<DenseMatrix<array_1d<std::string, 0>>>>;
    using data_type_1_trait = DataTypeTraits<data_type_1>;
    static_assert(!data_type_1_trait::IsContiguous);
    static_assert(data_type_1_trait::IsDynamic);
    static_assert(data_type_1_trait::Dimension == 6);
    static_assert(data_type_1_trait::IsDimensionDynamic<0>());
    static_assert(data_type_1_trait::IsDimensionDynamic<1>());
    static_assert(data_type_1_trait::IsDimensionDynamic<2>());
    static_assert(data_type_1_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_1_trait::IsDimensionDynamic<4>());
    static_assert(data_type_1_trait::IsDimensionDynamic<5>());
    data_type_1 test_data_type_1;
    KRATOS_EXPECT_EQ(data_type_1_trait::Size(test_data_type_1), 0);
    std::vector<unsigned int> data_type_1_shape{0, 0, 0, 0, 0, 0};
    KRATOS_EXPECT_EQ(data_type_1_trait::Shape(test_data_type_1), data_type_1_shape);

    using data_type_2 = std::vector<DenseVector<DenseMatrix<array_1d<std::string, 3>>>>;
    using data_type_2_trait = DataTypeTraits<data_type_2>;
    static_assert(!data_type_2_trait::IsContiguous);
    static_assert(data_type_2_trait::IsDynamic);
    static_assert(data_type_2_trait::Dimension == 6);
    static_assert(data_type_2_trait::IsDimensionDynamic<0>());
    static_assert(data_type_2_trait::IsDimensionDynamic<1>());
    static_assert(data_type_2_trait::IsDimensionDynamic<2>());
    static_assert(data_type_2_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_2_trait::IsDimensionDynamic<4>());
    static_assert(data_type_2_trait::IsDimensionDynamic<5>());
    data_type_2 test_data_type_2;
    KRATOS_EXPECT_EQ(data_type_2_trait::Size(test_data_type_2), 0);
    std::vector<unsigned int> data_type_2_shape{0, 0, 0, 0, 3, 0};
    KRATOS_EXPECT_EQ(data_type_2_trait::Shape(test_data_type_2), data_type_2_shape);

    using data_type_3 = std::vector<array_1d<array_1d<array_1d<array_1d<double, 4>, 0>, 4>, 0>>;
    using data_type_3_trait = DataTypeTraits<data_type_3>;
    static_assert(!data_type_3_trait::IsContiguous);
    static_assert(data_type_3_trait::IsDynamic);
    static_assert(data_type_3_trait::Dimension == 5);
    static_assert(data_type_3_trait::IsDimensionDynamic<0>());
    static_assert(!data_type_3_trait::IsDimensionDynamic<1>());
    static_assert(!data_type_3_trait::IsDimensionDynamic<2>());
    static_assert(!data_type_3_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_3_trait::IsDimensionDynamic<4>());
    data_type_3 test_data_type_3;
    KRATOS_EXPECT_EQ(data_type_3_trait::Size(test_data_type_3), 0);
    std::vector<unsigned int> data_type_3_shape{0, 0, 4, 0, 4};
    KRATOS_EXPECT_EQ(data_type_3_trait::Shape(test_data_type_3), data_type_3_shape);

    using data_type_4 = DenseVector<array_1d<array_1d<array_1d<array_1d<double, 4>, 0>, 4>, 0>>;
    using data_type_4_trait = DataTypeTraits<data_type_4>;
    static_assert(data_type_4_trait::IsContiguous);
    static_assert(data_type_4_trait::IsDynamic);
    static_assert(data_type_4_trait::Dimension == 5);
    static_assert(data_type_4_trait::IsDimensionDynamic<0>());
    static_assert(!data_type_4_trait::IsDimensionDynamic<1>());
    static_assert(!data_type_4_trait::IsDimensionDynamic<2>());
    static_assert(!data_type_4_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_4_trait::IsDimensionDynamic<4>());
    data_type_4 test_data_type_4;
    KRATOS_EXPECT_EQ(data_type_4_trait::Size(test_data_type_4), 0);
    std::vector<unsigned int> data_type_4_shape{0, 0, 4, 0, 4};
    KRATOS_EXPECT_EQ(data_type_4_trait::Shape(test_data_type_4), data_type_4_shape);
    KRATOS_EXPECT_TRUE(data_type_4_trait::GetContiguousData(test_data_type_4) == nullptr);

    using data_type_5 = DenseMatrix<array_1d<array_1d<array_1d<array_1d<double, 4>, 0>, 4>, 0>>;
    using data_type_5_trait = DataTypeTraits<data_type_5>;
    static_assert(data_type_5_trait::IsContiguous);
    static_assert(data_type_5_trait::IsDynamic);
    static_assert(data_type_5_trait::Dimension == 6);
    static_assert(data_type_5_trait::IsDimensionDynamic<0>());
    static_assert(data_type_5_trait::IsDimensionDynamic<1>());
    static_assert(!data_type_5_trait::IsDimensionDynamic<2>());
    static_assert(!data_type_5_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_5_trait::IsDimensionDynamic<4>());
    static_assert(!data_type_5_trait::IsDimensionDynamic<5>());
    data_type_5 test_data_type_5;
    KRATOS_EXPECT_EQ(data_type_5_trait::Size(test_data_type_5), 0);
    std::vector<unsigned int> data_type_5_shape{0, 0, 0, 4, 0, 4};
    KRATOS_EXPECT_EQ(data_type_5_trait::Shape(test_data_type_5), data_type_5_shape);
    KRATOS_EXPECT_TRUE(data_type_5_trait::GetContiguousData(test_data_type_5) == nullptr);

    using data_type_6 = array_1d<array_1d<array_1d<array_1d<array_1d<double, 4>, 0>, 4>, 0>, 2>;
    using data_type_6_trait = DataTypeTraits<data_type_6>;
    static_assert(data_type_6_trait::IsContiguous);
    static_assert(!data_type_6_trait::IsDynamic);
    static_assert(data_type_6_trait::Dimension == 5);
    static_assert(!data_type_6_trait::IsDimensionDynamic<0>());
    static_assert(!data_type_6_trait::IsDimensionDynamic<1>());
    static_assert(!data_type_6_trait::IsDimensionDynamic<2>());
    static_assert(!data_type_6_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_6_trait::IsDimensionDynamic<4>());
    data_type_6 test_data_type_6;
    KRATOS_EXPECT_EQ(data_type_6_trait::Size(test_data_type_6), 0);
    std::vector<unsigned int> data_type_6_shape{2, 0, 4, 0, 4};
    KRATOS_EXPECT_EQ(data_type_6_trait::Shape(test_data_type_6), data_type_6_shape);

    using data_type_7 = array_1d<array_1d<DenseVector<array_1d<array_1d<double, 4>, 0>>, 0>, 2>;
    using data_type_7_trait = DataTypeTraits<data_type_7>;
    static_assert(!data_type_7_trait::IsContiguous);
    static_assert(data_type_7_trait::IsDynamic);
    static_assert(data_type_7_trait::Dimension == 5);
    static_assert(!data_type_7_trait::IsDimensionDynamic<0>());
    static_assert(!data_type_7_trait::IsDimensionDynamic<1>());
    static_assert(data_type_7_trait::IsDimensionDynamic<2>());
    static_assert(!data_type_7_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_7_trait::IsDimensionDynamic<4>());
    data_type_7 test_data_type_7;
    KRATOS_EXPECT_EQ(data_type_7_trait::Size(test_data_type_7), 0);
    std::vector<unsigned int> data_type_7_shape{2, 0, 0, 0, 4};
    KRATOS_EXPECT_EQ(data_type_7_trait::Shape(test_data_type_7), data_type_7_shape);

    using data_type_8 = array_1d<array_1d<DenseMatrix<array_1d<array_1d<double, 4>, 0>>, 0>, 2>;
    using data_type_8_trait = DataTypeTraits<data_type_8>;
    static_assert(!data_type_8_trait::IsContiguous);
    static_assert(data_type_8_trait::IsDynamic);
    static_assert(data_type_8_trait::Dimension == 6);
    static_assert(!data_type_8_trait::IsDimensionDynamic<0>());
    static_assert(!data_type_8_trait::IsDimensionDynamic<1>());
    static_assert(data_type_8_trait::IsDimensionDynamic<2>());
    static_assert(data_type_8_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_8_trait::IsDimensionDynamic<4>());
    static_assert(!data_type_8_trait::IsDimensionDynamic<5>());
    data_type_8 test_data_type_8;
    KRATOS_EXPECT_EQ(data_type_8_trait::Size(test_data_type_8), 0);
    std::vector<unsigned int> data_type_8_shape{2, 0, 0, 0, 0, 4};
    KRATOS_EXPECT_EQ(data_type_8_trait::Shape(test_data_type_8), data_type_8_shape);

    using data_type_9 = array_1d<array_1d<std::vector<array_1d<array_1d<double, 4>, 0>>, 0>, 2>;
    using data_type_9_trait = DataTypeTraits<data_type_9>;
    static_assert(!data_type_9_trait::IsContiguous);
    static_assert(data_type_9_trait::IsDynamic);
    static_assert(data_type_9_trait::Dimension == 5);
    static_assert(!data_type_9_trait::IsDimensionDynamic<0>());
    static_assert(!data_type_9_trait::IsDimensionDynamic<1>());
    static_assert(data_type_9_trait::IsDimensionDynamic<2>());
    static_assert(!data_type_9_trait::IsDimensionDynamic<3>());
    static_assert(!data_type_9_trait::IsDimensionDynamic<4>());
    data_type_9 test_data_type_9;
    KRATOS_EXPECT_EQ(data_type_9_trait::Size(test_data_type_9), 0);
    std::vector<unsigned int> data_type_9_shape{2, 0, 0, 0, 4};
    KRATOS_EXPECT_EQ(data_type_9_trait::Shape(test_data_type_9), data_type_9_shape);
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsNestedDenseMatrixSizeWithShape, KratosCoreFastSuite)
{
    using data_type_traits = DataTypeTraits<std::vector<DenseMatrix<array_1d<std::vector<DenseVector<array_1d<double, 3>>>, 10>>>>;

    std::vector<unsigned int> shape(7);
    shape[0] = 10;
    shape[1] = 4;
    shape[2] = 8;
    shape[3] = 12;
    shape[4] = 20;
    shape[5] = 4;
    shape[6] = 7;

    KRATOS_EXPECT_EQ(data_type_traits::Dimension, shape.size());
    KRATOS_EXPECT_EQ(data_type_traits::Size(shape.begin(), shape.end()), 2150400);
    KRATOS_EXPECT_FALSE(data_type_traits::IsValidShape(shape.data(), shape.data() + shape.size()));

    shape[3] = 10;
    shape[6] = 3;
    KRATOS_EXPECT_TRUE(data_type_traits::IsValidShape(shape.data(), shape.data() + shape.size()));
    KRATOS_EXPECT_EQ(data_type_traits::Size(shape.begin(), shape.end()), 768000);
    KRATOS_EXPECT_FALSE(data_type_traits::IsValidShape(shape.data(), shape.data() + shape.size() - 1));
    KRATOS_EXPECT_FALSE(data_type_traits::IsValidShape(shape.data() + 1, shape.data() + shape.size()));

    shape[3] = 5;
    shape[6] = 2;
    KRATOS_EXPECT_TRUE(data_type_traits::IsValidShape(shape.data(), shape.data() + shape.size()));
    KRATOS_EXPECT_EQ(data_type_traits::Size(shape.begin(), shape.end()), 256000);
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsNestedDenseMatrixCopyToContiguousDataWithShape, KratosCoreFastSuite)
{
    using data_type_traits = DataTypeTraits<std::vector<DenseMatrix<array_1d<std::vector<DenseVector<array_1d<int, 3>>>, 10>>>>;

    std::vector<unsigned int> shape(7);
    shape[0] = 3;       // vector dim
    shape[1] = 2;       // matrix 1dim
    shape[2] = 2;       // matrix 2dim
    shape[3] = 5;       // array_1d dim
    shape[4] = 2;       // vector dim
    shape[5] = 3;       // dense vector dim
    shape[6] = 2;       // array_1d dim

    KRATOS_EXPECT_EQ(data_type_traits::Dimension, shape.size());
    KRATOS_EXPECT_TRUE(data_type_traits::IsValidShape(shape.data(), shape.data() + shape.size()));
    KRATOS_EXPECT_EQ(data_type_traits::Size(shape.begin(), shape.end()), 720);

    // generate the data
    unsigned int value = 0;
    using vec_type_1 = data_type_traits::ContainerType;
    vec_type_1 data;
    for (unsigned int vec_i = 0; vec_i < 5; ++vec_i) {
        using mat_type_1 = DataTypeTraits<vec_type_1>::ValueType;
        mat_type_1 mat_1(3, 4);
        for (unsigned int mat_i = 0; mat_i < 3; ++mat_i) {
            for (unsigned int mat_j = 0; mat_j < 4; ++mat_j) {
                using arr_type_1 = DataTypeTraits<mat_type_1>::ValueType;
                arr_type_1 arr_1;
                for (unsigned int arr_i = 0; arr_i < 10; ++arr_i) {
                    using vec_type_2 = DataTypeTraits<arr_type_1>::ValueType;
                    vec_type_2 vec_2;
                    for (unsigned int vec_j = 0; vec_j < 5; ++vec_j) {
                        using dvec_type_1 = DataTypeTraits<vec_type_2>::ValueType;
                        dvec_type_1 dvec(4);
                        for (unsigned int dvec_i = 0; dvec_i < 4; ++dvec_i) {
                            using arr_type_2 = DataTypeTraits<dvec_type_1>::ValueType;
                            arr_type_2 arr_2;
                            for (unsigned int arr_j = 0; arr_j < 3; ++arr_j) {
                                arr_2[arr_j] = ++value;
                            }
                            dvec[dvec_i] = arr_2;
                        }
                        vec_2.push_back(dvec);
                    }
                    arr_1[arr_i] = vec_2;
                }
                mat_1(mat_i, mat_j) = arr_1;
            }
        }
        data.push_back(mat_1);
    }

    // first copy all the data
    std::vector<int> all_data(data_type_traits::Size(data));
    data_type_traits::CopyToContiguousData(all_data.data(), data);
    for (unsigned int i = 0; i < all_data.size(); ++i) {
        KRATOS_EXPECT_EQ(all_data[i], i + 1);
    }

    const auto& r_origin_shape = data_type_traits::Shape(data);

    std::vector<int> copied_data(720);
    data_type_traits::CopyToContiguousData(copied_data.data(), data, shape.begin(), shape.end());

    std::vector<unsigned int> copy_dim_lengths;
    for (unsigned int i = 0; i < shape.size() - 1; ++i) {
        copy_dim_lengths.push_back(std::accumulate(shape.begin() + i + 1, shape.end(), 1, std::multiplies<int>{}));
    }

    std::vector<unsigned int> orig_dim_lengths;
    for (unsigned int i = 0; i < r_origin_shape.size() - 1; ++i) {
        orig_dim_lengths.push_back(std::accumulate(r_origin_shape.begin() + i + 1, r_origin_shape.end(), 1, std::multiplies<int>{}));
    }

    std::vector<unsigned int> current_copy_dim_index(copy_dim_lengths.size());
    for (unsigned int i = 0; i < 720; ++i) {
        unsigned int index_remainder = i;
        for (unsigned int j = 0; j < copy_dim_lengths.size(); ++j) {
            current_copy_dim_index[j] = index_remainder / copy_dim_lengths[j];
            index_remainder -= current_copy_dim_index[j] * copy_dim_lengths[j];
        }

        unsigned int orig_i = 0;
        for (unsigned int j = 0; j < orig_dim_lengths.size(); ++j) {
            orig_i += current_copy_dim_index[j] * orig_dim_lengths[j];
        }

        KRATOS_EXPECT_EQ(copied_data[i], all_data[orig_i + index_remainder]);
    }

}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsDenseMatrixCopyToFromContiguousDataWithShape, KratosCoreFastSuite)
{
    DenseMatrix<double> input(5, 2);
    for (unsigned int i = 0; i < 10; ++i) {
        *(input.data().begin() + i) = i;
    }

    std::vector<unsigned int> shape = {5, 2};
    std::vector<double> contiguous_data(10);
    DataTypeTraits<DenseMatrix<double>>::CopyToContiguousData(contiguous_data.data(), input, shape.begin(), shape.end());

    for (unsigned int i = 0; i < 10; ++i) {
        KRATOS_EXPECT_EQ(contiguous_data[i], i);
    }

    DenseMatrix<double> output(5, 2);
    DataTypeTraits<DenseMatrix<double>>::CopyFromContiguousData(output, contiguous_data.data(), shape.begin(), shape.end());

    for (unsigned int i = 0; i < 10; ++i) {
        KRATOS_EXPECT_EQ(*(input.data().begin() + i), *(output.data().begin() + i));
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataTypeTraitsNestedDenseMatrixCopyFromContiguousDataWithShape, KratosCoreFastSuite)
{
    using data_type_traits = DataTypeTraits<std::vector<DenseMatrix<array_1d<std::vector<DenseVector<array_1d<int, 3>>>, 10>>>>;

    std::vector<unsigned int> shape(7);
    shape[0] = 3;       // vector dim
    shape[1] = 2;       // matrix 1dim
    shape[2] = 2;       // matrix 2dim
    shape[3] = 5;       // array_1d dim
    shape[4] = 2;       // vector dim
    shape[5] = 3;       // dense vector dim
    shape[6] = 2;       // array_1d dim

    KRATOS_EXPECT_EQ(data_type_traits::Dimension, shape.size());
    KRATOS_EXPECT_TRUE(data_type_traits::IsValidShape(shape.data(), shape.data() + shape.size()));
    KRATOS_EXPECT_EQ(data_type_traits::Size(shape.begin(), shape.end()), 720);

    // generate the data
    unsigned int value = 0;
    using vec_type_1 = data_type_traits::ContainerType;
    vec_type_1 data;
    for (unsigned int vec_i = 0; vec_i < 5; ++vec_i) {
        using mat_type_1 = DataTypeTraits<vec_type_1>::ValueType;
        mat_type_1 mat_1(3, 4);
        for (unsigned int mat_i = 0; mat_i < 3; ++mat_i) {
            for (unsigned int mat_j = 0; mat_j < 4; ++mat_j) {
                using arr_type_1 = DataTypeTraits<mat_type_1>::ValueType;
                arr_type_1 arr_1;
                for (unsigned int arr_i = 0; arr_i < 10; ++arr_i) {
                    using vec_type_2 = DataTypeTraits<arr_type_1>::ValueType;
                    vec_type_2 vec_2;
                    for (unsigned int vec_j = 0; vec_j < 5; ++vec_j) {
                        using dvec_type_1 = DataTypeTraits<vec_type_2>::ValueType;
                        dvec_type_1 dvec(4);
                        for (unsigned int dvec_i = 0; dvec_i < 4; ++dvec_i) {
                            using arr_type_2 = DataTypeTraits<dvec_type_1>::ValueType;
                            arr_type_2 arr_2;
                            for (unsigned int arr_j = 0; arr_j < 3; ++arr_j) {
                                arr_2[arr_j] = ++value;
                            }
                            dvec[dvec_i] = arr_2;
                        }
                        vec_2.push_back(dvec);
                    }
                    arr_1[arr_i] = vec_2;
                }
                mat_1(mat_i, mat_j) = arr_1;
            }
        }
        data.push_back(mat_1);
    }

    // first copy all the data
    std::vector<int> all_data(data_type_traits::Size(shape.begin(), shape.end()));
    data_type_traits::CopyToContiguousData(all_data.data(), data, shape.begin(), shape.end());

    data_type_traits::ContainerType check_data;
    data_type_traits::Reshape(check_data, shape);
    data_type_traits::CopyFromContiguousData(check_data, all_data.data(), shape.begin(), shape.end());

    std::vector<int> check_flat_data(data_type_traits::Size(shape.begin(), shape.end()));
    data_type_traits::CopyToContiguousData(check_flat_data.data(), check_data, shape.begin(), shape.end());

    KRATOS_CHECK_VECTOR_EQUAL(check_flat_data, all_data);

}

} // namespace Kratos::Testing