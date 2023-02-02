//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifdef KRATOS_USE_AMATRIX
#ifdef KRATOS_DEBUG

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "includes/amatrix_interface.h"

namespace Kratos {
namespace Testing {

TEST(TestAMatrixAliasProd, KratosCoreFastSuite)
{
    Matrix A = ScalarMatrix(2,2,1.0);
    Matrix B = ScalarMatrix(2,2,2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(A) = prod(A,B), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(A,B), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(B,A), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = A * B, "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = B * A, "Aliasing found in assigning Matrix" );
}

TEST(TestAMatrixAliasTranspose, KratosCoreFastSuite)
{
    Matrix A = ScalarMatrix(2,2,1.0);
    Matrix B = ScalarMatrix(2,2,2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(A) = trans(A), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = trans(A), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(B,trans(A)), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = trans(A) * B, "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = B * trans(A), "Aliasing found in assigning Matrix" );
}

TEST(TestAMatrixAliasProdBounded, KratosCoreFastSuite)
{
    BoundedMatrix<double,2,2> A = ScalarMatrix(2,2,1.0);
    BoundedMatrix<double,2,2> B = ScalarMatrix(2,2,2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(A) = prod(A,B), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(A,B), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(B,A), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = A * B, "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = B * A, "Aliasing found in assigning Matrix" );
}

TEST(TestAMatrixAliasTransposeBounded, KratosCoreFastSuite)
{
    BoundedMatrix<double,2,2> A = ScalarMatrix(2,2,1.0);
    BoundedMatrix<double,2,2> B = ScalarMatrix(2,2,2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(A) = trans(A), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = trans(A), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(B,trans(A)), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = trans(A) * B, "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = B * trans(A), "Aliasing found in assigning Matrix" );
}


TEST(TestAMatrixAliasProdVector, KratosCoreFastSuite)
{
    Matrix A = ScalarMatrix(2,2,1.0);
    Vector b = ScalarMatrix(2,1,2.0);
    Vector c = ScalarMatrix(2,1,2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(b) = prod(A,b), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( b = prod(A,b), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( b = A * b, "Aliasing found in assigning Matrix" );
    c = A * b; // This one has no alias, should not throw
}

TEST(TestAMatrixAliasProdVectorBounded, KratosCoreFastSuite)
{
    BoundedMatrix<double,2,2> A = ScalarMatrix(2,2,1.0);
    BoundedVector<double,2> b = ScalarMatrix(2,1,2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(b) = prod(A,b), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( b = prod(A,b), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( b = A * b, "Aliasing found in assigning Matrix" );
}

TEST(TestAMatrixAliasProdProd, KratosCoreFastSuite)
{
    Matrix A = ScalarMatrix(2,2,1.0);
    Matrix B = ScalarMatrix(2,2,2.0);
    Matrix C = ScalarMatrix(2,2,3.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(A) = prod(B,prod(C,A)), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(B,prod(C,A)), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = B * C * A, "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = prod(B,prod(C,trans(A))), "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = B * C * trans(A), "Aliasing found in assigning Matrix" );
}

TEST(TestAMatrixAliasSum, KratosCoreFastSuite)
{
    Matrix A = ScalarMatrix(2,2,1.0);
    Matrix B = ScalarMatrix(2,2,2.0);

    // These should not throw
    A = A + B;
    A += A;
    A -= A;
}

TEST(TestAMatrixAliasSumTranspose, KratosCoreFastSuite)
{
    Matrix A = ScalarMatrix(2,2,1.0);
    Matrix B = ScalarMatrix(2,2,2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN( noalias(A) = trans(A) + B, "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A = trans(A) + B, "Aliasing found in assigning Matrix" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A += trans(A), "Aliasing found in += operator" );
    KRATOS_EXPECT_EXCEPTION_IS_THROWN( A -= trans(A), "Aliasing found in -= operator" );
}

}   // namespace Testing.
}  // namespace Kratos.

#endif // KRATOS_DEBUG
#endif // KRATOS_USE_AMATRIX
