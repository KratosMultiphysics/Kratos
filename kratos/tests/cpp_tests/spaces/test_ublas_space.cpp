//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "testing/testing.h"
#include "spaces/ublas_space.h"


namespace Kratos {
namespace Testing {


KRATOS_TEST_CASE_IN_SUITE(UblasSpaceNormSparseMatrix, KratosCoreFastSuite)
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;

    const std::size_t size = 10;
    SparseSpaceType::MatrixType mat(size, size);

    for (std::size_t i=0; i<mat.size1(); ++i)	{
        if (i>=1) {mat.push_back(i, i-1, -1.123);}
        mat.push_back(i, i, 4.5);
        if (i+1<mat.size2()) {mat.push_back(i, i+1, 2.336);}
    }

    KRATOS_CHECK_NEAR(16.216110045260546, SparseSpaceType::TwoNorm(mat), 1e-12);
    KRATOS_CHECK_NEAR(31.131, SparseSpaceType::JacobiNorm(mat), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(UblasSpaceNormDenseMatrix, KratosCoreFastSuite)
{
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    const std::size_t size = 10;
    LocalSpaceType::MatrixType mat(size, size, 0.0);

    for (std::size_t i=0; i<mat.size1(); ++i)	{
        mat(i,i) = 4.5;

        if (i>=1) {mat(i,i-1) = -1.123;}

        if (i+1<mat.size2()) {mat(i,i+1) = 2.336;}
    }

    KRATOS_CHECK_NEAR(16.216110045260546, LocalSpaceType::TwoNorm(mat), 1e-12);
    KRATOS_CHECK_NEAR(31.131, LocalSpaceType::JacobiNorm(mat), 1e-12);
}

} // namespace Testing
} // namespace Kratos.
