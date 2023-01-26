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
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "containers/model.h"

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

KRATOS_TEST_CASE_IN_SUITE(CheckAndCorrectZeroDiagonalValues, KratosCoreFastSuite)
{
    /// The sparse matrix type
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 1.0);
    
    SparseMatrixType matrix12x12(12, 12);
    for (IndexType i = 0; i < 12; ++i) {
        matrix12x12.push_back(i, i, static_cast<double>(i));
    }
    Vector vector12 = ZeroVector(12);

    // Test the norm of the matrix
    double norm = 0.0;
    norm = SparseSpaceType::CheckAndCorrectZeroDiagonalValues(r_process_info, matrix12x12, vector12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
    KRATOS_CHECK_DOUBLE_EQUAL(matrix12x12(0, 0), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(GetScaleNorm, KratosCoreFastSuite)
{
    /// The sparse matrix type
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 3.0);
    
    SparseMatrixType matrix12x12(12, 12);
    for (IndexType i = 0; i < 12; ++i) {
        matrix12x12.push_back(i, i, static_cast<double>(i+1));
    }

    // Test the norm of the matrix
    double norm = 0.0;
    norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
    norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 3.0);
    norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL);
    KRATOS_CHECK_NEAR(norm, 2.124591464, 1.0e-6);
    norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 12.0);
    norm = SparseSpaceType::GetAveragevalueDiagonal(matrix12x12);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 6.5);
    norm = SparseSpaceType::GetMinDiagonal(matrix12x12);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
}

} // namespace Testing
} // namespace Kratos.
