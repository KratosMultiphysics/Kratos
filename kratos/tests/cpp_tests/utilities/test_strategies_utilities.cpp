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

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "utilities/strategies_utilities.h"

namespace Kratos {
namespace Testing {

/// The sparse matrix type
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef typename SparseSpaceType::MatrixType SparseMatrixType;

KRATOS_TEST_CASE_IN_SUITE(CheckAndCorrectZeroDiagonalValues, KratosCoreFastSuite)
{
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
    norm = StrategiesUtilities::CheckAndCorrectZeroDiagonalValues<SparseSpaceType>(r_process_info, matrix12x12, vector12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
    KRATOS_CHECK_DOUBLE_EQUAL(matrix12x12(0, 0), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(GetScaleNorm, KratosCoreFastSuite)
{
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
    norm = StrategiesUtilities::GetScaleNorm<SparseSpaceType>(r_process_info, matrix12x12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
    norm = StrategiesUtilities::GetScaleNorm<SparseSpaceType>(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 3.0);
    norm = StrategiesUtilities::GetScaleNorm<SparseSpaceType>(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL);
    KRATOS_CHECK_NEAR(norm, 2.124591464, 1.0e-6);
    norm = StrategiesUtilities::GetScaleNorm<SparseSpaceType>(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 12.0);
    norm = StrategiesUtilities::GetAveragevalueDiagonal<SparseSpaceType>(matrix12x12);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 6.5);
    norm = StrategiesUtilities::GetMinDiagonal<SparseSpaceType>(matrix12x12);
    KRATOS_CHECK_DOUBLE_EQUAL(norm, 1.0);
}

}   // namespace Testing
}  // namespace Kratos.
