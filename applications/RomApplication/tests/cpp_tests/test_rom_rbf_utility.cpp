//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

// System includes
#include <cmath>

/* External includes */

/* Project includes */
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "custom_utilities/rom_rbf_utility.h"

namespace Kratos::Testing {
namespace RomRBFUtilityTestingInternal {

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using RomRBFUtilityType = RomRBFUtility<SparseSpaceType, LocalSpaceType>;

constexpr std::size_t NumberOfDofs = 6;
constexpr std::size_t NumberOfInfModes = 2;
constexpr std::size_t NumberOfSupModes = 3;
constexpr std::size_t NumberOfCenters = 4;

Matrix GetMatrix(const std::size_t Size1, const std::size_t Size2, const double Seed)
{
    Matrix matrix(Size1, Size2);
    for (std::size_t i = 0; i < Size1; ++i) {
        for (std::size_t j = 0; j < Size2; ++j) {
            matrix(i, j) = std::sin(Seed + 1.3*i + 0.7*j);
        }
    }
    return matrix;
}

/// Checks the decoder, its gradient (against finite differences) and that the centers are left untouched
void CheckDecoder(const IndexType KernelType)
{
    vector<Matrix> svd_phi_matrices(3);
    svd_phi_matrices[0] = GetMatrix(NumberOfDofs, NumberOfInfModes, 0.1);
    svd_phi_matrices[1] = GetMatrix(NumberOfDofs, NumberOfSupModes, 0.5);
    svd_phi_matrices[2] = ZeroMatrix(NumberOfInfModes, NumberOfInfModes);
    svd_phi_matrices[2](0, 0) = 0.5;
    svd_phi_matrices[2](1, 1) = 2.0;
    const Matrix weights = GetMatrix(NumberOfCenters, NumberOfSupModes, 0.9);
    const Matrix centers = GetMatrix(NumberOfCenters, NumberOfInfModes, 1.7);
    const Matrix original_centers = centers;
    const double kernel_eps = 0.8;
    Vector ref_snapshot(NumberOfDofs);
    for (std::size_t i = 0; i < NumberOfDofs; ++i) ref_snapshot[i] = 0.1*i;
    Vector rom_unknowns(NumberOfInfModes);
    rom_unknowns[0] = 0.3;
    rom_unknowns[1] = -0.2;

    Vector x(NumberOfDofs);
    Matrix phi_global(NumberOfDofs, NumberOfInfModes);
    RomRBFUtilityType::GetXAndDecoderGradient(rom_unknowns, x, phi_global, svd_phi_matrices, weights, centers, KernelType, kernel_eps, ref_snapshot);

    // The decoder and its gradient give the same results when called separately
    Vector x_decoder(NumberOfDofs);
    Matrix phi_global_gradient(NumberOfDofs, NumberOfInfModes);
    RomRBFUtilityType::GetXFromDecoder(rom_unknowns, x_decoder, svd_phi_matrices, weights, centers, KernelType, kernel_eps, ref_snapshot);
    RomRBFUtilityType::GetDecoderGradient(rom_unknowns, phi_global_gradient, svd_phi_matrices, weights, centers, KernelType, kernel_eps);
    KRATOS_EXPECT_VECTOR_NEAR(x, x_decoder, 1.0e-14);
    KRATOS_EXPECT_MATRIX_NEAR(phi_global, phi_global_gradient, 1.0e-14);

    // The gradient matches the central finite differences of the decoder
    const double delta = 1.0e-6;
    for (std::size_t j = 0; j < NumberOfInfModes; ++j) {
        Vector rom_unknowns_plus = rom_unknowns;
        Vector rom_unknowns_minus = rom_unknowns;
        rom_unknowns_plus[j] += delta;
        rom_unknowns_minus[j] -= delta;
        Vector x_plus(NumberOfDofs);
        Vector x_minus(NumberOfDofs);
        RomRBFUtilityType::GetXFromDecoder(rom_unknowns_plus, x_plus, svd_phi_matrices, weights, centers, KernelType, kernel_eps, ref_snapshot);
        RomRBFUtilityType::GetXFromDecoder(rom_unknowns_minus, x_minus, svd_phi_matrices, weights, centers, KernelType, kernel_eps, ref_snapshot);
        for (std::size_t i = 0; i < NumberOfDofs; ++i) {
            KRATOS_EXPECT_NEAR(phi_global(i, j), (x_plus[i] - x_minus[i])/(2.0*delta), 1.0e-8);
        }
    }

    // At a center the kernel of that center is one, so with a single center the decoder is known
    const Matrix single_center = GetMatrix(1, NumberOfInfModes, 1.7);
    const Matrix single_weights = GetMatrix(1, NumberOfSupModes, 0.9);
    Vector rom_unknowns_at_center(NumberOfInfModes);
    rom_unknowns_at_center[0] = single_center(0, 0)/svd_phi_matrices[2](0, 0);
    rom_unknowns_at_center[1] = single_center(0, 1)/svd_phi_matrices[2](1, 1);
    Vector x_at_center(NumberOfDofs);
    RomRBFUtilityType::GetXFromDecoder(rom_unknowns_at_center, x_at_center, svd_phi_matrices, single_weights, single_center, KernelType, kernel_eps, ref_snapshot);
    for (std::size_t i = 0; i < NumberOfDofs; ++i) {
        double expected = ref_snapshot[i];
        for (std::size_t j = 0; j < NumberOfInfModes; ++j) expected += svd_phi_matrices[0](i, j)*rom_unknowns_at_center[j];
        for (std::size_t j = 0; j < NumberOfSupModes; ++j) expected += svd_phi_matrices[1](i, j)*single_weights(0, j);
        KRATOS_EXPECT_NEAR(x_at_center[i], expected, 1.0e-12);
    }

    // The centers are not modified by the calls
    KRATOS_EXPECT_MATRIX_NEAR(centers, original_centers, 0.0);
}

} // namespace RomRBFUtilityTestingInternal

KRATOS_TEST_CASE_IN_SUITE(RomRBFUtilityGaussianKernel, RomApplicationFastSuite)
{
    RomRBFUtilityTestingInternal::CheckDecoder(0);
}

KRATOS_TEST_CASE_IN_SUITE(RomRBFUtilityInverseMultiquadricKernel, RomApplicationFastSuite)
{
    RomRBFUtilityTestingInternal::CheckDecoder(1);
}

KRATOS_TEST_CASE_IN_SUITE(RomRBFUtilityUnknownKernel, RomApplicationFastSuite)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(RomRBFUtilityTestingInternal::CheckDecoder(2), "Unknown RBF kernel type 2");
}

} // namespace Kratos::Testing
