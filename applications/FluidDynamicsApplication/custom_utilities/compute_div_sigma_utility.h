//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "utilities/dense_householder_qr_decomposition.h"

namespace Kratos
{

class ComputeDivSigmaUtility
{
public:
    using DenseSpace = UblasSpace<double, Matrix, Vector>;

    ComputeDivSigmaUtility() = default;
    virtual ~ComputeDivSigmaUtility() = default;

    Matrix ComputeDivergence(
        const Matrix& rSigmaValues,
        const Matrix& rShapeFunctionValues,
        const Matrix& rShapeFunctionValuesDx,
        const Matrix& rShapeFunctionValuesDy) const
    {
        return ComputeDivergenceInternal(
            rSigmaValues,
            rShapeFunctionValues,
            rShapeFunctionValuesDx,
            rShapeFunctionValuesDy,
            nullptr);
    }

    Matrix ComputeDivergence(
        const Matrix& rSigmaValues,
        const Matrix& rShapeFunctionValues,
        const Matrix& rShapeFunctionValuesDx,
        const Matrix& rShapeFunctionValuesDy,
        const Matrix& rShapeFunctionValuesDz) const
    {
        return ComputeDivergenceInternal(
            rSigmaValues,
            rShapeFunctionValues,
            rShapeFunctionValuesDx,
            rShapeFunctionValuesDy,
            &rShapeFunctionValuesDz);
    }

private:
    Matrix ComputeDivergenceInternal(
        const Matrix& rSigmaValues,
        const Matrix& rShapeFunctionValues,
        const Matrix& rShapeFunctionValuesDx,
        const Matrix& rShapeFunctionValuesDy,
        const Matrix* pShapeFunctionValuesDz) const
    {
        const Matrix coefficients = ComputeCoefficients(
            rSigmaValues,
            rShapeFunctionValues,
            rShapeFunctionValuesDx,
            rShapeFunctionValuesDy,
            pShapeFunctionValuesDz);

        const std::size_t gauss_point_number = rSigmaValues.size1();
        const bool is_3d = pShapeFunctionValuesDz != nullptr;
        Matrix divergence_results(gauss_point_number, is_3d ? 3 : 2);

        for (std::size_t i = 0; i < gauss_point_number; ++i) {
            double d_sigma_xx_dx = 0.0;
            double d_sigma_yy_dy = 0.0;
            double d_sigma_xy_dx = 0.0;
            double d_sigma_xy_dy = 0.0;
            double d_sigma_xz_dz = 0.0;
            double d_sigma_xz_dx = 0.0;
            double d_sigma_yz_dy = 0.0;
            double d_sigma_yz_dz = 0.0;
            double d_sigma_zz_dz = 0.0;

            for (std::size_t j = 0; j < gauss_point_number; ++j) {
                d_sigma_xx_dx += coefficients(0, j) * rShapeFunctionValuesDx(i, j);
                d_sigma_yy_dy += coefficients(1, j) * rShapeFunctionValuesDy(i, j);
                d_sigma_xy_dx += coefficients(2, j) * rShapeFunctionValuesDx(i, j);
                d_sigma_xy_dy += coefficients(2, j) * rShapeFunctionValuesDy(i, j);
                if (is_3d) {
                    d_sigma_zz_dz += coefficients(3, j) * (*pShapeFunctionValuesDz)(i, j);
                    d_sigma_yz_dy += coefficients(4, j) * rShapeFunctionValuesDy(i, j);
                    d_sigma_yz_dz += coefficients(4, j) * (*pShapeFunctionValuesDz)(i, j);
                    d_sigma_xz_dx += coefficients(5, j) * rShapeFunctionValuesDx(i, j);
                    d_sigma_xz_dz += coefficients(5, j) * (*pShapeFunctionValuesDz)(i, j);
                }
            }

            double div_sigma_x = d_sigma_xx_dx + d_sigma_xy_dy;
            double div_sigma_y = d_sigma_xy_dx + d_sigma_yy_dy;
            if (is_3d) {
                const double div_sigma_z = d_sigma_xz_dx + d_sigma_yz_dy + d_sigma_zz_dz;
                div_sigma_x += d_sigma_xz_dz;
                div_sigma_y += d_sigma_yz_dz;
                divergence_results(i, 0) = div_sigma_x;
                divergence_results(i, 1) = div_sigma_y;
                divergence_results(i, 2) = div_sigma_z;
            } else {
                divergence_results(i, 0) = div_sigma_x;
                divergence_results(i, 1) = div_sigma_y;
            }
        }

        return divergence_results;
    }

    Matrix ComputeCoefficients(
        const Matrix& rSigmaValues,
        const Matrix& rShapeFunctionValues,
        const Matrix& rShapeFunctionValuesDx,
        const Matrix& rShapeFunctionValuesDy,
        const Matrix* pShapeFunctionValuesDz) const
    {
        const std::size_t gauss_point_number = rSigmaValues.size1();
        KRATOS_ERROR_IF(gauss_point_number != rShapeFunctionValues.size1())
            << "Mismatch between number of sigma values and N!" << std::endl;
        KRATOS_ERROR_IF(rShapeFunctionValues.size2() != gauss_point_number)
            << "Reconstruction matrix must be square. Got " << rShapeFunctionValues.size1()
            << "x" << rShapeFunctionValues.size2() << "." << std::endl;
        KRATOS_ERROR_IF(rShapeFunctionValuesDx.size1() != rShapeFunctionValues.size1() ||
                        rShapeFunctionValuesDx.size2() != rShapeFunctionValues.size2())
            << "Mismatch between Dx values and N." << std::endl;
        KRATOS_ERROR_IF(rShapeFunctionValuesDy.size1() != rShapeFunctionValues.size1() ||
                        rShapeFunctionValuesDy.size2() != rShapeFunctionValues.size2())
            << "Mismatch between Dy values and N." << std::endl;

        const std::size_t sigma_size = rSigmaValues.size2();
        const bool is_3d = (sigma_size == 6);
        KRATOS_ERROR_IF(!(sigma_size == 3 || sigma_size == 6))
            << "Unsupported sigma size " << sigma_size << ". Expected 3 (2D) or 6 (3D)." << std::endl;
        KRATOS_ERROR_IF(is_3d != (pShapeFunctionValuesDz != nullptr))
            << "Dz shape-function data must be provided if and only if sigma size is 6." << std::endl;
        KRATOS_ERROR_IF(is_3d &&
                        (pShapeFunctionValuesDz->size1() != rShapeFunctionValues.size1() ||
                         pShapeFunctionValuesDz->size2() != rShapeFunctionValues.size2()))
            << "Mismatch between number of Dz values and N for 3D." << std::endl;

        Matrix qr_matrix = rShapeFunctionValues;
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(qr_matrix);

        Matrix r_matrix;
        qr_decomposition.MatrixR(r_matrix);

        double max_abs_r_diag = 0.0;
        double min_abs_r_diag = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < gauss_point_number; ++i) {
            const double abs_r_ii = std::abs(r_matrix(i, i));
            max_abs_r_diag = std::max(max_abs_r_diag, abs_r_ii);
            min_abs_r_diag = std::min(min_abs_r_diag, abs_r_ii);
        }

        KRATOS_ERROR_IF(max_abs_r_diag <= std::numeric_limits<double>::epsilon())
            << "Degenerate shape-function reconstruction matrix detected in ComputeDivSigmaUtility."
            << std::endl;

        const double rank_tolerance =
            1.0e3 * std::numeric_limits<double>::epsilon()
            * static_cast<double>(gauss_point_number) * max_abs_r_diag;

        KRATOS_ERROR_IF(min_abs_r_diag <= rank_tolerance)
            << "Near-singular shape-function reconstruction matrix detected in ComputeDivSigmaUtility. "
            << "min(|R_ii|) = " << min_abs_r_diag
            << ", tolerance = " << rank_tolerance
            << ", max(|R_ii|) = " << max_abs_r_diag << std::endl;

        Vector b_xx(gauss_point_number), b_yy(gauss_point_number), b_xy(gauss_point_number);
        for (std::size_t i = 0; i < gauss_point_number; ++i) {
            b_xx(i) = rSigmaValues(i, 0);
            b_yy(i) = rSigmaValues(i, 1);
            b_xy(i) = is_3d ? rSigmaValues(i, 3) : rSigmaValues(i, 2);
        }

        Vector c_xx(gauss_point_number), c_yy(gauss_point_number), c_xy(gauss_point_number);
        qr_decomposition.Solve(b_xx, c_xx);
        qr_decomposition.Solve(b_yy, c_yy);
        qr_decomposition.Solve(b_xy, c_xy);

        Vector c_zz, c_yz, c_xz;
        Vector b_zz, b_yz, b_xz;
        if (is_3d) {
            c_zz.resize(gauss_point_number, false);
            c_yz.resize(gauss_point_number, false);
            c_xz.resize(gauss_point_number, false);
            b_zz.resize(gauss_point_number, false);
            b_yz.resize(gauss_point_number, false);
            b_xz.resize(gauss_point_number, false);
            for (std::size_t i = 0; i < gauss_point_number; ++i) {
                b_zz(i) = rSigmaValues(i, 2);
                b_yz(i) = rSigmaValues(i, 4);
                b_xz(i) = rSigmaValues(i, 5);
            }
            qr_decomposition.Solve(b_zz, c_zz);
            qr_decomposition.Solve(b_yz, c_yz);
            qr_decomposition.Solve(b_xz, c_xz);
        }

#ifdef KRATOS_DEBUG
        CheckResidualNorm(rShapeFunctionValues, b_xx, c_xx, "sigma_xx");
        CheckResidualNorm(rShapeFunctionValues, b_yy, c_yy, "sigma_yy");
        CheckResidualNorm(rShapeFunctionValues, b_xy, c_xy, "sigma_xy");
        if (is_3d) {
            CheckResidualNorm(rShapeFunctionValues, b_zz, c_zz, "sigma_zz");
            CheckResidualNorm(rShapeFunctionValues, b_yz, c_yz, "sigma_yz");
            CheckResidualNorm(rShapeFunctionValues, b_xz, c_xz, "sigma_xz");
        }
#endif

        Matrix coefficients(is_3d ? 6 : 3, gauss_point_number);
        for (std::size_t j = 0; j < gauss_point_number; ++j) {
            coefficients(0, j) = c_xx(j);
            coefficients(1, j) = c_yy(j);
            coefficients(2, j) = c_xy(j);
            if (is_3d) {
                coefficients(3, j) = c_zz(j);
                coefficients(4, j) = c_yz(j);
                coefficients(5, j) = c_xz(j);
            }
        }

        return coefficients;
    }

    void CheckResidualNorm(
        const Matrix& rA,
        const Vector& rB,
        const Vector& rX,
        const char* pLabel) const
    {
        const Vector residual = rB - prod(rA, rX);
        const double norm_residual = norm_2(residual);
        constexpr double threshold = 1e-9;
        if (norm_residual > threshold) {
            std::cerr << "Warning: High residual norm detected for " << pLabel << " fit!" << std::endl;
            std::cerr << "Residual norm (" << pLabel << "): " << norm_residual << std::endl;
        }
    }
};

} // namespace Kratos
