//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

/* System includes */
#include <vector>
// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/SVD>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "geometries/geometry.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
using VectorType = Vector;
using CoordinateType = array_1d<double, 3>;
using CoordinateArrayType = std::vector<CoordinateType>;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


class ComputeDivSigmaUtility
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ComputeDivSigmaUtility() = default;

    /// Destructor
    virtual ~ComputeDivSigmaUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Sets the input data for the utility
     * @param rSigmaValues A vector of stress vectors (sigma)
     * @param rCoordinates A vector of coordinates corresponding to each sigma
     */

    void SetInputData(
        const std::vector<std::vector<double>>& rSigmaValues,
        const CoordinateArrayType& rCoordinates,
        const std::vector<std::vector<double>>& rShapeFunctionValues,
        const std::vector<std::vector<double>>& rShapeFunctionValuesDx,
        const std::vector<std::vector<double>>& rShapeFunctionValuesDy)
    {
        mSigmaValues = rSigmaValues;
        mCoordinates = rCoordinates;
        mShapeFunctionValues = rShapeFunctionValues;
        mShapeFunctionValuesDx = rShapeFunctionValuesDx;
        mShapeFunctionValuesDy = rShapeFunctionValuesDy;
        mShapeFunctionValuesDz.clear();
    }

    void SetInputData(
        const std::vector<std::vector<double>>& rSigmaValues,
        const CoordinateArrayType& rCoordinates,
        const std::vector<std::vector<double>>& rShapeFunctionValues,
        const std::vector<std::vector<double>>& rShapeFunctionValuesDx,
        const std::vector<std::vector<double>>& rShapeFunctionValuesDy,
        const std::vector<std::vector<double>>& rShapeFunctionValuesDz)
    {
        mSigmaValues = rSigmaValues;
        mCoordinates = rCoordinates;
        mShapeFunctionValues = rShapeFunctionValues;
        mShapeFunctionValuesDx = rShapeFunctionValuesDx;
        mShapeFunctionValuesDy = rShapeFunctionValuesDy;
        mShapeFunctionValuesDz = rShapeFunctionValuesDz;
    }


    std::vector<std::vector<double>> ComputeDivergence()
    {
        std::vector<std::vector<double>> collected_coefficients = ComputeCoefficients();

        std::vector<std::vector<double>> divergence_results;
        const std::size_t sigma_size = mSigmaValues.empty() ? 0 : mSigmaValues[0].size();
        const bool is_3d = (sigma_size == 6);

        const std::vector<double>& coefficients_xx = collected_coefficients[0];
        const std::vector<double>& coefficients_yy = collected_coefficients[1];
        const std::vector<double>& coefficients_xy = collected_coefficients[2];

        for (size_t i = 0; i < mShapeFunctionValuesDx.size(); ++i)
        {
            double d_sigma_xx_dx = 0.0;
            double d_sigma_yy_dy = 0.0;
            double d_sigma_xy_dx = 0.0;
            double d_sigma_xy_dy = 0.0;
            double d_sigma_xz_dz = 0.0;
            double d_sigma_xz_dx = 0.0;
            double d_sigma_yz_dy = 0.0;
            double d_sigma_yz_dz = 0.0;
            double d_sigma_zz_dz = 0.0;

            const std::vector<double>* p_coefficients_zz = nullptr;
            const std::vector<double>* p_coefficients_yz = nullptr;
            const std::vector<double>* p_coefficients_xz = nullptr;
            if (is_3d) {
                p_coefficients_zz = &collected_coefficients[3];
                p_coefficients_yz = &collected_coefficients[4];
                p_coefficients_xz = &collected_coefficients[5];
            }

            for (size_t j = 0; j < coefficients_xx.size(); ++j)
            {
                d_sigma_xx_dx += coefficients_xx[j] * mShapeFunctionValuesDx[i][j];
                d_sigma_yy_dy += coefficients_yy[j] * mShapeFunctionValuesDy[i][j];
                d_sigma_xy_dx += coefficients_xy[j] * mShapeFunctionValuesDx[i][j];
                d_sigma_xy_dy += coefficients_xy[j] * mShapeFunctionValuesDy[i][j];
                if (is_3d) {
                    d_sigma_zz_dz += (*p_coefficients_zz)[j] * mShapeFunctionValuesDz[i][j];
                    d_sigma_yz_dy += (*p_coefficients_yz)[j] * mShapeFunctionValuesDy[i][j];
                    d_sigma_yz_dz += (*p_coefficients_yz)[j] * mShapeFunctionValuesDz[i][j];
                    d_sigma_xz_dx += (*p_coefficients_xz)[j] * mShapeFunctionValuesDx[i][j];
                    d_sigma_xz_dz += (*p_coefficients_xz)[j] * mShapeFunctionValuesDz[i][j];
                }
            }
            double div_sigma_x = d_sigma_xx_dx + d_sigma_xy_dy;
            double div_sigma_y = d_sigma_xy_dx + d_sigma_yy_dy;
            if (is_3d) {
                double div_sigma_z = d_sigma_xz_dx + d_sigma_yz_dy + d_sigma_zz_dz;
                div_sigma_x += d_sigma_xz_dz;
                div_sigma_y += d_sigma_yz_dz;
                divergence_results.push_back({div_sigma_x, div_sigma_y, div_sigma_z});
            } else {
                divergence_results.push_back({div_sigma_x, div_sigma_y});
            }
        }

        return divergence_results;
    }


    std::vector<std::vector<double>> ComputeCoefficients()
    {
        KRATOS_ERROR_IF(mSigmaValues.size() != mShapeFunctionValues.size()) 
            << "Mismatch between number of sigma values and N!" << std::endl;
        const std::size_t sigma_size = mSigmaValues.empty() ? 0 : mSigmaValues[0].size();
        const bool is_3d = (sigma_size == 6);
        KRATOS_ERROR_IF(!(sigma_size == 3 || sigma_size == 6))
            << "Unsupported sigma size " << sigma_size << ". Expected 3 (2D) or 6 (3D)." << std::endl;
        KRATOS_ERROR_IF(is_3d && mShapeFunctionValuesDz.size() != mShapeFunctionValues.size())
            << "Mismatch between number of Dz values and N for 3D." << std::endl;
        
        std::vector<std::vector<double>> coefficients_results;
        size_t gauss_point_number = mSigmaValues.size();
        // Construct design matrix A and vector b for quadratic fit
        Eigen::MatrixXd A(gauss_point_number, gauss_point_number); // 9 Gauss points, 9 shape function values per Gauss point
        Eigen::VectorXd b_xx(gauss_point_number), b_yy(gauss_point_number), b_xy(gauss_point_number);
        Eigen::VectorXd b_zz(gauss_point_number), b_yz(gauss_point_number), b_xz(gauss_point_number);

        for (size_t i = 0; i < gauss_point_number; ++i)
        {
            // Use the shape function values directly to build matrix A
            const std::vector<double>& shape_function_values = mShapeFunctionValues[i];

            for (size_t j = 0; j < gauss_point_number; ++j)
            {
                A(i, j) = shape_function_values[j];
            }

            // Prepare the three b vectors for sigma_xx, sigma_yy, and sigma_xy
            b_xx(i) = mSigmaValues[i][0];
            b_yy(i) = mSigmaValues[i][1];
            b_xy(i) = mSigmaValues[i][2];
            if (is_3d) {
                b_zz(i) = mSigmaValues[i][2];
                b_xy(i) = mSigmaValues[i][3];
                b_yz(i) = mSigmaValues[i][4];
                b_xz(i) = mSigmaValues[i][5];
            }
        }

        // Solve for coefficients using least squares
        Eigen::VectorXd c_xx = A.colPivHouseholderQr().solve(b_xx);
        Eigen::VectorXd c_yy = A.colPivHouseholderQr().solve(b_yy);
        Eigen::VectorXd c_xy = A.colPivHouseholderQr().solve(b_xy);
        Eigen::VectorXd c_zz, c_yz, c_xz;
        if (is_3d) {
            c_zz = A.colPivHouseholderQr().solve(b_zz);
            c_yz = A.colPivHouseholderQr().solve(b_yz);
            c_xz = A.colPivHouseholderQr().solve(b_xz);
        }

        Eigen::VectorXd residual_xx = b_xx - A * c_xx;
        Eigen::VectorXd residual_yy = b_yy - A * c_yy;
        Eigen::VectorXd residual_xy = b_xy - A * c_xy;
        Eigen::VectorXd residual_zz, residual_yz, residual_xz;
        if (is_3d) {
            residual_zz = b_zz - A * c_zz;
            residual_yz = b_yz - A * c_yz;
            residual_xz = b_xz - A * c_xz;
        }

        // Compute norms of residuals
        double norm_residual_xx = residual_xx.norm();
        double norm_residual_yy = residual_yy.norm();
        double norm_residual_xy = residual_xy.norm();
        double norm_residual_zz = 0.0;
        double norm_residual_yz = 0.0;
        double norm_residual_xz = 0.0;
        if (is_3d) {
            norm_residual_zz = residual_zz.norm();
            norm_residual_yz = residual_yz.norm();
            norm_residual_xz = residual_xz.norm();
        }

        // Define a threshold for acceptable residuals
        const double threshold = 1e-9; // Adjust this value as necessary

        if (norm_residual_xx > threshold) {
            std::cerr << "Warning: High residual norm detected for sigma_xx fit!" << std::endl;
            std::cerr << "Residual norm (sigma_xx): " << norm_residual_xx << std::endl;
            std::cerr << "Residual vector: " << residual_xx.transpose() << std::endl;
        }

        if (norm_residual_yy > threshold) {
            std::cerr << "Warning: High residual norm detected for sigma_yy fit!" << std::endl;
            std::cerr << "Residual norm (sigma_yy): " << norm_residual_yy << std::endl;
            std::cerr << "Residual vector: " << residual_yy.transpose() << std::endl;
        }

        if (norm_residual_xy > threshold) {
            std::cerr << "Warning: High residual norm detected for sigma_xy fit!" << std::endl;
            std::cerr << "Residual norm (sigma_xy): " << norm_residual_xy << std::endl;
            std::cerr << "Residual vector: " << residual_xy.transpose() << std::endl;
        }
        if (is_3d && norm_residual_zz > threshold) {
            std::cerr << "Warning: High residual norm detected for sigma_zz fit!" << std::endl;
            std::cerr << "Residual norm (sigma_zz): " << norm_residual_zz << std::endl;
            std::cerr << "Residual vector: " << residual_zz.transpose() << std::endl;
        }
        if (is_3d && norm_residual_yz > threshold) {
            std::cerr << "Warning: High residual norm detected for sigma_yz fit!" << std::endl;
            std::cerr << "Residual norm (sigma_yz): " << norm_residual_yz << std::endl;
            std::cerr << "Residual vector: " << residual_yz.transpose() << std::endl;
        }
        if (is_3d && norm_residual_xz > threshold) {
            std::cerr << "Warning: High residual norm detected for sigma_xz fit!" << std::endl;
            std::cerr << "Residual norm (sigma_xz): " << norm_residual_xz << std::endl;
            std::cerr << "Residual vector: " << residual_xz.transpose() << std::endl;
        }

        // Convert Eigen vectors to std::vector and store results
        std::vector<double> coefficients_xx(c_xx.data(), c_xx.data() + c_xx.size());
        std::vector<double> coefficients_yy(c_yy.data(), c_yy.data() + c_yy.size());
        std::vector<double> coefficients_xy(c_xy.data(), c_xy.data() + c_xy.size());
        std::vector<double> coefficients_zz;
        std::vector<double> coefficients_yz;
        std::vector<double> coefficients_xz;
        if (is_3d) {
            coefficients_zz.assign(c_zz.data(), c_zz.data() + c_zz.size());
            coefficients_yz.assign(c_yz.data(), c_yz.data() + c_yz.size());
            coefficients_xz.assign(c_xz.data(), c_xz.data() + c_xz.size());
        }

        // Collect the coefficients
        coefficients_results.push_back(coefficients_xx);
        coefficients_results.push_back(coefficients_yy);
        coefficients_results.push_back(coefficients_xy);
        if (is_3d) {
            coefficients_results.push_back(coefficients_zz);
            coefficients_results.push_back(coefficients_yz);
            coefficients_results.push_back(coefficients_xz);
        }

        return coefficients_results;
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

private:

    ///@name Private static Member Variables
    ///@{

    std::vector<std::vector<double>> mSigmaValues;
    CoordinateArrayType mCoordinates;
    std::vector<std::vector<double>> mShapeFunctionValues;
    std::vector<std::vector<double>> mShapeFunctionValuesDx;
    std::vector<std::vector<double>> mShapeFunctionValuesDy;
    std::vector<std::vector<double>> mShapeFunctionValuesDz;
    std::vector<Matrix> mCollected_D_constitutive_matrix;

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; /* Class ComputeDivSigmaUtility */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/
