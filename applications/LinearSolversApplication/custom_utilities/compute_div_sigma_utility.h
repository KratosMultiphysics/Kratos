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
    }


    std::vector<std::vector<double>> ComputeDivergence()
    {
        std::vector<std::vector<double>> collected_coefficients = ComputeCoefficients();

        std::vector<std::vector<double>> divergence_results;
        const std::vector<double>& coefficients_xx = collected_coefficients[0];
        const std::vector<double>& coefficients_yy = collected_coefficients[1];
        const std::vector<double>& coefficients_xy = collected_coefficients[2];

        for (size_t i = 0; i < mShapeFunctionValuesDx.size(); ++i)
        {
            double d_sigma_xx_dx = 0.0;
            double d_sigma_yy_dy = 0.0;
            double d_sigma_xy_dx = 0.0;
            double d_sigma_xy_dy = 0.0;
            for (size_t j = 0; j < coefficients_xx.size(); ++j)
            {
                d_sigma_xx_dx += coefficients_xx[j] * mShapeFunctionValuesDx[i][j];
                d_sigma_yy_dy += coefficients_yy[j] * mShapeFunctionValuesDy[i][j];
                d_sigma_xy_dx += coefficients_xy[j] * mShapeFunctionValuesDx[i][j];
                d_sigma_xy_dy += coefficients_xy[j] * mShapeFunctionValuesDy[i][j];
            }
            double div_sigma_x = d_sigma_xx_dx + d_sigma_xy_dy;
            double div_sigma_y = d_sigma_xy_dx + d_sigma_yy_dy;
            
            divergence_results.push_back({div_sigma_x, div_sigma_y});
        }

        return divergence_results;
    }


    std::vector<std::vector<double>> ComputeCoefficients()
    {
        KRATOS_ERROR_IF(mSigmaValues.size() != mShapeFunctionValues.size()) 
            << "Mismatch between number of sigma values and N!" << std::endl;
        
        std::vector<std::vector<double>> coefficients_results;
        size_t gauss_point_number = mSigmaValues.size();
        // Construct design matrix A and vector b for quadratic fit
        Eigen::MatrixXd A(gauss_point_number, gauss_point_number); // 9 Gauss points, 9 shape function values per Gauss point
        Eigen::VectorXd b_xx(gauss_point_number), b_yy(gauss_point_number), b_xy(gauss_point_number);

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
        }

        // Solve for coefficients using least squares
        Eigen::VectorXd c_xx = A.colPivHouseholderQr().solve(b_xx);
        Eigen::VectorXd c_yy = A.colPivHouseholderQr().solve(b_yy);
        Eigen::VectorXd c_xy = A.colPivHouseholderQr().solve(b_xy);

        Eigen::VectorXd residual_xx = b_xx - A * c_xx;
        Eigen::VectorXd residual_yy = b_yy - A * c_yy;
        Eigen::VectorXd residual_xy = b_xy - A * c_xy;

        // Compute norms of residuals
        double norm_residual_xx = residual_xx.norm();
        double norm_residual_yy = residual_yy.norm();
        double norm_residual_xy = residual_xy.norm();

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

        // Convert Eigen vectors to std::vector and store results
        std::vector<double> coefficients_xx(c_xx.data(), c_xx.data() + c_xx.size());
        std::vector<double> coefficients_yy(c_yy.data(), c_yy.data() + c_yy.size());
        std::vector<double> coefficients_xy(c_xy.data(), c_xy.data() + c_xy.size());

        // Collect the coefficients
        coefficients_results.push_back(coefficients_xx);
        coefficients_results.push_back(coefficients_yy);
        coefficients_results.push_back(coefficients_xy);

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

