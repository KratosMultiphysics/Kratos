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
        const CoordinateArrayType& rCoordinates)
    {
        mSigmaValues = rSigmaValues;
        mCoordinates = rCoordinates;
    }

    void SetInputData2(
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
        std::vector<std::vector<double>> collected_derivatives = ComputeDerivatives();

        // Vector to store the computed divergence at each Gauss point
        std::vector<std::vector<double>> divergence_results;

        // Iterate over each Gauss point
        for (const auto& derivatives : collected_derivatives)
        {
            // The expected structure of `derivatives`:
            // {d_sigma_xx_dx, d_sigma_yy_dy, d_sigma_xy_dx, d_sigma_xy_dy}

            double d_sigma_xx_dx = derivatives[0];
            double d_sigma_yy_dy = derivatives[1];
            double d_sigma_xy_dx = derivatives[2];
            double d_sigma_xy_dy = derivatives[3];

            // Compute the divergence components
            double div_sigma_x = d_sigma_xx_dx + d_sigma_xy_dy;
            double div_sigma_y = d_sigma_xy_dx + d_sigma_yy_dy;

            // Store the results
            divergence_results.push_back({div_sigma_x, div_sigma_y});
        }

        return divergence_results;
    }

    std::vector<std::vector<double>> ComputeDivergence2()
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
                // d_sigma_xx_dx += coefficients_xx[j] * mShapeFunctionValues[i][j];
                // d_sigma_yy_dy += coefficients_yy[j] * mShapeFunctionValues[i][j];
                // d_sigma_xy_dx += coefficients_xy[j] * mShapeFunctionValues[i][j];
                
            }
            // KRATOS_WATCH(d_sigma_xx_dx)
            // KRATOS_WATCH(d_sigma_yy_dy)
            // KRATOS_WATCH(d_sigma_xy_dx)
            // KRATOS_WATCH(mSigmaValues[i])

            double div_sigma_x = d_sigma_xx_dx + d_sigma_xy_dy;
            double div_sigma_y = d_sigma_xy_dx + d_sigma_yy_dy;
            
            divergence_results.push_back({div_sigma_x, div_sigma_y});
        }

        return divergence_results;
    }

    /**
     * @brief Perform computations on the provided sigma and coordinate values
     */
    std::vector<std::vector<double>> ComputeDerivatives()
    {
        KRATOS_ERROR_IF(mSigmaValues.size() != mCoordinates.size()) 
            << "Mismatch between number of sigma values and coordinates!" << std::endl;
        
        std::vector<std::vector<double>> divergence_results;

        for (std::size_t i = 0; i < mSigmaValues.size(); ++i)
        {
            const std::vector<double>& sigma = mSigmaValues[i];
            auto coord = mCoordinates[i];

            KRATOS_ERROR_IF(sigma.size() != 3) 
                << "Sigma vector does not have 3 components at index " << i << std::endl;
        }

        if (mCoordinates.size() == 4) // p=1 case
        {
            // Construct design matrix A and vector b for linear fit
            Eigen::MatrixXd A(4, 3); // 4 Gauss points, 3 coefficients
            Eigen::VectorXd b_xx(4), b_yy(4), b_xy(4);
            for (size_t i = 0; i < 4; ++i)
            {
                const auto& coord = mCoordinates[i];
                double x = coord[0];
                double y = coord[1];
                // Prepare the matrix A
                A(i, 0) = 1.0;
                A(i, 1) = x;
                A(i, 2) = y;
                // Prepare the three b vectors for sigma_xx, sigma__yy and sigma_xy
                b_xx(i) = mSigmaValues[i][0];
                b_yy(i) = mSigmaValues[i][1];
                b_xy(i) = mSigmaValues[i][2];
            }
            // Solve for coefficients using least squares
            Eigen::VectorXd c_xx = A.colPivHouseholderQr().solve(b_xx);
            Eigen::VectorXd c_yy = A.colPivHouseholderQr().solve(b_yy);
            Eigen::VectorXd c_xy = A.colPivHouseholderQr().solve(b_xy);
            // Derivative is simply the coefficients a1 (d/dx) and a2 (d/dy)
            for (size_t i = 0; i < 4; ++i)
            {
                const auto& coord = mCoordinates[i];
                double x = coord[0];
                double y = coord[1];
                // Derivatives
                double d_sigma_xx_dx = c_xx(1);
                double d_sigma_yy_dy = c_yy(2);
                double d_sigma_xy_dx = c_xy(1);
                double d_sigma_xy_dy = c_xy(2);
                // Store in the result
                divergence_results.push_back({d_sigma_xx_dx, d_sigma_yy_dy, d_sigma_xy_dx, d_sigma_xy_dy});
            }
        }
        else if (mCoordinates.size() == 9) // p=2 case
        {
            // Construct design matrix A and vector b for quadratic fit
            Eigen::MatrixXd A(9, 6); // 9 Gauss points, 6 coefficients for a quadratic polynomial
            Eigen::VectorXd b_xx(9), b_yy(9), b_xy(9);
            for (size_t i = 0; i < 9; ++i)
            {
                const auto& coord = mCoordinates[i];
                double x = coord[0];
                double y = coord[1];
                // Prepare the matrix A for quadratic terms
                A(i, 0) = 1.0;
                A(i, 1) = x;
                A(i, 2) = y;
                A(i, 3) = x * x;
                A(i, 4) = x * y;
                A(i, 5) = y * y;
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
            // KRATOS_WATCH(norm_residual_xx)
            // KRATOS_WATCH(norm_residual_yy)
            // KRATOS_WATCH(norm_residual_xy)

            // Define a threshold for acceptable residuals
            const double threshold = 1e-5; // Adjust this value as necessary

            // Check if residuals are too high and issue warnings or errors
            if (norm_residual_xx > threshold) {
                std::cerr << "Warning: High residual detected for sigma_xx fit! Norm: " << norm_residual_xx << std::endl;
            }
            if (norm_residual_yy > threshold) {
                std::cerr << "Warning: High residual detected for sigma_yy fit! Norm: " << norm_residual_yy << std::endl;
            }
            if (norm_residual_xy > threshold) {
                std::cerr << "Warning: High residual detected for sigma_xy fit! Norm: " << norm_residual_xy << std::endl;
            }
            // Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(A);
            // Eigen::VectorXd c_xx = cod.solve(b_xx);
            // Eigen::VectorXd c_yy = cod.solve(b_yy);
            // Eigen::VectorXd c_xy = cod.solve(b_xy);
            // Derivatives
            for (size_t i = 0; i < 9; ++i)
            {
                double d_sigma_xx_dx = c_xx(1) + 2 * c_xx(3) * mCoordinates[i][0] + c_xx(4) * mCoordinates[i][1];
                double d_sigma_yy_dy = c_yy(2) + 2 * c_yy(5) * mCoordinates[i][1] + c_yy(4) * mCoordinates[i][0];
                double d_sigma_xy_dx = c_xy(1) + 2 * c_xy(3) * mCoordinates[i][0] + c_xy(4) * mCoordinates[i][1];
                double d_sigma_xy_dy = c_xy(2) + 2 * c_xy(5) * mCoordinates[i][1] + c_xy(4) * mCoordinates[i][0];
                divergence_results.push_back({d_sigma_xx_dx, d_sigma_yy_dy, d_sigma_xy_dx, d_sigma_xy_dy});
            }

            // // Construct design matrix A and vector b for quadratic fit
            // Eigen::MatrixXd A(9, 3); // 9 Gauss points, 6 coefficients for a quadratic polynomial
            // Eigen::VectorXd b_xx(9), b_yy(9), b_xy(9);
            // for (size_t i = 0; i < 9; ++i)
            // {
            //     const auto& coord = mCoordinates[i];
            //     double x = coord[0];
            //     double y = coord[1];
            //     // Prepare the matrix A for quadratic terms
            //     A(i, 0) = 1.0;
            //     A(i, 1) = x;
            //     A(i, 2) = y;
            //     // Prepare the three b vectors for sigma_xx, sigma_yy, and sigma_xy
            //     b_xx(i) = mSigmaValues[i][0];
            //     b_yy(i) = mSigmaValues[i][1];
            //     b_xy(i) = mSigmaValues[i][2];
            // }

            // // Solve for coefficients using least squares
            // Eigen::VectorXd c_xx = A.colPivHouseholderQr().solve(b_xx);
            // Eigen::VectorXd c_yy = A.colPivHouseholderQr().solve(b_yy);
            // Eigen::VectorXd c_xy = A.colPivHouseholderQr().solve(b_xy);
            // // Derivatives
            // for (size_t i = 0; i < 9; ++i)
            // {
            //     double d_sigma_xx_dx = c_xx(1);
            //     double d_sigma_yy_dy = c_yy(2);
            //     double d_sigma_xy_dx = c_xy(1);
            //     double d_sigma_xy_dy = c_xy(2);
            //     divergence_results.push_back({d_sigma_xx_dx, d_sigma_yy_dy, d_sigma_xy_dx, d_sigma_xy_dy});
            // }
        }
        else if (mCoordinates.size() == 16) // p=3 case
        {
            // Construct design matrix A and vector b for cubic fit
            Eigen::MatrixXd A(16, 10); // 16 Gauss points, 10 coefficients for a cubic polynomial
            Eigen::VectorXd b_xx(16), b_yy(16), b_xy(16);

            for (size_t i = 0; i < 16; ++i)
            {
                const auto& coord = mCoordinates[i];
                double x = coord[0];
                double y = coord[1];
                // Prepare the matrix A for cubic terms
                A(i, 0) = 1.0;
                A(i, 1) = x;
                A(i, 2) = y;
                A(i, 3) = x * x;
                A(i, 4) = x * y;
                A(i, 5) = y * y;
                A(i, 6) = x * x * x;
                A(i, 7) = x * x * y;
                A(i, 8) = x * y * y;
                A(i, 9) = y * y * y;
                // Prepare the three b vectors for sigma_xx, sigma_yy, and sigma_xy
                b_xx(i) = mSigmaValues[i][0];
                b_yy(i) = mSigmaValues[i][1];
                b_xy(i) = mSigmaValues[i][2];
            }

            // Solve for coefficients using least squares
            Eigen::VectorXd c_xx = A.colPivHouseholderQr().solve(b_xx);
            Eigen::VectorXd c_yy = A.colPivHouseholderQr().solve(b_yy);
            Eigen::VectorXd c_xy = A.colPivHouseholderQr().solve(b_xy);

            // Derivatives
            for (size_t i = 0; i < 16; ++i)
            {
                double d_sigma_xx_dx = c_xx(1) + 2 * c_xx(3) * mCoordinates[i][0] + c_xx(4) * mCoordinates[i][1] + 3 * c_xx(6) * mCoordinates[i][0] * mCoordinates[i][0] + 2 * c_xx(7) * mCoordinates[i][0] * mCoordinates[i][1] + c_xx(8) * mCoordinates[i][1] * mCoordinates[i][1];
                double d_sigma_yy_dy = c_yy(2) + 2 * c_yy(5) * mCoordinates[i][1] + c_yy(4) * mCoordinates[i][0] + 3 * c_yy(9) * mCoordinates[i][1] * mCoordinates[i][1] + 2 * c_yy(8) * mCoordinates[i][1] * mCoordinates[i][0] + c_yy(7) * mCoordinates[i][0] * mCoordinates[i][0];
                double d_sigma_xy_dx = c_xy(1) + 2 * c_xy(3) * mCoordinates[i][0] + c_xy(4) * mCoordinates[i][1] + 3 * c_xy(6) * mCoordinates[i][0] * mCoordinates[i][0] + 2 * c_xy(7) * mCoordinates[i][0] * mCoordinates[i][1] + c_xy(8) * mCoordinates[i][1] * mCoordinates[i][1];
                double d_sigma_xy_dy = c_xy(2) + 2 * c_xy(5) * mCoordinates[i][1] + c_xy(4) * mCoordinates[i][0] + 3 * c_xy(9) * mCoordinates[i][1] * mCoordinates[i][1] + 2 * c_xy(8) * mCoordinates[i][1] * mCoordinates[i][0] + c_xy(7) * mCoordinates[i][0] * mCoordinates[i][0];

                divergence_results.push_back({d_sigma_xx_dx, d_sigma_yy_dy, d_sigma_xy_dx, d_sigma_xy_dy});
            }
        } else { KRATOS_ERROR << "p > 3, still need to be implemented";}


        return divergence_results;
    }













    std::vector<std::vector<double>> ComputeCoefficients()
    {
        KRATOS_ERROR_IF(mSigmaValues.size() != mShapeFunctionValues.size()) 
            << "Mismatch between number of sigma values and N!" << std::endl;
        
        std::vector<std::vector<double>> coefficients_results;
        int gauss_point_number = mSigmaValues.size();

        // if (mCoordinates.size() == 9) // p=2 case
        // {
            // Construct design matrix A and vector b for quadratic fit
            Eigen::MatrixXd A(gauss_point_number, gauss_point_number); // 9 Gauss points, 9 shape function values per Gauss point
            Eigen::VectorXd b_xx(gauss_point_number), b_yy(gauss_point_number), b_xy(gauss_point_number);

            for (size_t i = 0; i < gauss_point_number; ++i)
            {
                // Use the shape function values directly to build matrix A
                const std::vector<double>& shape_function_values = mShapeFunctionValues[i];

                // Assuming mShapeFunctionValues[i] has 9 values for each Gauss point
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
            // KRATOS_WATCH(norm_residual_xx)
            // KRATOS_WATCH(norm_residual_yy)
            // KRATOS_WATCH(norm_residual_xy)

            // Define a threshold for acceptable residuals
            const double threshold = 1e-5; // Adjust this value as necessary

            // Check if residuals are too high and issue warnings or errors
            if (norm_residual_xx > threshold) {
                std::cerr << "Warning: High residual detected for sigma_xx fit! Norm: " << norm_residual_xx << std::endl;
            }
            if (norm_residual_yy > threshold) {
                std::cerr << "Warning: High residual detected for sigma_yy fit! Norm: " << norm_residual_yy << std::endl;
            }
            if (norm_residual_xy > threshold) {
                std::cerr << "Warning: High residual detected for sigma_xy fit! Norm: " << norm_residual_xy << std::endl;
            }

            // Convert Eigen vectors to std::vector and store results
            std::vector<double> coefficients_xx(c_xx.data(), c_xx.data() + c_xx.size());
            std::vector<double> coefficients_yy(c_yy.data(), c_yy.data() + c_yy.size());
            std::vector<double> coefficients_xy(c_xy.data(), c_xy.data() + c_xy.size());

            // Collect the coefficients
            coefficients_results.push_back(coefficients_xx);
            coefficients_results.push_back(coefficients_yy);
            coefficients_results.push_back(coefficients_xy);

        // }
        // else { KRATOS_ERROR << "p > 3, still need to be implemented";}


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

