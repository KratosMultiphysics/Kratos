//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sebastian Ares de Parga, Juan I. Camarotti
//

#if !defined(KRATOS_rbf_shape_functions_utility_H_INCLUDED)
#define  KRATOS_rbf_shape_functions_utility_H_INCLUDED


// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "utilities/dense_householder_qr_decomposition.h"

namespace Kratos
{

/**
 * @brief Radial Basis Functions utility to calculate shape function values
 * This class uses Gaussian raidal basis functions to calculate
 * the shape function values for a given value (i.e. norm of a point) with partition of unity.
 */
class KRATOS_API(KRATOS_CORE) RBFShapeFunctionsUtility
{

public:

    ///@name Type Definitions
    ///@{

    /// Dense space definition for the QR decomposition using in the solve
    using DenseSpace = UblasSpace<double, Matrix, Vector>;

    /// QR decomposition pointer definition
    using DenseQRPointerType = typename DenseQRDecomposition<DenseSpace>::Pointer;

    /// Kratos core QR decomposition type
    using KratosCoreQRType = DenseHouseholderQRDecomposition<DenseSpace>;

    enum class RBFType {
        InverseMultiquadric,
        Multiquadric,
        Gaussian,
        ThinPlateSpline,
        WendlandC2
    };

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Inverse Multiquadric RBF
    struct InverseMultiquadric{
        double h;                      
        double operator()(double r) const {
            const double q = h * r;
            return 1.0 / std::sqrt(1.0 + q * q);
        }
    };

    /// Multiquadric RBF
    struct Multiquadric{
        double h;                      
        double operator()(double r) const {
            const double q = r / h;
            return std::sqrt(1.0 + q * q);
        }
    };

    /// Gaussian RBF
    struct Gaussian{
        double h;                      
        double operator()(double r) const {
            const double q = r / h;
            return std::exp(-0.5 * q * q);
        }
    };

    /// Thin Plate Spline RBF
    struct ThinPlateSpline {
        double operator()(double r) const {
            if (r < 1.0e-12)
                return 0.0; 
            const double r2 = r * r;
            return r2 * std::log(r2);
        }
    };

    /// Wendland C2 RBF
    struct WendlandC2 {
        double h;                      
        double operator()(double r) const {
            const double q = r / h;
            if (q >= 1.0)
                return 0.0;
            return std::pow(1.0 - q, 4) * (4.0 * q + 1.0); // (1-q)^4 * (4q+1)
        }
    };

    /**
     * @brief Calculates the RBF shape function values
     * This method calculates the RBF shape function values in one point using as
     * support the given cloud of points.
     * @param rPoints Matrix containing the coordinates of the support cloud of points
     * @param rX Coordinates where the shape functions are to be computed
     * @param h RBF shape parameter
     * @param rN Shape functions container
     */
    static void CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        DenseQRPointerType pDenseQR = nullptr,
        const RBFType RBFType = RBFType::InverseMultiquadric);

    /**
     * @brief Calculates the RBF shape function values
     * This method calculates the RBF shape function values in one point using as
     * support the given cloud of points.
     * @param rPoints Matrix containing the coordinates of the support cloud of points
     * @param rX Coordinates where the shape functions are to be computed
     * @param rN Shape functions container
     */
    static void CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        Vector& rN,
        DenseQRPointerType pDenseQR = nullptr,
        const RBFType RBFType = RBFType::InverseMultiquadric);

    /**
     * @brief Calculates the RBF shape function values
     * This method calculates the RBF shape function values in one point using as
     * support the given cloud of points.
     * @param rPoints Matrix containing the coordinates of the support cloud of points
     * @param rX Coordinates where the shape functions are to be computed
     * @param h RBF shape parameter
     * @param rN Shape functions container
     * @param Y Function to interpolate (if values of RHS are known apriori)
     */
    static double CalculateShapeFunctionsAndInterpolation(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Vector& rY,
        const RBFType RBFType = RBFType::InverseMultiquadric);

    static double CalculateInverseMultiquadricShapeParameter(const Matrix& rPoints);

    static double CalculateWendlandC2SupportRadius(const Matrix& rPoints, const double k);

    ///@}
private:
    ///@name Unaccessible methods
    ///@{

    RBFShapeFunctionsUtility(){};

    static double CalculateShapeParameter(const Matrix& rPoints, const RBFType rbf_type)
    {
        double h = 1.0;

        switch (rbf_type) {
            case RBFType::InverseMultiquadric:
            {
                h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(rPoints);
                break;
            }
            case RBFType::Multiquadric:
            {
                // Same shape parameter as the InverseMultiquadric could be used
                h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(rPoints);
                break;
            }
            case RBFType::Gaussian:
            {
                // Same shape parameter as the InverseMultiquadric could be used
                h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(rPoints);
                break;
            }
            case RBFType::ThinPlateSpline:
            {
                // TPS does not require any shape parameter
                break;
            }
            case RBFType::WendlandC2:
            {
                const double k = 1.0;
                h = RBFShapeFunctionsUtility::CalculateWendlandC2SupportRadius(rPoints, k);
                break;
            }
            default:
            {
                h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(rPoints);
                break;
            }
        }

        return h;
    }
    ///@}
};

}  // namespace Kratos.

#endif // KRATOS_rbf_shape_functions_utility_H_INCLUDED  defined