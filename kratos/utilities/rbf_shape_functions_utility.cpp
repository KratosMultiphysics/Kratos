//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sebastian Ares de Parga
//

// System includes

// External includes

// Project includes
#include "includes/global_variables.h"
#include "utilities/math_utils.h"
#include "rbf_shape_functions_utility.h"

namespace Kratos
{
    // 1) Matrix-only version
    template <class TOperation>
    void BuildRBFSystem(
        Matrix& rMatrixTarget,
        const TOperation& rOperation,
        const Matrix& rPoints) noexcept
    {
        const std::size_t n_points = rPoints.size1();
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                const double norm_xij = norm_2(row(rPoints,i_pt) - row(rPoints,j_pt));
                rMatrixTarget(i_pt,j_pt) = rOperation(norm_xij);
            }
        }
    }

    // 2) Matrix + vector version
    template <class TOperation>
    void BuildRBFSystem(
        Matrix& rMatrixTarget,
        Vector& rVectorTarget,
        const TOperation& rOperation,
        const array_1d<double,3>& rX,
        const Matrix& rPoints) noexcept
    {
        const std::size_t n_points = rPoints.size1();
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                const double norm_xij = norm_2(row(rPoints,i_pt) - row(rPoints,j_pt));
                rMatrixTarget(i_pt,j_pt) = rOperation(norm_xij);
            }

            const double normX = norm_2(rX - row(rPoints,i_pt));
            rVectorTarget[i_pt] = rOperation(normX);
        }
    }

    template <class TOperation>
    void ExecuteRBFInterpolation(double& rInterpolation,
                            const Vector& rN,
                            const TOperation& rOperation,
                            const array_1d<double,3>& rX,
                            const Matrix& rPoints) noexcept
    {
        const std::size_t n_points = rPoints.size1();
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            const double norm_xi = norm_2(rX-row(rPoints,i_pt));
            const double Phi = rOperation(norm_xi);
            rInterpolation += Phi*rN[i_pt];
        }
    }

    void RBFShapeFunctionsUtility::CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        DenseQRPointerType pDenseQR,
        const RBFType RBFType)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set RBF shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }

        // Initialize the RBF interpolation matrix and RBF interpolated vector
        Matrix A = ZeroMatrix(n_points,n_points);
        Vector Phi = ZeroVector(n_points);

        switch (RBFType) {
            case RBFType::InverseMultiquadric: BuildRBFSystem(A, Phi, RBFShapeFunctionsUtility::InverseMultiquadric {h}, rX, rPoints);    break;
            case RBFType::Multiquadric:        BuildRBFSystem(A, Phi, RBFShapeFunctionsUtility::Multiquadric {h}, rX, rPoints);           break;
            case RBFType::Gaussian:            BuildRBFSystem(A, Phi, RBFShapeFunctionsUtility::Gaussian {h}, rX, rPoints);               break;
            case RBFType::ThinPlateSpline:     BuildRBFSystem(A, Phi, RBFShapeFunctionsUtility::ThinPlateSpline(), rX, rPoints);          break;
            case RBFType::WendlandC2:          BuildRBFSystem(A, Phi, RBFShapeFunctionsUtility::WendlandC2 {h}, rX, rPoints);             break;
            default: KRATOS_ERROR << "Unhandled RBF operation type";
        }

        // Obtain the RBF shape functions (N)
        // Note that we do a QR solve if a QR decomposition pointer is provided
        // Otherwise we solve the system with standard LU factorization
        if (pDenseQR) {
            pDenseQR->Compute(A);
            pDenseQR->Solve(Phi, rN);
        } else {
            MathUtils<double>::Solve(A, rN, Phi);
        }

        // Partition of unity
        noalias(rN) = rN/sum(rN);

        KRATOS_CATCH("");
    }

    void RBFShapeFunctionsUtility::CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        Vector& rN,
        DenseQRPointerType pDenseQR,
        const RBFType RBFType)
    {
        KRATOS_TRY;

        double h = CalculateShapeParameter(rPoints, RBFType);

        CalculateShapeFunctions(rPoints, rX, h, rN, pDenseQR);

        KRATOS_CATCH("");
    }

    double RBFShapeFunctionsUtility::CalculateShapeFunctionsAndInterpolation(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Vector& rY,
        const RBFType RBFType)
    {
        double interpolation = 0;
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set RBF shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }

        // Initialize the RBF interpolation matrix and RBF interpolated vector
        Matrix A = ZeroMatrix(n_points,n_points);

        // Build the RBF interpolation matrix and RBF interpolated vector
        switch (RBFType) {
            case RBFType::InverseMultiquadric: BuildRBFSystem(A, RBFShapeFunctionsUtility::InverseMultiquadric {h}, rPoints);    break;
            case RBFType::Multiquadric:        BuildRBFSystem(A, RBFShapeFunctionsUtility::Multiquadric {h}, rPoints);           break;
            case RBFType::Gaussian:            BuildRBFSystem(A, RBFShapeFunctionsUtility::Gaussian {h}, rPoints);               break;
            case RBFType::ThinPlateSpline:     BuildRBFSystem(A, RBFShapeFunctionsUtility::ThinPlateSpline(), rPoints);          break;
            case RBFType::WendlandC2:          BuildRBFSystem(A, RBFShapeFunctionsUtility::WendlandC2 {h}, rPoints);             break;
            default: KRATOS_ERROR << "Unhandled RBF operation type";
        }

        // Obtain the RBF shape functions (N)
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(A);
        qr_decomposition.Solve(rY, rN);

        // Interpolate solution
        switch (RBFType) {
            case RBFType::InverseMultiquadric: ExecuteRBFInterpolation(interpolation, rN, RBFShapeFunctionsUtility::InverseMultiquadric {h}, rX, rPoints);    break;
            case RBFType::Multiquadric:        ExecuteRBFInterpolation(interpolation, rN, RBFShapeFunctionsUtility::Multiquadric {h}, rX, rPoints);           break;
            case RBFType::Gaussian:            ExecuteRBFInterpolation(interpolation, rN, RBFShapeFunctionsUtility::Gaussian {h}, rX, rPoints);               break;
            case RBFType::ThinPlateSpline:     ExecuteRBFInterpolation(interpolation, rN, RBFShapeFunctionsUtility::ThinPlateSpline(), rX, rPoints);     break;
            case RBFType::WendlandC2:          ExecuteRBFInterpolation(interpolation, rN, RBFShapeFunctionsUtility::WendlandC2 {h}, rX, rPoints);        break;
            default: KRATOS_ERROR << "Unhandled RBF operation type";
        }

        KRATOS_CATCH("");

        return interpolation;
    }

    double RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(const Matrix& rPoints)
    {
        // Find shape parameter http://www.math.iit.edu/~fass/Dolomites.pdf [Hardy]
        double d = 0; // Total distance of nearest x_i neighbors
        const std::size_t n_points = rPoints.size1();
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            double d_nearest_neighbor = std::numeric_limits<double>::max(); // Initialize distance to nearest x_i neighbor
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                if (i_pt!=j_pt){ // Avoid measuring distance between same points
                    const double norm_xij = norm_2(row(rPoints,i_pt)-row(rPoints,j_pt));
                    if (norm_xij < d_nearest_neighbor){
                        d_nearest_neighbor = norm_xij;
                    }
                }

            }
            d += d_nearest_neighbor;
        }
        d /= n_points;
        KRATOS_ERROR_IF(d < 1.0e-12) << "Nearest neighbours distance is close to zero. Check that the cloud points are not overlapping." << std::endl;
        return 1/(0.815*d);// Shape parameter (Only for inverted multiquadratic)
    }

    // Extracted from: Review of coupling methods for non-matching meshes (https://doi.org/10.1016/j.cma.2006.03.017)
    double RBFShapeFunctionsUtility::CalculateWendlandC2SupportRadius(const Matrix& rPoints, const double k = 2.5)
    {
        const std::size_t n_points = rPoints.size1();
        KRATOS_ERROR_IF(n_points < 2) << "At least two points are required to estimate spacing." << std::endl;

        double total_distance = 0.0;
        std::size_t count = 0;

        for (std::size_t i = 0; i < n_points; ++i) {
            for (std::size_t j = i + 1; j < n_points; ++j) {
                const double distance = norm_2(row(rPoints, i) - row(rPoints, j));
                total_distance += distance;
                ++count;
            }
        }

        const double average_spacing = total_distance / count;
        return k * average_spacing; // Support radius for Wendland C2
    }

}  // namespace Kratos.