//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//                   Zhiming Guo
//

// System includes

// External includes

// Project includes
#include "includes/global_variables.h"
#include "spaces/ublas_space.h"
#include "utilities/dense_householder_qr_decomposition.h"
#include "mls_shape_functions_utility.h"

namespace Kratos
{

namespace
{
    using DenseSpace = UblasSpace<double, Matrix, Vector>;
}

    double MLSShapeFunctionsUtility::CalculateKernel(
        const array_1d<double,3>& rRadVect,
        const double h)
    {
        const double c_0 = 1.0;
        const double q_squared = (rRadVect(0)*rRadVect(0) + rRadVect(1)*rRadVect(1) + rRadVect(2)*rRadVect(2)) / (h*h);
        return std::exp(-c_0*q_squared/2.0);
    }

    template<std::size_t TDim>
    void MLSShapeFunctionsUtility::CalculateKernelDerivative(
        const array_1d<double,3>& rRadVect,
        const double h,
        array_1d<double,TDim>& rKernelDerivative)
    {
        const double rad = norm_2(rRadVect);

        if (rad < 1.0e-12) {
#ifdef KRATOS_DEBUG
            KRATOS_WARNING("MLSShapeFunctionsUtility") << "Radius is close to zero. Assuming null kernel derivative." << std::endl;
#endif
            noalias(rKernelDerivative) = ZeroVector(TDim);
        } else {
            // const double q = rad / h;
            // const double der_kernel_value = (std::exp(-std::pow(q,2)) * (-2.0*q/h)) / (Globals::Pi*std::pow(h,2));
            const double der_kernel_value = (-2.0 * rad) / (std::exp(rad * rad / h / h) * Globals::Pi * h * h * h * h);
            const double rel_der_kernel_value = der_kernel_value / rad;
            for (std::size_t d = 0; d < TDim; ++d) {
                rKernelDerivative[d] = rel_der_kernel_value * rRadVect[d];
            }
        }
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double,3>& rBasis)
    {
        rBasis[0] = 1.0;
        rBasis[1] = rX[0];
        rBasis[2] = rX[1];
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double,6>& rBasis)
    {
        rBasis[0] = 1.0;
        rBasis[1] = rX[0];
        rBasis[2] = rX[1];
        rBasis[3] = rX[0]*rX[1];
        rBasis[4] = std::pow(rX[0],2);
        rBasis[5] = std::pow(rX[1],2);
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double,4>& rBasis)
    {
        rBasis[0] = 1.0;
        rBasis[1] = rX[0];
        rBasis[2] = rX[1];
        rBasis[3] = rX[2];
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double,10>& rBasis)
    {
        rBasis[0] = 1.0;
        rBasis[1] = rX[0];
        rBasis[2] = rX[1];
        rBasis[3] = rX[2];
        rBasis[4] = rX[0]*rX[1];
        rBasis[5] = rX[0]*rX[2];
        rBasis[6] = rX[1]*rX[2];
        rBasis[7] = std::pow(rX[0],2);
        rBasis[8] = std::pow(rX[1],2);
        rBasis[9] = std::pow(rX[2],2);
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double,2,3>& rBasisDerivatives)
    {
        rBasisDerivatives(0,0) = 0.0; rBasisDerivatives(0,1) = 1.0; rBasisDerivatives(0,2) = 0.0;
        rBasisDerivatives(1,0) = 0.0; rBasisDerivatives(1,1) = 0.0; rBasisDerivatives(1,2) = 1.0;
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double,3,4>& rBasisDerivatives)
    {
        rBasisDerivatives(0,0) = 0.0; rBasisDerivatives(0,1) = 1.0; rBasisDerivatives(0,2) = 0.0; rBasisDerivatives(0,3) = 0.0;
        rBasisDerivatives(1,0) = 0.0; rBasisDerivatives(1,1) = 0.0; rBasisDerivatives(1,2) = 1.0; rBasisDerivatives(1,3) = 0.0;
        rBasisDerivatives(2,0) = 0.0; rBasisDerivatives(2,1) = 0.0; rBasisDerivatives(2,2) = 0.0; rBasisDerivatives(2,3) = 1.0;
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double,2,6>& rBasisDerivatives)
    {
        rBasisDerivatives(0,0) = 0.0; rBasisDerivatives(0,1) = 1.0; rBasisDerivatives(0,2) = 0.0; rBasisDerivatives(0,3) = rX[1]; rBasisDerivatives(0,4) = 2.0*rX[0]; rBasisDerivatives(0,5) = 0.0;
        rBasisDerivatives(1,0) = 0.0; rBasisDerivatives(1,1) = 0.0; rBasisDerivatives(1,2) = 1.0; rBasisDerivatives(1,3) = rX[0]; rBasisDerivatives(1,4) = 0.0; rBasisDerivatives(1,5) = 2.0*rX[1];
    }

    void MLSShapeFunctionsUtility::EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double,3,10>& rBasisDerivatives)
    {
        rBasisDerivatives(0,0) = 0.0; rBasisDerivatives(0,1) = 1.0; rBasisDerivatives(0,2) = 0.0; rBasisDerivatives(0,3) = 0.0; rBasisDerivatives(0,4) = rX[1]; rBasisDerivatives(0,5) = rX[2]; rBasisDerivatives(0,6) = 0.0; rBasisDerivatives(0,7) = 2.0*rX[0]; rBasisDerivatives(0,8) = 0.0; rBasisDerivatives(0,9) = 0.0;
        rBasisDerivatives(1,0) = 0.0; rBasisDerivatives(1,1) = 0.0; rBasisDerivatives(1,2) = 1.0; rBasisDerivatives(1,3) = 0.0; rBasisDerivatives(1,4) = rX[0]; rBasisDerivatives(1,5) = 0.0; rBasisDerivatives(1,6) = rX[2]; rBasisDerivatives(1,7) = 0.0; rBasisDerivatives(1,8) = 2.0*rX[1]; rBasisDerivatives(1,9) = 0.0;
        rBasisDerivatives(2,0) = 0.0; rBasisDerivatives(2,1) = 0.0; rBasisDerivatives(2,2) = 0.0; rBasisDerivatives(2,3) = 1.0; rBasisDerivatives(2,4) = 0.0; rBasisDerivatives(2,5) = rX[0]; rBasisDerivatives(2,6) = rX[1]; rBasisDerivatives(2,7) = 0.0; rBasisDerivatives(2,8) = 0.0; rBasisDerivatives(2,9) = 2.0*rX[2];
    }

    template<std::size_t TDim, std::size_t TOrder>
    void MLSShapeFunctionsUtility::CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set MLS shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }

        // Set the auxiliary arrays for the L2-norm problem minimization
        static constexpr std::size_t BaseSize = TDim == 2 ? (TOrder+1)*(TOrder+2)/2 : (TOrder+1)*(TOrder+2)*(TOrder+3)/6;
        Matrix W = ZeroMatrix(n_points,n_points);
        Matrix A = ZeroMatrix(n_points,BaseSize);

        // Evaluate the L2-norm minimization problem
        array_1d<double,BaseSize> p;
        array_1d<double,3> rad_vect;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(rPoints, i_pt);
            noalias(rad_vect) = rX - r_i_pt_coords;

            // Calculate kernel values
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h);

            // Evaluate the current point polynomial basis
            EvaluatePolynomialBasis(r_i_pt_coords, p);

            // Add current point data
            W(i_pt,i_pt) = kernel;
            for (std::size_t j = 0; j < BaseSize; ++j) {
                A(i_pt, j) = kernel * p[j];
            }
        }

        // QR problem resolution
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(A);
        Matrix aux_sol(BaseSize,n_points);
        qr_decomposition.Solve(W, aux_sol);

        // Set the polynomial basis values at the point of interest
        array_1d<double,BaseSize> p0;
        EvaluatePolynomialBasis(rX, p0);

        // Get the solution for each node to obtain the corresponding MLS shape function
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            rN[i_pt] = 0.0;
            for (std::size_t j = 0; j < BaseSize; ++j) {
                rN[i_pt] += p0[j] * aux_sol(j, i_pt);
            }
        }

        KRATOS_CATCH("");
    }

    template<>
    void KRATOS_API(KRATOS_CORE) MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,1>(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Matrix& rDNDX)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set MLS shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }
        if (rDNDX.size1() != n_points || rDNDX.size2() != 2) {
            rDNDX.resize(n_points, 2, false);
        }

        // Set the auxiliary arrays for the L2-norm problem minimization
        Vector W(n_points);
        Vector DW_DX(n_points);
        Vector DW_DY(n_points);
        Matrix A = ZeroMatrix(n_points,3);
        Matrix DA_DX = ZeroMatrix(n_points,3);
        Matrix DA_DY = ZeroMatrix(n_points,3);

        // Evaluate the L2-norm minimization problem
        array_1d<double,3> p;
        array_1d<double,3> rad_vect;
        array_1d<double,2> kernel_der;
        BoundedMatrix<double,2,3> Dp_Dx;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(rPoints, i_pt);
            noalias(rad_vect) = rX - r_i_pt_coords;

            // Calculate kernel values
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(rad_vect, h, kernel_der);

            // Evaluate the current point basis
            EvaluatePolynomialBasis(r_i_pt_coords, p);
            EvaluatePolynomialBasisDerivatives(r_i_pt_coords, Dp_Dx);

            // Add current point data
            W(i_pt) = kernel;
            DW_DX(i_pt) = kernel_der[0];
            DW_DY(i_pt) = kernel_der[1];
            for (std::size_t j = 0; j < 3; ++j) {
                A(i_pt, j) = kernel * p[j];
                DA_DX(i_pt, j) = kernel_der[0] * p[j];
                DA_DY(i_pt, j) = kernel_der[1] * p[j];
            }
        }

        // QR problem resolution
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(A);

        // Set the polynomial basis values at the point of interest
        array_1d<double,3> p0;
        BoundedMatrix<double,2,3> Dp0_DX;
        EvaluatePolynomialBasis(rX, p0);
        EvaluatePolynomialBasisDerivatives(rX, Dp0_DX);

        // Do the solve for each node to obtain the corresponding MLS shape function
        Vector aux_RHS(n_points);
        Vector aux_RHS_dx_1(n_points);
        Vector aux_RHS_dy_1(n_points);
        Vector aux_RHS_dx_2(n_points);
        Vector aux_RHS_dy_2(n_points);
        Vector i_pt_sol(3);
        Vector i_pt_sol_dx_1(3);
        Vector i_pt_sol_dy_1(3);
        Vector i_pt_sol_dx_2(3);
        Vector i_pt_sol_dy_2(3);
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            aux_RHS = ZeroVector(n_points);
            aux_RHS(i_pt) = W(i_pt);
            qr_decomposition.Solve(aux_RHS, i_pt_sol);

            aux_RHS_dx_1 = prod(DA_DX, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dx_1, i_pt_sol_dx_1);

            aux_RHS_dy_1 = prod(DA_DY, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dy_1, i_pt_sol_dy_1);

            aux_RHS_dx_2 = ZeroVector(n_points);
            aux_RHS_dx_2(i_pt) = DW_DX(i_pt);
            qr_decomposition.Solve(aux_RHS_dx_2, i_pt_sol_dx_2);

            aux_RHS_dy_2 = ZeroVector(n_points);
            aux_RHS_dy_2(i_pt) = DW_DY(i_pt);
            qr_decomposition.Solve(aux_RHS_dy_2, i_pt_sol_dy_2);

            rN[i_pt] = inner_prod(p0, i_pt_sol);
            rDNDX(i_pt,0) = inner_prod(row(Dp0_DX,0), i_pt_sol) - inner_prod(p0, i_pt_sol_dx_1) + inner_prod(p0, i_pt_sol_dx_2);
            rDNDX(i_pt,1) = inner_prod(row(Dp0_DX,1), i_pt_sol) - inner_prod(p0, i_pt_sol_dy_1) + inner_prod(p0, i_pt_sol_dy_2);
        }

        KRATOS_CATCH("");
    }

    template<>
    void KRATOS_API(KRATOS_CORE) MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,2>(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Matrix& rDNDX)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set MLS shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }
        if (rDNDX.size1() != n_points || rDNDX.size2() != 2) {
            rDNDX.resize(n_points, 2, false);
        }

        // Set the auxiliary arrays for the L2-norm problem minimization
        Vector W(n_points);
        Vector DW_DX(n_points);
        Vector DW_DY(n_points);
        Matrix A = ZeroMatrix(n_points,6);
        Matrix DA_DX = ZeroMatrix(n_points,6);
        Matrix DA_DY = ZeroMatrix(n_points,6);

        // Evaluate the L2-norm minimization problem
        array_1d<double,6> p;
        array_1d<double,3> rad_vect;
        array_1d<double,2> kernel_der;
        BoundedMatrix<double,2,6> Dp_Dx;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(rPoints, i_pt);
            noalias(rad_vect) = rX - r_i_pt_coords;

            // Calculate kernel values
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(rad_vect, h, kernel_der);

            // Evaluate the current point basis
            EvaluatePolynomialBasis(r_i_pt_coords, p);
            EvaluatePolynomialBasisDerivatives(r_i_pt_coords, Dp_Dx);

            // Add current point data
            W(i_pt) = kernel;
            DW_DX(i_pt) = kernel_der[0];
            DW_DY(i_pt) = kernel_der[1];
            for (std::size_t j = 0; j < 6; ++j) {
                A(i_pt, j) = kernel * p[j];
                DA_DX(i_pt, j) = kernel_der[0] * p[j];
                DA_DY(i_pt, j) = kernel_der[1] * p[j];
            }
        }

        // QR problem resolution
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(A);

        // Set the polynomial basis values at the point of interest
        array_1d<double,6> p0;
        BoundedMatrix<double,2,6> Dp0_DX;
        EvaluatePolynomialBasis(rX, p0);
        EvaluatePolynomialBasisDerivatives(rX, Dp0_DX);

        // Do the solve for each node to obtain the corresponding MLS shape function
        Vector aux_RHS(n_points);
        Vector aux_RHS_dx_1(n_points);
        Vector aux_RHS_dy_1(n_points);
        Vector aux_RHS_dx_2(n_points);
        Vector aux_RHS_dy_2(n_points);
        Vector i_pt_sol(6);
        Vector i_pt_sol_dx_1(6);
        Vector i_pt_sol_dy_1(6);
        Vector i_pt_sol_dx_2(6);
        Vector i_pt_sol_dy_2(6);
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            aux_RHS = ZeroVector(n_points);
            aux_RHS(i_pt) = W(i_pt);
            qr_decomposition.Solve(aux_RHS, i_pt_sol);

            aux_RHS_dx_1 = prod(DA_DX, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dx_1, i_pt_sol_dx_1);

            aux_RHS_dy_1 = prod(DA_DY, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dy_1, i_pt_sol_dy_1);

            aux_RHS_dx_2 = ZeroVector(n_points);
            aux_RHS_dx_2(i_pt) = DW_DX(i_pt);
            qr_decomposition.Solve(aux_RHS_dx_2, i_pt_sol_dx_2);

            aux_RHS_dy_2 = ZeroVector(n_points);
            aux_RHS_dy_2(i_pt) = DW_DY(i_pt);
            qr_decomposition.Solve(aux_RHS_dy_2, i_pt_sol_dy_2);

            rN[i_pt] = inner_prod(p0, i_pt_sol);
            rDNDX(i_pt,0) = inner_prod(row(Dp0_DX,0), i_pt_sol) - inner_prod(p0, i_pt_sol_dx_1) + inner_prod(p0, i_pt_sol_dx_2);
            rDNDX(i_pt,1) = inner_prod(row(Dp0_DX,1), i_pt_sol) - inner_prod(p0, i_pt_sol_dy_1) + inner_prod(p0, i_pt_sol_dy_2);
        }

        KRATOS_CATCH("");
    }

    template<>
    void KRATOS_API(KRATOS_CORE) MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3,1>(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Matrix& rDNDX)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set MLS shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }
        if (rDNDX.size1() != n_points || rDNDX.size2() != 3) {
            rDNDX.resize(n_points, 3, false);
        }

        // Set the auxiliary arrays for the L2-norm problem minimization
        Vector W(n_points);
        Vector DW_DX(n_points);
        Vector DW_DY(n_points);
        Vector DW_DZ(n_points);
        Matrix A = ZeroMatrix(n_points,4);
        Matrix DA_DX = ZeroMatrix(n_points,4);
        Matrix DA_DY = ZeroMatrix(n_points,4);
        Matrix DA_DZ = ZeroMatrix(n_points,4);

        // Evaluate the L2-norm minimization problem
        array_1d<double,4> p;
        array_1d<double,3> rad_vect;
        array_1d<double,3> kernel_der;
        BoundedMatrix<double,3,4> Dp_Dx;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(rPoints, i_pt);
            noalias(rad_vect) = rX - r_i_pt_coords;

            // Calculate kernel values
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<3>(rad_vect, h, kernel_der);

            // Evaluate the current point basis
            EvaluatePolynomialBasis(r_i_pt_coords, p);
            EvaluatePolynomialBasisDerivatives(r_i_pt_coords, Dp_Dx);

            // Add current point data
            W(i_pt) = kernel;
            DW_DX(i_pt) = kernel_der[0];
            DW_DY(i_pt) = kernel_der[1];
            DW_DZ(i_pt) = kernel_der[2];
            for (std::size_t j = 0; j < 4; ++j) {
                A(i_pt, j) = kernel * p[j];
                DA_DX(i_pt, j) = kernel_der[0] * p[j];
                DA_DY(i_pt, j) = kernel_der[1] * p[j];
                DA_DZ(i_pt, j) = kernel_der[2] * p[j];
            }
        }

        // QR problem resolution
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(A);

        // Set the polynomial basis values at the point of interest
        array_1d<double,4> p0;
        BoundedMatrix<double,3,4> Dp0_DX;
        EvaluatePolynomialBasis(rX, p0);
        EvaluatePolynomialBasisDerivatives(rX, Dp0_DX);

        // Do the solve for each node to obtain the corresponding MLS shape function
        Vector aux_RHS(n_points);
        Vector aux_RHS_dx_1(n_points);
        Vector aux_RHS_dy_1(n_points);
        Vector aux_RHS_dz_1(n_points);
        Vector aux_RHS_dx_2(n_points);
        Vector aux_RHS_dy_2(n_points);
        Vector aux_RHS_dz_2(n_points);
        Vector i_pt_sol(4);
        Vector i_pt_sol_dx_1(4);
        Vector i_pt_sol_dy_1(4);
        Vector i_pt_sol_dz_1(4);
        Vector i_pt_sol_dx_2(4);
        Vector i_pt_sol_dy_2(4);
        Vector i_pt_sol_dz_2(4);
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            aux_RHS = ZeroVector(n_points);
            aux_RHS(i_pt) = W(i_pt);
            qr_decomposition.Solve(aux_RHS, i_pt_sol);

            aux_RHS_dx_1 = prod(DA_DX, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dx_1, i_pt_sol_dx_1);

            aux_RHS_dy_1 = prod(DA_DY, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dy_1, i_pt_sol_dy_1);

            aux_RHS_dz_1 = prod(DA_DZ, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dz_1, i_pt_sol_dz_1);

            aux_RHS_dx_2 = ZeroVector(n_points);
            aux_RHS_dx_2(i_pt) = DW_DX(i_pt);
            qr_decomposition.Solve(aux_RHS_dx_2, i_pt_sol_dx_2);

            aux_RHS_dy_2 = ZeroVector(n_points);
            aux_RHS_dy_2(i_pt) = DW_DY(i_pt);
            qr_decomposition.Solve(aux_RHS_dy_2, i_pt_sol_dy_2);

            aux_RHS_dz_2 = ZeroVector(n_points);
            aux_RHS_dz_2(i_pt) = DW_DZ(i_pt);
            qr_decomposition.Solve(aux_RHS_dz_2, i_pt_sol_dz_2);

            rN[i_pt] = inner_prod(p0, i_pt_sol);
            rDNDX(i_pt,0) = inner_prod(row(Dp0_DX,0), i_pt_sol) - inner_prod(p0, i_pt_sol_dx_1) + inner_prod(p0, i_pt_sol_dx_2);
            rDNDX(i_pt,1) = inner_prod(row(Dp0_DX,1), i_pt_sol) - inner_prod(p0, i_pt_sol_dy_1) + inner_prod(p0, i_pt_sol_dy_2);
            rDNDX(i_pt,2) = inner_prod(row(Dp0_DX,2), i_pt_sol) - inner_prod(p0, i_pt_sol_dz_1) + inner_prod(p0, i_pt_sol_dz_2);
        }

        KRATOS_CATCH("");
    }

    template<>
    void KRATOS_API(KRATOS_CORE) MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3,2>(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Matrix& rDNDX)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set MLS shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }
        if (rDNDX.size1() != n_points || rDNDX.size2() != 3) {
            rDNDX.resize(n_points, 3, false);
        }

        // Set the auxiliary arrays for the L2-norm problem minimization
        Vector W(n_points);
        Vector DW_DX(n_points);
        Vector DW_DY(n_points);
        Vector DW_DZ(n_points);
        Matrix A = ZeroMatrix(n_points,10);
        Matrix DA_DX = ZeroMatrix(n_points,10);
        Matrix DA_DY = ZeroMatrix(n_points,10);
        Matrix DA_DZ = ZeroMatrix(n_points,10);

        // Evaluate the L2-norm minimization problem
        array_1d<double,10> p;
        array_1d<double,3> rad_vect;
        array_1d<double,3> kernel_der;
        BoundedMatrix<double,3,10> Dp_Dx;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(rPoints, i_pt);
            noalias(rad_vect) = rX - r_i_pt_coords;

            // Calculate kernel values
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<3>(rad_vect, h, kernel_der);

            // Evaluate the current point basis
            EvaluatePolynomialBasis(r_i_pt_coords, p);
            EvaluatePolynomialBasisDerivatives(r_i_pt_coords, Dp_Dx);

            // Add current point data
            W(i_pt) = kernel;
            DW_DX(i_pt) = kernel_der[0];
            DW_DY(i_pt) = kernel_der[1];
            DW_DZ(i_pt) = kernel_der[2];
            for (std::size_t j = 0; j < 10; ++j) {
                A(i_pt, j) = kernel * p[j];
                DA_DX(i_pt, j) = kernel_der[0] * p[j];
                DA_DY(i_pt, j) = kernel_der[1] * p[j];
                DA_DZ(i_pt, j) = kernel_der[2] * p[j];
            }
        }

        // QR problem resolution
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(A);

        // Set the polynomial basis values at the point of interest
        array_1d<double,10> p0;
        BoundedMatrix<double,3,10> Dp0_DX;
        EvaluatePolynomialBasis(rX, p0);
        EvaluatePolynomialBasisDerivatives(rX, Dp0_DX);

        // Do the solve for each node to obtain the corresponding MLS shape function
        Vector aux_RHS(n_points);
        Vector aux_RHS_dx_1(n_points);
        Vector aux_RHS_dy_1(n_points);
        Vector aux_RHS_dz_1(n_points);
        Vector aux_RHS_dx_2(n_points);
        Vector aux_RHS_dy_2(n_points);
        Vector aux_RHS_dz_2(n_points);
        Vector i_pt_sol(10);
        Vector i_pt_sol_dx_1(10);
        Vector i_pt_sol_dy_1(10);
        Vector i_pt_sol_dz_1(10);
        Vector i_pt_sol_dx_2(10);
        Vector i_pt_sol_dy_2(10);
        Vector i_pt_sol_dz_2(10);
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            aux_RHS = ZeroVector(n_points);
            aux_RHS(i_pt) = W(i_pt);
            qr_decomposition.Solve(aux_RHS, i_pt_sol);

            aux_RHS_dx_1 = prod(DA_DX, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dx_1, i_pt_sol_dx_1);

            aux_RHS_dy_1 = prod(DA_DY, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dy_1, i_pt_sol_dy_1);

            aux_RHS_dz_1 = prod(DA_DZ, i_pt_sol);
            qr_decomposition.Solve(aux_RHS_dz_1, i_pt_sol_dz_1);

            aux_RHS_dx_2 = ZeroVector(n_points);
            aux_RHS_dx_2(i_pt) = DW_DX(i_pt);
            qr_decomposition.Solve(aux_RHS_dx_2, i_pt_sol_dx_2);

            aux_RHS_dy_2 = ZeroVector(n_points);
            aux_RHS_dy_2(i_pt) = DW_DY(i_pt);
            qr_decomposition.Solve(aux_RHS_dy_2, i_pt_sol_dy_2);

            aux_RHS_dz_2 = ZeroVector(n_points);
            aux_RHS_dz_2(i_pt) = DW_DZ(i_pt);
            qr_decomposition.Solve(aux_RHS_dz_2, i_pt_sol_dz_2);

            rN[i_pt] = inner_prod(p0, i_pt_sol);
            rDNDX(i_pt,0) = inner_prod(row(Dp0_DX,0), i_pt_sol) - inner_prod(p0, i_pt_sol_dx_1) + inner_prod(p0, i_pt_sol_dx_2);
            rDNDX(i_pt,1) = inner_prod(row(Dp0_DX,1), i_pt_sol) - inner_prod(p0, i_pt_sol_dy_1) + inner_prod(p0, i_pt_sol_dy_2);
            rDNDX(i_pt,2) = inner_prod(row(Dp0_DX,2), i_pt_sol) - inner_prod(p0, i_pt_sol_dz_1) + inner_prod(p0, i_pt_sol_dz_2);
        }

        KRATOS_CATCH("");
    }

    template KRATOS_API(KRATOS_CORE) void MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(const array_1d<double,3>& rRadVect, const double h, array_1d<double,2>& rKernelDerivative);
    template KRATOS_API(KRATOS_CORE) void MLSShapeFunctionsUtility::CalculateKernelDerivative<3>(const array_1d<double,3>& rRadVect, const double h, array_1d<double,3>& rKernelDerivative);
    template KRATOS_API(KRATOS_CORE) void MLSShapeFunctionsUtility::CalculateShapeFunctions<2,1>(const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN);
    template KRATOS_API(KRATOS_CORE) void MLSShapeFunctionsUtility::CalculateShapeFunctions<2,2>(const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN);
    template KRATOS_API(KRATOS_CORE) void MLSShapeFunctionsUtility::CalculateShapeFunctions<3,1>(const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN);
    template KRATOS_API(KRATOS_CORE) void MLSShapeFunctionsUtility::CalculateShapeFunctions<3,2>(const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN);

}  // namespace Kratos.
