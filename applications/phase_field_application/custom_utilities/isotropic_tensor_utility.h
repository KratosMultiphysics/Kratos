//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 11 Sep 2015 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ISOTROPIC_TENSOR_UTILITY_H_INCLUDED )
#define  KRATOS_ISOTROPIC_TENSOR_UTILITY_H_INCLUDED



// System includes
#include <iostream>
#include <string>
#include <vector>
#include <cmath>


// External includes


// Project includes
#include "includes/ublas_interface.h"
#include "includes/define.h"


#define ENABLE_CHECK


// LAPACK subroutine to compute eigenvalues and eigenvectors
extern "C" void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA, double* WR, double* WI,
        int* VL, int* LDVL, int* VR, int* LDVR, double* WORK, int* LWORK, int* INFO);

extern "C" void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA,
        double* W, double* WORK, int* LWORK, int* INFO);

namespace Kratos
{

/**
 * Class representing an isotropic tensor function Y = (y1(x1,x2,x3), y2(x1,x2,x3), y3(x1,x2,x3))
 * REF:
 *  +   E. A. de Souza Neto, Computational Plasticity
 *  +   E. A. de Souza Neto, On general isotropic tensor functions of one tensor
 */
template<std::size_t TDimension>
class GeneralIsotropicTensorFunction
{
public:
    GeneralIsotropicTensorFunction() {}
    virtual ~GeneralIsotropicTensorFunction() {}

    /// Computes Y = (y1(x1,x2,x3), y2(x1,x2,x3), y3(x1,x2,x3))
    virtual void Values(const Vector& X, Vector& Y)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class of", __FUNCTION__)
    }

    /// Computes dY = [d_y1/d_x1    d_y1/d_x2    d_y1/d_x3
    ///                d_y2/d_x1    d_y2/d_x2    d_y2/d_x3
    ///                d_y3/d_x1    d_y3/d_x2    d_y3/d_x3]
    virtual void Derivatives(const Vector& X, Matrix& dY)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class of", __FUNCTION__)
    }

    ///////////////////////////////////////////////////////////////////////////

    /// Computes dY = [d_ya/d_xa    d_ya/d_xb    d_ya/d_xc
    ///                d_yb/d_xa    d_yb/d_xb    d_yb/d_xc
    ///                d_yc/d_xa    d_yc/d_xb    d_yc/d_xc]
    /// With (a,b,c) is a permutation of (0,1,2)
    void Derivatives(std::vector<int> permutation, const Vector& X, Matrix& dY)
    {
        Matrix Aux;
        this->Derivatives(X, Aux);

        if(dY.size1() != Aux.size1() || dY.size2() != Aux.size2())
            dY.resize(Aux.size1(), Aux.size2());

        for(unsigned int i = 0; i < permutation.size(); ++i)
            for(unsigned int j = 0; j < permutation.size(); ++j)
                dY(i, j) = Aux(permutation[i], permutation[j]);
    }
};

template<std::size_t TDimension>
class SingleIsotropicTensorFunction : public GeneralIsotropicTensorFunction<TDimension>
{
public:
    SingleIsotropicTensorFunction() {}
    virtual ~SingleIsotropicTensorFunction() {}

    /// Compute y(x1)
    virtual double Value(const double& x)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class of", __FUNCTION__)
    }

    /// Compute y(x1,x2)
    virtual double Value(const double& x, const double& y)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class of", __FUNCTION__)
    }

    /// Compute y(x1,x2,x3)
    virtual double Value(const double& x, const double& y, const double& z)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class of", __FUNCTION__)
    }

    /// Compute d(y(x1),(x1))
    virtual void Derivative(const double& x, Vector& Y)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class of", __FUNCTION__)
    }

    /// Compute d(y(x1,x2),(x1,x2))
    virtual void Derivative(const double& x, const double& y, Vector& Y)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class of", __FUNCTION__)
    }

    /// Compute d(y(x1,x2,x3),(x1,x2,x3))
    virtual void Derivative(const double& x, const double& y, const double& z, Vector& Y)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class of", __FUNCTION__)
    }

    /// Computes Y = (y1(x1,x2,x3), y2(x1,x2,x3), y3(x1,x2,x3))
    /// Where y1(x1,x2,x3) = y(x1,x2,x3)
    ///       y2(x1,x2,x3) = y(x2,x3,x1)
    ///       y2(x1,x2,x3) = y(x3,x1,x2)
    virtual void Values(const Vector& X, Vector& Y)
    {
        if(TDimension == 2)
        {
            Y[0] = Value(X[0], X[1]);
            Y[1] = Value(X[1], X[0]);
        }
        else if(TDimension == 3)
        {
            Y[0] = Value(X[0], X[1], X[2]);
            Y[1] = Value(X[1], X[2], X[0]);
            Y[2] = Value(X[2], X[0], X[1]);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Not supported Dimension", TDimension)
    }

    /// Computes dY = [d_y1/d_x1    d_y1/d_x2    d_y1/d_x3
    ///                d_y2/d_x1    d_y2/d_x2    d_y2/d_x3
    ///                d_y3/d_x1    d_y3/d_x2    d_y3/d_x3]
    virtual void Derivatives(const Vector& X, Matrix& dY)
    {
        if(TDimension == 2)
        {
            Vector Aux(2);
            Derivative(X[0], X[1], Aux);
            dY(0, 0) = Aux(0);
            dY(0, 1) = Aux(1);
            Derivative(X[1], X[0], Aux);
            dY(1, 0) = Aux(1);
            dY(1, 1) = Aux(0);
        }
        else if(TDimension == 3)
        {
            Vector Aux(3);
            Derivative(X[0], X[1], X[2], Aux);
            dY(0, 0) = Aux(0);
            dY(0, 1) = Aux(1);
            dY(0, 2) = Aux(2);
            Derivative(X[1], X[2], X[0], Aux);
            dY(1, 0) = Aux(2);
            dY(1, 1) = Aux(0);
            dY(1, 2) = Aux(1);
            Derivative(X[2], X[0], X[1], Aux);
            dY(2, 0) = Aux(1);
            dY(2, 1) = Aux(2);
            dY(2, 2) = Aux(0);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Not supported Dimension", TDimension)
    }
};

/// Structure for return type of spectral decomposition
struct SDParameters
{
    int p;
    std::vector<int> permutation;
    int flag;
};


/*
 * Utility to compute the spectral decomposition and derivative of isotropic tensor function
 * REF: Souza de Neto, Computational Plasticity, Appendix A
 */
template<std::size_t TDimension>
class IsotropicTensorUtility
{
public:

    #ifdef BOOST_NO_CXX11_CONSTEXPR
    static const double TOL10 = 1.0e-10;
    static const double TOL16 = 1.0e-16;
    static const double _PI_ = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
    #else
    constexpr static double TOL10 = 1.0e-10;
    constexpr static double TOL16 = 1.0e-16;
    constexpr static double _PI_ = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
    #endif

    KRATOS_CLASS_POINTER_DEFINITION(IsotropicTensorUtility);

    typedef std::size_t IndexType;
    typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;

    IsotropicTensorUtility()
    {
    }

    virtual ~IsotropicTensorUtility()
    {
    }


    /// Compute the spectral decomposition of symmetric tensor X, the eigenvalues are returned in e and the eigenprojection tensors are returned in E. The subroutine returns the number of distinct eigenvalues of tensor X.
    /// The output has this properties:
    /// +   e only contains distinct eigenvalues. E contains the respective eigenprojection tensors of the eigenvalues.
    /// +   X = sum{e[i] * E[i]}
    /// +   I = sum{E[i]}
    /// +   E[i] : E[j] = \delta_[ij]
    /// return_params.flag is nonzero if something wrong happen
    static SDParameters SpectralDecomposition(const Matrix& X, Vector& e, std::vector<Matrix>& E)
    {
        Matrix I = IdentityMatrix(TDimension);

        #ifdef ENABLE_CHECK
        bool is_symmetry = CheckSymmetry(X);
        if(!is_symmetry)
            KRATOS_THROW_ERROR(std::logic_error, "Tensor X is not symmetric:", X)
        #endif

        SDParameters return_params;
        return_params.flag = 0;

        if(TDimension == 2)
        {
            double I1 = X(0, 0) + X(1, 1);
            double I2 = X(0, 0) * X(1, 1) - X(0, 1) * X(1, 0);

            if(e.size() != 2)
                e.resize(2);

            e[0] = 0.5 * (I1 + sqrt(pow(I1, 2) - 4 * I2));
            e[1] = 0.5 * (I1 - sqrt(pow(I1, 2) - 4 * I2));

            if(!IsSame(e[0], e[1]))
            {
                if(E.size() != 2)
                    E.resize(2);

                if(E[0].size1() != 2 || E[0].size2() != 2)
                    E[0].resize(2, 2);
                noalias(E[0]) = 1.0 / (e[0] - e[1]) * (X - e[1] * I);

                if(E[1].size1() != 2 || E[1].size2() != 2)
                    E[1].resize(2, 2);
                noalias(E[1]) = 1.0 / (e[1] - e[0]) * (X - e[0] * I);

                return_params.p = 2;
            }
            else
            {
                e.resize(1);

                if(E.size() != 1)
                    E.resize(1);

                if(E[0].size1() != 2 || E[0].size2() != 2)
                    E[0].resize(2, 2);
                noalias(E[0]) = I;

                return_params.p = 1;
            }
        }
        else if(TDimension == 3)
        {
            double I1 = Trace(X);
            double I2 = 0.5 * (pow(I1, 2) - Trace(prod(X, X)));
            double I3 = Det3(X);
            double R = (-2.0 * pow(I1, 3) + 9 * I1 * I2 - 27 * I3) / 54;
            double Q = (pow(I1, 2) - 3 * I2) / 9;

            std::vector<int> permutation;
            int p;
            if(fabs(Q) < TOL16)
            {
                e.resize(1);
                e[0] = X(0, 0);
                p = 1;
            }
            else
            {
                double aux = R / sqrt(pow(Q, 3));
                double theta;
                if(IsSame(aux, 1.0))
                    theta = 0.0;
                else if(IsSame(aux, -1.0))
                    theta = _PI_;
                else
                    theta = acos(aux);

                e.resize(3);
                e[0] = -2.0 * sqrt(Q) * cos(theta / 3) + I1 / 3;
                e[1] = -2.0 * sqrt(Q) * cos((theta + 2.0 * _PI_) / 3) + I1 / 3;
                e[2] = -2.0 * sqrt(Q) * cos((theta - 2.0 * _PI_) / 3) + I1 / 3;
                p = Check3(e, permutation);
            }

            return_params.p = p;
            if(p == 1)
            {
                double tmp = e[0];
                e.resize(1);
                e[0] = tmp;

                if(E.size() != 1)
                    E.resize(1);

                E[0] = I;
            }
            else if(p == 2)
            {
                double ea = e[permutation[0]];
                double eb = e[permutation[1]];
                double ec = e[permutation[2]];

                // variant 1
//                double aux = 2.0 * pow(ea, 3) - I1 * pow(ea, 2) + I3;
//                Matrix Ea = (ea * (prod(X, X) - (I1 - ea) * X) + I3 * I) / aux;
//                // TODO: handle the case that denominator is vanished
// //                if(fabs(aux) < TOL16)
// //                    return_params.flag = 32;
//                Matrix Eb = I - Ea;

                // variant 2
                Matrix Ea = (X - eb * I) / (ea - eb);
                Matrix Eb = (X - ea * I) / (eb - ea);

                e.resize(2);
                e[0] = ea;
                e[1] = eb;

                if(E.size() != 2)
                    E.resize(2);
                E[0] = Ea;
                E[1] = Eb;

                return_params.permutation = permutation;
            }
            else if(p == 3)
            {
                if(E.size() != 3)
                    E.resize(3);

                // variant 1
//                for(unsigned int i = 0; i < 3; ++i)
//                {
//                    if(E[i].size1() != 3 || E[i].size2() != 3)
//                        E[i].resize(3, 3);
//                    double aux;
//                    if(fabs(I3) > TOL16)
//                    {
//                        aux = 2.0 * pow(e[i], 3) - I1 * pow(e[i], 2) + I3;
//                        noalias(E[i]) = (e[i] * (prod(X, X) - (I1 - e[i]) * X) + I3 * I) / aux;
// //                        if(fabs(aux) < TOL16)
// //                        {
// //                            KRATOS_WATCH(e[i])
// //                            KRATOS_WATCH(I1)
// //                            KRATOS_WATCH(I3)
// //                            KRATOS_WATCH(aux)
// //                            KRATOS_WATCH(E[i])
// //                            return_params.flag = 331;
// //                        }
//                    }
//                    else
//                    {
//                        aux = 2.0 * pow(e[i], 2) - I1;
//                        noalias(E[i]) = (prod(X, X) - (I1 - e[i]) * X) / aux;
// //                        if(fabs(aux) < TOL16)
// //                            return_params.flag = 332;
//                    }
//                }

                // variant 2
                for(unsigned int i = 0; i < 3; ++i)
                {
                    if(E[i].size1() != 3 || E[i].size2() != 3)
                        E[i].resize(3, 3);
                }
                noalias(E[0]) = prod(X - e[1] * I, X - e[2] * I) / ( (e[0] - e[1]) * (e[0] - e[2]) );
                noalias(E[1]) = prod(X - e[0] * I, X - e[2] * I) / ( (e[1] - e[0]) * (e[1] - e[2]) );
                noalias(E[2]) = prod(X - e[0] * I, X - e[1] * I) / ( (e[2] - e[0]) * (e[2] - e[1]) );
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This dimension is not supported yet:", TDimension)

        return return_params;
    }

    /// SpectralDecomposition for nonsymmetric matrix based on LAPACK
    /// this is not required by now
//    static SDParameters SpectralDecompositionNonSymLAPACK(const Matrix& X, Vector& e, std::vector<Matrix>& E)
//    {
//        char JOBVL = 'N';
//        char JOBVR = 'V';
//        int N = X.size1();
//        double A[N * N];
//        int LDA = N;
//        double WR[N]; // real part of the eigenvalues
//        double WI[N]; // imaginary part of eigenvalues
//        int LDVL = 1;
//        int LDVR = N;
//        double VL[LDVL * N]; // left eigenvectors
//        double VR[LDVR * N]; // right eigenvectors
//        int LWORK = 6 * N;
//        double WORK[LWORK];
//        int INFO;
// 
//        for(int i = 0; i < N; ++i)
//            for(int j = 0; j < N; ++j)
//                A(N * j + i) = X(i, j);
// 
//        dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, &LDVL, &LDVR, VL, VR, &LWORK, WORK, &INFO);
// 
//        if(e.size() != N)
//            e.resize(N);
//        for(int i = 0; i < N; ++i)
//            e(i) = WR[i];
// 
//        if(E.size() != N)
//            E.resize(N);
//        
//    }

    /// SpectralDecomposition for symmetric matrix based on LAPACK
    /// Use DSYEV:
    ///     http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen.html#ga442c43fca5493590f8f26cf42fed4044
    /// The output is descending order
    /// Currently this function does not produce the eigenprojection tensor correctly for the repeated eigenvalues case
    static SDParameters SpectralDecompositionSymLAPACK(const Matrix& X, Vector& e, std::vector<Matrix>& E)
    {
        char JOBZ = 'V';
        char UPLO = 'L';
        int N = X.size1();
        int LDA = N;
        double A[LDA * N];
        double W[N]; // eigenvalues
        int LWORK = 5 * N; // assume block size is 3
        double WORK[LWORK];
        int INFO;

        for(int i = 0; i < N; ++i)
            for(int j = 0; j < N; ++j)
                A[N * i + j] = X(i, j);

        dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);

        // TODO compute the eigenprojection tensor

        if(e.size() != N)
            e.resize(N);
        for(int i = 0; i < N; ++i)
            e(i) = W[N-1-i];

        if(E.size() != N)
            E.resize(N);
        for(int i = 0; i < N; ++i)
        {
            if(E[i].size1() != N || E[i].size2() != N)
                E[i].resize(N, N);

            Vector V(N);
            for(int j = 0; j < N; ++j)
                V(j) = A[N * (N-1-i) + j];
            noalias(E[i]) = outer_prod(V, V);
        }

        SDParameters return_params;
        return_params.p = 3;
        return_params.flag = INFO;

        return return_params;
    }

    /// Compute the value of isotropic tensor function Phi w.r.t tensor X. The output is a fourth order tensor.
    static void Value(const Matrix& X, const Vector& e, const std::vector<Matrix>& E, const SDParameters& decomp_params, GeneralIsotropicTensorFunction<TDimension>& Phi, Matrix& V)
    {
        if(V.size1() != TDimension || V.size2() != TDimension)
            V.resize(TDimension, TDimension);
        noalias(V) = ZeroMatrix(TDimension, TDimension);

        if(TDimension == 2)
        {
            if(decomp_params.p == 1)
            {
                Vector e_full(2);
                e_full[0] = e[0];
                e_full[1] = e[0];

                Vector P(2);
                Phi.Values(e_full, P);

                noalias(V) += P(0) * E[0];
            }
            else if(decomp_params.p == 2)
            {
                Vector P(2);
                Phi.Values(e, P);

                noalias(V) += P(0) * E[0] + P(1) * E[1];
            }
        }
        else if(TDimension == 3)
        {
            if(decomp_params.p == 1)
            {
                Vector e_full(2);
                e_full[0] = e[0];
                e_full[1] = e[0];
                e_full[2] = e[0];

                Vector P(3);
                Phi.Values(e_full, P);

                noalias(V) += P(0) * E[0];
            }
            else if(decomp_params.p == 2)
            {
                Vector P(3);

                int a = decomp_params.permutation[0];
                int b = decomp_params.permutation[1];
                int c = decomp_params.permutation[2];

                Vector e_full(3);
                e_full[a] = e[0];
                e_full[b] = e[1];
                e_full[c] = e[1];

                Phi.Values(e_full, P);

                noalias(V) += P(a) * E[0] + P(b) * E[1]; // TODO: need check
            }
            else if(decomp_params.p == 3)
            {
                Vector P(3);
                Phi.Values(e, P);

                noalias(V) += P(0) * E[0] + P(1) * E[1] + P(2) * E[2];
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This dimension is not supported yet:", TDimension)
    }

    /// Compute the derivative of isotropic tensor function Phi w.r.t tensor X. The output is a fourth order tensor.
    /// Input:
    /// +   X: input vector
    /// +   e, E, decomp_params: output of SpectralDecomposition on X
    /// +   Phi: isotropic tensor function need to compute the derivative
    /// Output:
    /// +   D: fourth order tensor
    static void Derivative(const Matrix& X, const Vector& e, const std::vector<Matrix>& E, const SDParameters& decomp_params, GeneralIsotropicTensorFunction<TDimension>& Phi, Fourth_Order_Tensor& D)
    {
        Matrix I = IdentityMatrix(TDimension);

        if(TDimension == 2)
        {
            if(decomp_params.p == 1)
            {
                Vector e_full(2);
                e_full[0] = e[0];
                e_full[1] = e[0];

                Matrix dP(2, 2);
                Phi.Derivatives(e_full, dP);

                AddFourthOrderSymmetricIdentityTensor(D, dP(0, 0) - dP(0, 1));
                ProductFourthOrderTensor(D, dP(0, 1), I, I);
            }
            else if(decomp_params.p == 2)
            {
                Vector P(2);
                Matrix dP(2, 2);
                Phi.Values(e, P);
                Phi.Derivatives(e, dP);

                AddFourthOrderSymmetricIdentityTensor(D, (P(0) - P(1)) / (e[0] - e[1]));
                ProductFourthOrderTensor(D, -(P(0) - P(1)) / (e[0] - e[1]), E[0], E[0]);
                ProductFourthOrderTensor(D, -(P(0) - P(1)) / (e[0] - e[1]), E[1], E[1]);

                for(unsigned int i = 0; i < 2; ++i)
                    for(unsigned int j = 0; j < 2; ++j)
                        ProductFourthOrderTensor(D, dP(i, j), E[i], E[j]);
            }
        }
        else if(TDimension == 3)
        {
            if(decomp_params.p == 1)
            {
                Vector e_full(3);
                e_full[0] = e[0];
                e_full[1] = e[0];
                e_full[2] = e[0];

                Matrix dP(3, 3);
                Phi.Derivatives(e_full, dP);

                AddFourthOrderSymmetricIdentityTensor(D, dP(0, 0) - dP(0, 1));
                ProductFourthOrderTensor(D, dP(0, 1), I, I);
            }
            else if(decomp_params.p == 2)
            {
                Vector P(3);
                Matrix dP(3, 3);

                int a = decomp_params.permutation[0];
                int b = decomp_params.permutation[1];
                int c = decomp_params.permutation[2];

                Vector e_full(3);
                e_full[a] = e[0];
                e_full[b] = e[1];
                e_full[c] = e[1];

                Phi.Values(e_full, P);
                Phi.Derivatives(e_full, dP);

                double s1 = (P[a] - P[c]) / pow(e_full[a] - e_full[c], 2) + (dP(c, b) - dP(c, c)) / (e_full[a] - e_full[c]);
                double s2 = 2.0 * e_full[c] * (P[a] - P[c]) / pow(e_full[a] - e_full[c], 2) + (dP(c, b) - dP(c, c)) * (e_full[a] + e_full[c]) / (e_full[a] - e_full[c]);
                double s3 = 2.0 * (P[a] - P[c]) / pow(e_full[a] - e_full[c], 3) + (dP(a, c) + dP(c, a) - dP(a, a) - dP(c, c)) / pow(e_full[a] - e_full[c], 2);
                double s4 = 2.0 * e_full[c] * (P[a] - P[c]) / pow(e_full[a] - e_full[c], 3) + (dP(a, c) - dP(c, b)) / (e_full[a] - e_full[c]) + e_full[c] * (dP(a, c) + dP(c, a) - dP(a, a) - dP(c, c)) / pow(e_full[a] - e_full[c], 2);
                double s5 = 2.0 * e_full[c] * (P[a] - P[c]) / pow(e_full[a] - e_full[c], 3) + (dP(c, a) - dP(c, b)) / (e_full[a] - e_full[c]) + e_full[c] * (dP(a, c) + dP(c, a) - dP(a, a) - dP(c, c)) / pow(e_full[a] - e_full[c], 2);
                double s6 = 2.0 * pow(e_full[c], 2) * (P[a] - P[c]) / pow(e_full[a] - e_full[c], 3) + e_full[a] * e_full[c] * (dP(a, c) + dP(c, a)) / pow(e_full[a] - e_full[c], 2) - pow(e_full[c], 2) * (dP(a, a) + dP(c, c)) / pow(e_full[a] - e_full[c], 2) - (e_full[a] + e_full[c]) * dP(c, b) / (e_full[a] - e_full[c]);

                AddFourthOrderTensorDX2DX(D, s1, X);
                AddFourthOrderSymmetricIdentityTensor(D, -s2);
                ProductFourthOrderTensor(D, -s3, X, X);
                ProductFourthOrderTensor(D, s4, X, I);
                ProductFourthOrderTensor(D, s5, I, X);
                ProductFourthOrderTensor(D, -s6, I, I);
            }
            else if(decomp_params.p == 3)
            {
                Vector P(3);
                Matrix dP(3, 3);
                Phi.Values(e, P);
                Phi.Derivatives(e, dP);
//                KRATOS_WATCH(P)
//                KRATOS_WATCH(dP)

                const int cyclic[][3] = { {0, 1, 2}, {1, 2, 0}, { 2, 0, 1} };
                for(unsigned int i = 0; i < 3; ++i)
                {
                    int a = cyclic[i][0], b = cyclic[i][1], c = cyclic[i][2];
                    double aux = P[a] / ((e[a] - e[b]) * (e[a] - e[c]));
                    AddFourthOrderTensorDX2DX(D, aux, X);
                    AddFourthOrderSymmetricIdentityTensor(D, -1.0 * aux * (e[b] + e[c]));
                    ProductFourthOrderTensor(D, -1.0 * aux * (2.0 * e[a] - e[b] - e[c]), E[a], E[a]);
                    ProductFourthOrderTensor(D, -1.0 * aux * (e[b] - e[c]), E[b], E[b]);
                    ProductFourthOrderTensor(D,        aux * (e[b] - e[c]), E[c], E[c]);
                }

                for(unsigned int i = 0; i < 3; ++i)
                    for(unsigned int j = 0; j < 3; ++j)
                        ProductFourthOrderTensor(D, dP(i, j), E[i], E[j]);
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This dimension is not supported yet:", TDimension)
    }

    ///////////////////////////////////////////////////////////////////////////
    //      Second Order Tensor subroutines                                  //
    ///////////////////////////////////////////////////////////////////////////
    static double TwoDotProduct(const Matrix& A, const Matrix& B)
    {
        double res = 0.0;
        for(unsigned int i = 0; i < TDimension; ++i)
            for(unsigned int j = 0; j < TDimension; ++j)
                res += A(i, j) * B(i, j);
        return res;
    }

    ///////////////////////////////////////////////////////////////////////////
    //      Fourth_Order_Tensor subroutines                                  //
    ///////////////////////////////////////////////////////////////////////////

    /// C += alpha A \odot B
    static void ProductFourthOrderTensor( Fourth_Order_Tensor& C, double alpha, const Matrix& A, const Matrix& B )
    {
        for(unsigned int i = 0; i < TDimension; ++i)
            for(unsigned int j = 0; j < TDimension; ++j)
                for(unsigned int k = 0; k < TDimension; ++k)
                    for(unsigned int l = 0; l < TDimension; ++l)
                        C[i][j](k, l) += alpha * A(i, j) * B(k, l);
    }

    /// C *= alpha
    static void ScaleFourthOrderTensor( Fourth_Order_Tensor& C, double alpha )
    {
        for(unsigned int i = 0; i < TDimension; ++i)
            for(unsigned int j = 0; j < TDimension; ++j)
                for(unsigned int k = 0; k < TDimension; ++k)
                    for(unsigned int l = 0; l < TDimension; ++l)
                        C[i][j](k,l) *= alpha;
    }

    /// Initialize a fourth order tensor to zero values
    static void InitializeFourthOrderTensor(Fourth_Order_Tensor& C)
    {
        if(C.size() != TDimension)
            C.resize(TDimension);
        for(unsigned int i = 0; i < TDimension; ++i)
        {
            if(C[i].size() != TDimension)
                C[i].resize(TDimension);
            for(unsigned int j = 0; j < TDimension; ++j)
            {
                if(C[i][j].size1() != TDimension || C[i][j].size2() != TDimension)
                    C[i][j].resize(TDimension, TDimension);
                noalias(C[i][j]) = ZeroMatrix(TDimension, TDimension);
            }
        }
    }

    /// Create a fourth order tensor and initialize it to zero
    static Fourth_Order_Tensor CreateFourthOrderTensor()
    {
        Fourth_Order_Tensor C;
        InitializeFourthOrderTensor(C);
        return C;
    }

    /// REF: eq (2.110) Souza de Neto, Computational Plasticity
    static void AddFourthOrderSymmetricIdentityTensor(Fourth_Order_Tensor& C, double alpha)
    {
        for(unsigned int i = 0; i < TDimension; ++i)
            for(unsigned int j = 0; j < TDimension; ++j)
                for(unsigned int k = 0; k < TDimension; ++k)
                    for(unsigned int l = 0; l < TDimension; ++l)
                        C[i][j](k,l) += alpha * 0.5 * (Kronecker(i, k) * Kronecker(j, l) + Kronecker(i, l) * Kronecker(j, k));
    }

    /// REF: eq (2.100) Souza de Neto, Computational Plasticity
    static void AddFourthOrderIdentityTensor(Fourth_Order_Tensor& C, double alpha)
    {
        for(unsigned int i = 0; i < TDimension; ++i)
            for(unsigned int j = 0; j < TDimension; ++j)
                for(unsigned int k = 0; k < TDimension; ++k)
                    for(unsigned int l = 0; l < TDimension; ++l)
                        C[i][j](k,l) += alpha * Kronecker(i, k) * Kronecker(j, l);
    }

    ///
    static void AddFourthOrderElasticTensor(Fourth_Order_Tensor& C, double alpha, const double E, const double NU)
    {
        double lambda = NU * E / ((1 + NU) * (1 - 2 * NU));
        double mu     = E / (2 * (1 + NU));
        for(unsigned int i = 0; i < TDimension; ++i)
            for(unsigned int j = 0; j < TDimension; ++j)
                for(unsigned int k = 0; k < TDimension; ++k)
                    for(unsigned int l = 0; l < TDimension; ++l)
                        C[i][j](k,l) += alpha * (lambda * Kronecker(i, j) * Kronecker(k, l) + mu * (Kronecker(i, k) * Kronecker(j, l) + Kronecker(i, l) * Kronecker(j, k)));
    }

    /// REF: eq (A.46) Souza de Neto, Computational Plasticity
    static void AddFourthOrderTensorDX2DX(Fourth_Order_Tensor& C, double alpha, const Matrix& X)
    {
        for(unsigned int i = 0; i < TDimension; ++i)
            for(unsigned int j = 0; j < TDimension; ++j)
                for(unsigned int k = 0; k < TDimension; ++k)
                    for(unsigned int l = 0; l < TDimension; ++l)
                        C[i][j](k,l) += alpha * 0.5 * (Kronecker(i, k) * X(l, j) + Kronecker(i, l) * X(k, j) + Kronecker(j, l) * X(i, k) + Kronecker(k, j) * X(i, l));
    }

    /// Transfer fourth order tensor to second order tensor
    static void FourthOrderTensorToMatrix(const Fourth_Order_Tensor& C, Matrix& M)
    {
        if(TDimension == 3)
        {
            if(M.size1() != 6 || M.size2() != 6)
                M.resize(6, 6);

            M(0, 0) = C[0][0](0, 0);
            M(0, 1) = C[0][0](1, 1);
            M(0, 2) = C[0][0](2, 2);
            M(0, 3) = C[0][0](0, 1);
            M(0, 4) = C[0][0](0, 2);
            M(0, 5) = C[0][0](1, 2);

            M(1, 0) = C[1][1](0, 0);
            M(1, 1) = C[1][1](1, 1);
            M(1, 2) = C[1][1](2, 2);
            M(1, 3) = C[1][1](0, 1);
            M(1, 4) = C[1][1](0, 2);
            M(1, 5) = C[1][1](1, 2);

            M(2, 0) = C[2][2](0, 0);
            M(2, 1) = C[2][2](1, 1);
            M(2, 2) = C[2][2](2, 2);
            M(2, 3) = C[2][2](0, 1);
            M(2, 4) = C[2][2](0, 2);
            M(2, 5) = C[2][2](1, 2);

            M(3, 0) = C[0][1](0, 0);
            M(3, 1) = C[0][1](1, 1);
            M(3, 2) = C[0][1](2, 2);
            M(3, 3) = C[0][1](0, 1);
            M(3, 4) = C[0][1](0, 2);
            M(3, 5) = C[0][1](1, 2);

            M(4, 0) = C[0][2](0, 0);
            M(4, 1) = C[0][2](1, 1);
            M(4, 2) = C[0][2](2, 2);
            M(4, 3) = C[0][2](0, 1);
            M(4, 4) = C[0][2](0, 2);
            M(4, 5) = C[0][2](1, 2);

            M(5, 0) = C[1][2](0, 0);
            M(5, 1) = C[1][2](1, 1);
            M(5, 2) = C[1][2](2, 2);
            M(5, 3) = C[1][2](0, 1);
            M(5, 4) = C[1][2](0, 2);
            M(5, 5) = C[1][2](1, 2);
        }
        else if(TDimension == 2)
        {
            if(M.size1() != 3 || M.size2() != 3)
                M.resize(3, 3);

            M(0, 0) = C[0][0](0, 0);
            M(0, 1) = C[0][0](1, 1);
            M(0, 2) = C[0][0](0, 1);
            M(1, 0) = C[1][1](0, 0);
            M(1, 1) = C[1][1](1, 1);
            M(1, 2) = C[1][1](0, 1);
            M(2, 0) = C[0][1](0, 0);
            M(2, 1) = C[0][1](1, 1);
            M(2, 2) = C[0][1](0, 1);
        }
    }

    static void StrainVectorToTensor(const Vector& StrainVector, Matrix& StrainTensor)
    {
        if(TDimension == 2)
        {
            if(StrainTensor.size1() != 2 || StrainTensor.size2() != 2)
                StrainTensor.resize(2, 2);

            StrainTensor(0, 0) = StrainVector(0);
            StrainTensor(0, 1) = 0.5 * StrainVector(2);
            StrainTensor(1, 0) = 0.5 * StrainVector(2);
            StrainTensor(1, 1) = StrainVector(1);
        }
        else if(TDimension == 3)
        {
            if(StrainTensor.size1() != 3 || StrainTensor.size2() != 3)
                StrainTensor.resize(3, 3);

            StrainTensor(0, 0) = StrainVector(0);
            StrainTensor(0, 1) = 0.5 * StrainVector(3);
            StrainTensor(0, 2) = 0.5 * StrainVector(5);
            StrainTensor(1, 0) = StrainTensor(0, 1);
            StrainTensor(1, 1) = StrainVector(1);
            StrainTensor(1, 2) = 0.5 * StrainVector(4);
            StrainTensor(2, 0) = StrainTensor(0, 2);
            StrainTensor(2, 1) = StrainTensor(1, 2);
            StrainTensor(2, 2) = StrainVector(2);
        }
    }

    static void StressVectorToTensor(const Vector& StressVector, Matrix& StressTensor)
    {
        if(TDimension == 2)
        {
            if(StressTensor.size1() != 2 || StressTensor.size2() != 2)
                StressTensor.resize(2, 2);

            StressTensor(0, 0) = StressVector(0);
            StressTensor(0, 1) = StressVector(2);
            StressTensor(1, 0) = StressVector(2);
            StressTensor(1, 1) = StressVector(1);
        }
        else if(TDimension == 3)
        {
            if(StressTensor.size1() != 3 || StressTensor.size2() != 3)
                StressTensor.resize(3, 3);

            StressTensor(0, 0) = StressVector(0);
            StressTensor(0, 1) = StressVector(3);
            StressTensor(0, 2) = StressVector(5);
            StressTensor(1, 0) = StressVector(3);
            StressTensor(1, 1) = StressVector(1);
            StressTensor(1, 2) = StressVector(4);
            StressTensor(2, 0) = StressVector(5);
            StressTensor(2, 1) = StressVector(4);
            StressTensor(2, 2) = StressVector(2);
        }
    }

    static void StressTensorToVector(const Matrix& StressTensor, Vector& StressVector)
    {
        if(TDimension == 2)
        {
            if(StressVector.size() != 3)
                StressVector.resize(3);

	        StressVector[0] = StressTensor(0, 0);
	        StressVector[1] = StressTensor(1, 1);
	        StressVector[2] = StressTensor(0, 1);
        }
        else if(TDimension == 3)
        {
            if(StressVector.size() != 6)
                StressVector.resize(6);

            StressVector[0] = StressTensor(0, 0);
            StressVector[1] = StressTensor(1, 1);
            StressVector[2] = StressTensor(2, 2);
            StressVector[3] = StressTensor(0, 1);
            StressVector[4] = StressTensor(1, 2);
            StressVector[5] = StressTensor(0, 2);
        }
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "IsotropicTensorUtility" << std::endl;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsotropicTensorUtility" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    /// Kronecker operator
    static inline double Kronecker(int i, int j)
    {
        return (i == j) ? 1.0 : 0.0;
    }

    /// Compute trace of matrix
    inline static double Trace(const Matrix& X)
    {
        if(X.size1() != X.size2())
            KRATOS_THROW_ERROR(std::logic_error, "This operation only works with square matrix", __FUNCTION__)

        double tr = 0.0;
        for(unsigned int i = 0; i < X.size1(); ++i)
            tr += X(i, i);
        return tr;
    }

    /// Compute determinant of 2x2 matrix
    inline static double Det2(const Matrix& A)
    {
        return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
    }

    /// Check if 2 number are the same within the tolerance
    inline static bool IsSame(const double& a, const double& b)
    {
        if(fabs(a) + fabs(b) > TOL16)
            return fabs(a - b) / (fabs(a) + fabs(b)) < TOL16;
        else
            return true;
    }

    /// Compute determinant of 3x3 matrix
    inline static double Det3(const Matrix& A)
    {
        double a = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
        double b = A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0);
        double c = A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0);
        return A(0, 0) * a - A(0, 1) * b + A(0, 2) * c;
    }

    /// Compute the number of distinct values of vector e of dimension 3. If it returns 2, the permutation (a,b,c) will be return such as e[a] <> e[b] == e[c]
    static int Check3(const Vector& e, std::vector<int>& permutation)
    {
        if(e.size() != 3)
            KRATOS_THROW_ERROR(std::logic_error, "Dimension is invalid", "")

        if(IsSame(e[0], e[1]))
        {
            if(IsSame(e[1], e[2]))
            {
                return 1;
            }
            else
            {
                if(permutation.size() != 3)
                    permutation.resize(3);
                permutation[0] = 2;
                permutation[1] = 0;
                permutation[2] = 1;
                return 2;
            }
        }
        else
        {
            if(IsSame(e[1], e[2]))
            {
                if(permutation.size() != 3)
                    permutation.resize(3);
                permutation[0] = 0;
                permutation[1] = 1;
                permutation[2] = 2;
                return 2;
            }
            else
            {
                if(IsSame(e[0], e[2]))
                {
                    if(permutation.size() != 3)
                        permutation.resize(3);
                    permutation[0] = 1;
                    permutation[1] = 0;
                    permutation[2] = 2;
                    return 2;
                }
                else
                    return 3;
            }
        }
    }

    static bool CheckSymmetry(const Matrix& X)
    {
        for(unsigned int i = 1; i < X.size1(); ++i)
            for(unsigned int j = 0; j < i; ++j)
            {
                if(!IsSame(X(i, j), X(j, i)))
                    return false;
            }
        return true;
    }

    static double GetMachineEpsilon()
    {
        double x = 1.0;
        while((1.0 + x/2.0) != 1.0)
            x /= 2.0;
        return x;
    }

};

}  // namespace Kratos.

#undef ENABLE_CHECK

#endif // KRATOS_ISOTROPIC_TENSOR_UTILITY_H_INCLUDED  defined 

