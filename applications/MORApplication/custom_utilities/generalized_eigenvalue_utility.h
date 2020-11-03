#if !defined(KRATOS_GENERALIZED_EIGENVALUE_UTILITY_H_INCLUDED)
#define KRATOS_GENERALIZED_EIGENVALUE_UTILITY_H_INCLUDED

// System includes

// External includes
#include <Eigen/Eigenvalues>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"

namespace Kratos
{
namespace GeneralizedEigenvalueUtility
{

    void KRATOS_API(MOR_APPLICATION) SolveGEP(const Eigen::MatrixXcd& rA, const Eigen::MatrixXcd& rB, Eigen::Map<Eigen::ArrayXcd>& rEV);

    void KRATOS_API(MOR_APPLICATION) SolveGEP(const Eigen::MatrixXd& rA, const Eigen::MatrixXd& rB, Eigen::Map<Eigen::ArrayXcd>& rEV);

    template <typename DenseSpaceType>
    void ComputePolynomial(typename DenseSpaceType::MatrixType& rA0,
        typename DenseSpaceType::MatrixType& rA1,
        typename DenseSpaceType::MatrixType& rA2,
        ComplexVector& rEigenvalues)
    {
        typedef typename DenseSpaceType::DataType ScalarType;
        typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixType;
        const size_t system_size = rA0.size1();

        if( rEigenvalues.size() != 2*system_size )
            rEigenvalues.resize(2*system_size);

        Eigen::Map<EigenMatrixType> A0(rA0.data().begin(), rA0.size1(), rA0.size2());
        Eigen::Map<EigenMatrixType> A1(rA1.data().begin(), rA1.size1(), rA1.size2());
        Eigen::Map<EigenMatrixType> A2(rA2.data().begin(), rA2.size1(), rA2.size2());

        // Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>> ev(rEigenvalues.data().begin(), rEigenvalues.size());
        Eigen::Map<Eigen::ArrayXcd> ev(rEigenvalues.data().begin(), rEigenvalues.size());

        // linearization (L1)
        EigenMatrixType I(EigenMatrixType::Identity(system_size, system_size));
        EigenMatrixType A(EigenMatrixType::Zero(2*system_size, 2*system_size));
        EigenMatrixType B(EigenMatrixType::Zero(2*system_size, 2*system_size));

        A.block(0, system_size, system_size, system_size) = I;
        A.block(system_size, 0, system_size, system_size) = -A0;
        A.block(system_size, system_size, system_size, system_size) = -A1;

        B.block(0, 0, system_size, system_size) = I;
        B.block(system_size, system_size, system_size, system_size) = A2;

        // linearization (L2)
        // A.block(0,0,system_size,system_size) = -A0;
        // A.block(system_size,system_size,system_size,system_size) = I;

        // B.block(0,0,system_size,system_size) = A1;
        // B.block(0, system_size, system_size, system_size) = A2;
        // B.block(system_size, 0, system_size, system_size) = I;

        // solve generalized eigenvalue problem
        // Eigen::GeneralizedEigenSolver<EigenMatrixType> ges;
        // ges.compute(A,B,false);

        // ev = ges.eigenvalues();
        // dispatch
        SolveGEP(A, B, ev);
    }

    /**
     * @brief Computes the generalized eigenvalue problem A*v = lambda*B*v
     * @param rA A square matrix
     * @param rB A square matrix
     * @param rEigenvalues vector of eigenvalues
     * @details This utility pairs computes the complex generalized eigenvalue problem for two matrices A and B.
     *      The gep is converted to a standard eigenvalue problem, because Eigen does not support complex generalized
     *      eigenvalue problems. Therefore B has to be invertible.
     * @author Quirin Aumann
     */
    template <typename DenseSpaceType>
    void Compute(typename DenseSpaceType::MatrixType& rA, typename DenseSpaceType::MatrixType& rB, ComplexVector& rEigenvalues)
    {
        typedef typename DenseSpaceType::DataType ScalarType;
        typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixType;
        const size_t system_size = rA.size1();

        if( rEigenvalues.size() != system_size )
            rEigenvalues.resize(system_size);

        Eigen::Map<EigenMatrixType> A(rA.data().begin(), rA.size1(), rA.size2());
        Eigen::Map<EigenMatrixType> B(rB.data().begin(), rB.size1(), rB.size2());
        Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>> ev(rEigenvalues.data().begin(), rEigenvalues.size());

        Eigen::ComplexEigenSolver<EigenMatrixType> ces;
        ces.compute(B.lu().solve(A));

        KRATOS_ERROR_IF(ces.eigenvalues().size() != ev.size() ) << "No eigenvalues could be obtained. Check input matrices." << std::endl;

        ev = ces.eigenvalues();
    }

}; // namespace GeneralizedEigenvalueUtility

} // namespace Kratos

#endif