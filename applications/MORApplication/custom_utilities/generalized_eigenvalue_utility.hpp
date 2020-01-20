#if !defined(KRATOS_GENERALIZED_EIGENVALUE_UTILITY_H_INCLUDED)
#define KRATOS_GENERALIZED_EIGENVALUE_UTILITY_H_INCLUDED

// System includes

// External includes
#include <Eigen/Eigenvalues>

// Project includes
#include "includes/define.h"
#include "includes/ublas_complex_interface.h"

namespace Kratos
{
namespace GeneralizedEigenvalueUtility
{

    template <typename DenseSpaceType>
    void Compute(typename DenseSpaceType::MatrixType& rA0, 
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

        Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>> ev(rEigenvalues.data().begin(), rEigenvalues.size());
        // Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1>> x(rX.data().begin(), rX.size());

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
        Eigen::GeneralizedEigenSolver<EigenMatrixType> ges;
        ges.compute(A,B,false);

        ev = ges.eigenvalues();

    }

} // namespace GeneralizedEigenvalueUtility

} // namespace Kratos

#endif