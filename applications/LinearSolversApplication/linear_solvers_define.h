/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Armin Geiser
*/



#if !defined(KRATOS_LINEARSOLVERS_DEFINE_H_INCLUDED)
#define KRATOS_LINEARSOLVERS_DEFINE_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Kratos
{

// IMPORTANT:
//
// The matrices in Kratos use row-major storage ordering.
// Eigen matrices, however, use column-major storage ordering by default.
//
// If these two worlds are combined, it is crucial that Eigen matrices are used
// with **row-major** storage order to match the Kratos matrices!
//
// To simplify things, the most commonly used types are defined below:

template<typename _Scalar> using EigenDynamicMatrix = Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
template<typename _Scalar> using EigenDynamicVector = Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>;

template<typename _Scalar> using EigenSparseMatrix = Eigen::SparseMatrix<_Scalar, Eigen::RowMajor, int>;

} // namespace Kratos

#endif // defined(KRATOS_LINEARSOLVERS_DEFINE_H_INCLUDED)
