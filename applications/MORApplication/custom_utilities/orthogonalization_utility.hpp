#if !defined(KRATOS_ORTHOGONALIZATION_UTILITY_H_INCLUDED)
#define KRATOS_ORTHOGONALIZATION_UTILITY_H_INCLUDED

// System includes

// External includes
#include <Eigen/QR>

// Project includes
#include "includes/define.h"

namespace Kratos
{
namespace OrthogonalizationUtility
{

    template <typename TSparseSpaceType>
    void OrthogonalizeQR(typename TSparseSpaceType::MatrixType& rA)
    {
        typedef typename TSparseSpaceType::DataType ScalarType;

        Eigen::Map<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> A(rA.data().begin(), rA.size1(), rA.size2());
        Eigen::HouseholderQR<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> qr(A);

        Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> thinQ;
        Eigen::MatrixXd I(Eigen::MatrixXd::Identity(rA.size1(),rA.size2()));

        thinQ = qr.householderQ() * I;
        A = thinQ;
    }

} // namespace OrthogonalizationUtility

} // namespace Kratos

#endif