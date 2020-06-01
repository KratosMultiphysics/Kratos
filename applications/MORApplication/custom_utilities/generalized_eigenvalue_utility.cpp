
// System includes

// External includes

// Project includes
#include "generalized_eigenvalue_utility.h"

namespace Kratos
{
namespace GeneralizedEigenvalueUtility
{
    void SolveGEP(const Eigen::MatrixXcd& rA, const Eigen::MatrixXcd& rB, Eigen::Map<Eigen::ArrayXcd>& rEV)
    {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
        ces.compute(rB.lu().solve(rA));

        rEV = ces.eigenvalues();
    }

    void SolveGEP(const Eigen::MatrixXd& rA, const Eigen::MatrixXd& rB, Eigen::Map<Eigen::ArrayXcd>& rEV)
    {
        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
        ges.compute(rA,rB,false);

        rEV = ges.eigenvalues();
    }

} // namespace GeneralizedEigenvalueUtility

} // namespace Kratos