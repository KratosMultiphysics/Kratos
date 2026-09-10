//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

// System includes

// External includes
#include <boost/numeric/ublas/matrix_proxy.hpp>

// Project includes
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{

GeometricalSensitivityUtility::GeometricalSensitivityUtility(const JacobianType& rJ, const ShapeFunctionsLocalGradientType& rDN_De)
: mrJ(rJ), mrDN_De(rDN_De)
{
    KRATOS_TRY;

    Initialize();

    KRATOS_CATCH("");
}

void GeometricalSensitivityUtility::Initialize()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrJ.size1() != mrJ.size2()) << "Non-square Jacobian matrix." << std::endl;

    KRATOS_ERROR_IF(mrJ.size2() != mrDN_De.size2())
        << "Jacobian local-coordinates size (" << mrJ.size2()
        << ") != shape function local-coordinates size(" << mrDN_De.size2()
        << ")." << std::endl;

    mCofactorJ = MathUtils<double>::CofactorMatrix(mrJ);
    mDetJ = MathUtils<double>::Det(mrJ);

    KRATOS_CATCH("");
}

void GeometricalSensitivityUtility::CalculateSensitivity(ShapeParameter Deriv, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const
{
    KRATOS_TRY;

    rDetJ_Deriv = CalculateDeterminantOfJacobianSensitivity(Deriv);

    MatrixType cofactorJ_deriv = CalculateCofactorOfJacobianSensitivity(Deriv);
    if (rDN_DX_Deriv.size1() != mrDN_De.size1() || rDN_DX_Deriv.size2() != mCofactorJ.size1())
        rDN_DX_Deriv.resize(mrDN_De.size1(), mCofactorJ.size1());
    noalias(rDN_DX_Deriv) = (1.0 / mDetJ) * prod(mrDN_De, trans(cofactorJ_deriv));
    rDN_DX_Deriv += -(rDetJ_Deriv / (mDetJ * mDetJ)) * prod(mrDN_De, trans(mCofactorJ));

    KRATOS_CATCH("");
}

double GeometricalSensitivityUtility::CalculateDeterminantOfJacobianSensitivity(ShapeParameter Deriv) const
{
    return inner_prod(
        row(mCofactorJ, Deriv.Direction),
        row(mrDN_De, Deriv.NodeIndex));
}

GeometricalSensitivityUtility::MatrixType GeometricalSensitivityUtility::CalculateCofactorOfJacobianSensitivity(
    ShapeParameter Deriv) const
{
    KRATOS_TRY;
    MatrixType result(mrJ.size1(), mrJ.size2());

    for (unsigned i = 0; i < mrJ.size1(); ++i)
    {
        if (i == Deriv.Direction)
        {
            // Here the derivative is automatically zero.
            for (unsigned j = 0; j < mrJ.size2(); ++j)
                result(i, j) = 0.0;
        }
        else
        {
            // Decrement the coordinate index if it's greater than the deleted row.
            IndexType i_coord_sub = (Deriv.Direction > i) ? Deriv.Direction - 1 : Deriv.Direction;
            for (unsigned j = 0; j < mrJ.size2(); ++j)
            {
                // Gather the first-minor submatrices explicitly (backend-generic:
                // works on the uBLAS and the Eigen-backed Matrix alike, without
                // uBLAS indirect proxies).
                MatrixType sub_jacobian(mrJ.size1() - 1, mrJ.size2() - 1);
                unsigned i_sub = 0;
                for (unsigned k = 0; k < mrJ.size1(); ++k) {
                    if (k == i) continue;
                    unsigned j_sub = 0;
                    for (unsigned l = 0; l < mrJ.size2(); ++l) {
                        if (l == j) continue;
                        sub_jacobian(i_sub, j_sub) = mrJ(k, l);
                        ++j_sub;
                    }
                    ++i_sub;
                }

                const MatrixType cofactor_sub_jacobian = MathUtils<double>::CofactorMatrix(sub_jacobian);

                // The corresponding shape function local gradients submatrix
                // (all rows, column j removed).
                MatrixType sub_DN_De(mrDN_De.size1(), mrDN_De.size2() - 1);
                for (unsigned k = 0; k < mrDN_De.size1(); ++k) {
                    unsigned j_sub = 0;
                    for (unsigned l = 0; l < mrDN_De.size2(); ++l) {
                        if (l == j) continue;
                        sub_DN_De(k, j_sub) = mrDN_De(k, l);
                        ++j_sub;
                    }
                }

                const double first_minor_deriv = inner_prod(
                    row(cofactor_sub_jacobian, i_coord_sub),
                    row(sub_DN_De, Deriv.NodeIndex));

                result(i, j) = ((i + j) % 2) ? -first_minor_deriv : first_minor_deriv;
            }
        }
    }

    return result;
    KRATOS_CATCH("");
}

} /* namespace Kratos.*/
