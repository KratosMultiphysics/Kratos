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
#include <type_traits>

// External includes
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// Project includes
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class GeometricalSensitivityUtility::Impl
{
public:
    ///@name Type Definitions
    ///@{

    using IndirectArrayType = boost::numeric::ublas::indirect_array<boost::numeric::ublas::vector<std::size_t>>;

    using SubMatrixType = boost::numeric::ublas::matrix_indirect<const MatrixType, IndirectArrayType>;

    template <class T>
    using matrix_row = boost::numeric::ublas::matrix_row<T>;

    ///@}
    ///@name Life Cycle
    ///@{

    Impl(const JacobianType& rJ, const ShapeFunctionsLocalGradientType& rDN_De)
        : mrJ(rJ), mrDN_De(rDN_De)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize();

    void CalculateSensitivity(IndexType iNode, IndexType iCoord, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    const JacobianType& mrJ;
    const ShapeFunctionsLocalGradientType& mrDN_De;
    MatrixType mCofactorJ;
    double mDetJ;

    ///@}
    ///@name Private Operations
    ///@{

    double CalculateDeterminantOfJacobianSensitivity(IndexType iNode, IndexType iCoord) const;

    MatrixType CalculateCofactorOfJacobianSensitivity(IndexType iNode, IndexType iCoord) const;

    ///@}
};

///@} // Kratos Classes

GeometricalSensitivityUtility::GeometricalSensitivityUtility(const JacobianType& rJ, const ShapeFunctionsLocalGradientType& rDN_De)
: mpImpl(new Impl(rJ, rDN_De))
{
    KRATOS_TRY;
    mpImpl->Initialize();
    KRATOS_CATCH("");
}

GeometricalSensitivityUtility::~GeometricalSensitivityUtility() = default;

void GeometricalSensitivityUtility::CalculateSensitivity(IndexType iNode, IndexType iCoord, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const
{
    KRATOS_TRY;
    mpImpl->CalculateSensitivity(iNode, iCoord, rDetJ_Deriv, rDN_DX_Deriv);
    KRATOS_CATCH("");
}

void GeometricalSensitivityUtility::Impl::Initialize()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrJ.size1() != mrJ.size2()) << "Non-square Jacobian matrix." << std::endl;

    KRATOS_ERROR_IF(mrJ.size2() != mrDN_De.size2())
        << "Jacobian local-coordinates size (" << mrJ.size2()
        << ") != shape function local-coordinates size(" << mrDN_De.size2()
        << ")." << std::endl;

    mCofactorJ = MathUtils<double>::CofactorMatrix(mrJ);
    mDetJ = MathUtils<double>::DetMat(mrJ);

    KRATOS_CATCH("");
}

void GeometricalSensitivityUtility::Impl::CalculateSensitivity(IndexType iNode, IndexType iCoord, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const
{
    KRATOS_TRY;

    rDetJ_Deriv = CalculateDeterminantOfJacobianSensitivity(iNode, iCoord);

    MatrixType cofactorJ_deriv = CalculateCofactorOfJacobianSensitivity(iNode, iCoord);
    if (rDN_DX_Deriv.size1() != mrDN_De.size1() || rDN_DX_Deriv.size2() != mCofactorJ.size1())
        rDN_DX_Deriv.resize(mrDN_De.size1(), mCofactorJ.size1());
    noalias(rDN_DX_Deriv) = (1.0 / mDetJ) * prod(mrDN_De, trans(cofactorJ_deriv));
    rDN_DX_Deriv += -(rDetJ_Deriv / (mDetJ * mDetJ)) * prod(mrDN_De, trans(mCofactorJ));

    KRATOS_CATCH("");
}

double GeometricalSensitivityUtility::Impl::CalculateDeterminantOfJacobianSensitivity(
    IndexType iNode, IndexType iCoord) const
{
    return inner_prod(matrix_row<const MatrixType>(mCofactorJ, iCoord),
                      matrix_row<const ShapeFunctionsLocalGradientType>(mrDN_De, iNode));
}

GeometricalSensitivityUtility::MatrixType GeometricalSensitivityUtility::Impl::CalculateCofactorOfJacobianSensitivity(
    IndexType iNode, IndexType iCoord) const
{
    KRATOS_TRY;
    MatrixType result(mrJ.size1(), mrJ.size2());

    IndirectArrayType ia3(mrDN_De.size1());
    for (std::size_t k = 0; k < ia3.size(); ++k)
        ia3[k] = k;

    for (unsigned i = 0; i < mrJ.size1(); ++i)
    {
        if (i == iCoord)
        {
            // Here the derivative is automatically zero.
            for (unsigned j = 0; j < mrJ.size2(); ++j)
                result(i, j) = 0.0;
        }
        else
        {
            // Decrement the coordinate index if it's greater than the deleted row.
            IndexType i_coord_sub = (iCoord > i) ? iCoord - 1 : iCoord;
            for (unsigned j = 0; j < mrJ.size2(); ++j)
            {
                IndirectArrayType ia1(mrJ.size1() - 1), ia2(mrJ.size2() - 1);

                // Construct the Jacobian submatrix structure for the first minor.
                unsigned i_sub = 0;
                for (unsigned k = 0; k < mrJ.size1(); ++k)
                    if (k != i)
                        ia1(i_sub++) = k;

                unsigned j_sub = 0;
                for (unsigned k = 0; k < mrJ.size2(); ++k)
                    if (k != j)
                        ia2(j_sub++) = k;

                const SubMatrixType sub_jacobian(mrJ, ia1, ia2);
                const MatrixType cofactor_sub_jacobian = MathUtils<double>::CofactorMatrix(sub_jacobian);

                // Construct the corresponding shape function local gradients
                // submatrix.
                const SubMatrixType sub_DN_De(mrDN_De, ia3, ia2);

                const double first_minor_deriv = inner_prod(
                    matrix_row<const MatrixType>(cofactor_sub_jacobian, i_coord_sub),
                    matrix_row<const SubMatrixType>(sub_DN_De, iNode));

                result(i, j) = ((i + j) % 2) ? -first_minor_deriv : first_minor_deriv;
            }
        }
    }

    return result;
    KRATOS_CATCH("");
}

} /* namespace Kratos.*/
