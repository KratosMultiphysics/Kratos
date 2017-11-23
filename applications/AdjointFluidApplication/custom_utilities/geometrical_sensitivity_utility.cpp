//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

// System includes
#include <type_traits>
#include <iomanip> // DEBUG

// External includes
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// Application includes
#include "custom_utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

class GeometricalSensitivityUtility::Impl
{
public:
    ///@name Type Definitions
    ///@{

    using IndirectArrayType = boost::numeric::ublas::indirect_array<boost::numeric::ublas::vector<std::size_t>>;

    using SubMatrixType = boost::numeric::ublas::matrix_indirect<const MatrixType, IndirectArrayType>;

    using SubSubMatrixType = boost::numeric::ublas::matrix_indirect<const SubMatrixType, IndirectArrayType>;
    
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

    template<class TMatrixType>
    double CalculateDeterminant(const TMatrixType& rMat) const;

    template<class TMatrixType>
    double CalculateCofactor(const TMatrixType& rMat, IndexType i, IndexType j) const;

    template<class TMatrixType>
    MatrixType CalculateCofactorMatrix(const TMatrixType& rMat) const;

    double CalculateDeterminantOfJacobianSensitivity(IndexType iNode, IndexType iCoord) const;

    MatrixType CalculateCofactorOfJacobianSensitivity(IndexType iNode, IndexType iCoord) const;

    ///@}
};

///@} // Kratos Classes
///@} // Adjoint Fluid Application group

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

    mCofactorJ = CalculateCofactorMatrix(mrJ);
    mDetJ = CalculateDeterminant(mrJ);

    if (mrJ.size1() == 3)
    {
        MatrixType cofactorJ_sensitivity = CalculateCofactorOfJacobianSensitivity(5, 2);

        std::cout << "mrJ = " << std::setprecision(12) << mrJ << std::endl;
        std::cout << "mCofactorJ = " << mCofactorJ << std::endl;
        std::cout << "cofactorJ_sensitivity = " << cofactorJ_sensitivity << std::endl;
    }
    // std::cout << "Det(J) = " << mDetJ << std::endl;
    // std::cout << "Cofactor(J) = " << mCofactorJ << std::endl;

    KRATOS_CATCH("");
}

void GeometricalSensitivityUtility::Impl::CalculateSensitivity(IndexType iNode, IndexType iCoord, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const
{
    KRATOS_TRY;

    rDetJ_Deriv = CalculateDeterminantOfJacobianSensitivity(iNode, iCoord);

    KRATOS_CATCH("");
}

template<class TMatrixType>
double GeometricalSensitivityUtility::Impl::CalculateDeterminant(const TMatrixType& rMat) const
{
    KRATOS_TRY;

    //KRATOS_ERROR_IF(rMat.size1() != rMat.size2())
    //    << "Non-square matrix detected." << std::endl;

    double det;
    if (rMat.size1() == 1)
        det = rMat(0, 0);
    else if (rMat.size1() == 2)
        det = rMat(0, 0) * rMat(1, 1) - rMat(1, 0) * rMat(0, 1);
    else if (rMat.size1() == 3)
        det = rMat(0, 0) * rMat(1, 1) * rMat(2, 2) +
              rMat(1, 0) * rMat(2, 1) * rMat(0, 2) +
              rMat(0, 1) * rMat(1, 2) * rMat(2, 0) -
              rMat(2, 0) * rMat(1, 1) * rMat(0, 2) -
              rMat(2, 1) * rMat(1, 2) * rMat(0, 0) -
              rMat(1, 0) * rMat(0, 1) * rMat(2, 2);
    else
        KRATOS_ERROR << "Matrix dimension = " << rMat.size1()
                     << " is not supported." << std::endl;

    return det;
    KRATOS_CATCH("");
}

template <class TMatrixType>
double GeometricalSensitivityUtility::Impl::CalculateCofactor(const TMatrixType& rMat, IndexType i, IndexType j) const
{
    static_assert(std::is_same<TMatrixType, MatrixType>::value ||
                      std::is_same<TMatrixType, SubMatrixType>::value,
                  "Invalid template parameter.");
    KRATOS_TRY;

    if (rMat.size1() == 1 && rMat.size2() == 1)
        return 1.0;

    IndirectArrayType ia1(rMat.size1() - 1), ia2(rMat.size2() - 1);

    // Construct the submatrix for the first minor.
    unsigned i_sub = 0;
    for (unsigned k = 0; k < rMat.size1(); ++k)
        if (k != i)
            ia1(i_sub++) = k;

    unsigned j_sub = 0;
    for (unsigned k = 0; k < rMat.size2(); ++k)
        if (k != j)
            ia2(j_sub++) = k;

    double first_minor;
    typename std::conditional<std::is_same<TMatrixType, MatrixType>::value, SubMatrixType, SubSubMatrixType>::type sub_mat(
        rMat, ia1, ia2);
    first_minor = CalculateDeterminant(sub_mat);

    return ((i + j) % 2) ? -first_minor : first_minor;
    KRATOS_CATCH("");
}

template <class TMatrixType>
GeometricalSensitivityUtility::MatrixType GeometricalSensitivityUtility::Impl::CalculateCofactorMatrix(const TMatrixType& rMat) const
{
    static_assert(std::is_same<TMatrixType, MatrixType>::value ||
                      std::is_same<TMatrixType, SubMatrixType>::value,
                  "Invalid template parameter.");
    KRATOS_TRY;

    MatrixType cofactor_matrix(rMat.size1(), rMat.size2());

    for (unsigned i = 0; i < rMat.size1(); ++i)
        for (unsigned j = 0; j < rMat.size2(); ++j)
            cofactor_matrix(i, j) = CalculateCofactor(rMat, i, j);

    return cofactor_matrix;

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
            IndexType i_coord_sub = (i > iCoord) ? iCoord : iCoord - 1;
            for (unsigned j = 0; j < mrJ.size2(); ++j)
            {
                IndirectArrayType ia1(mrJ.size1() - 1), ia2(mrJ.size2() - 1);

                // Construct the Jacobian sub-matrix for the first minor.
                unsigned i_sub = 0;
                for (unsigned k = 0; k < mrJ.size1(); ++k)
                    if (k != i)
                        ia1(i_sub++) = k;

                unsigned j_sub = 0;
                for (unsigned k = 0; k < mrJ.size2(); ++k)
                    if (k != j)
                        ia2(j_sub++) = k;

                const SubMatrixType sub_jacobian(mrJ, ia1, ia2);
                const MatrixType cofactor_sub_jacobian =
                    CalculateCofactorMatrix(sub_jacobian);
                    
                // Construct the corresponding shape function local gradients
                // sub-matrix.
                const SubMatrixType sub_DN_De(mrDN_De, ia3, ia2);

                // if (i == 1 && j == 2)
                // {
                //     std::cout << "sub_jacobian = " << sub_jacobian << std::endl;
                //     std::cout << "cofactor_sub_jacobian = " << cofactor_sub_jacobian << std::endl;
                //     std::cout << "mrDN_De = " << mrDN_De << std::endl;
                //     std::cout << "sub_DN_De = " << sub_DN_De << std::endl;
                // }

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
