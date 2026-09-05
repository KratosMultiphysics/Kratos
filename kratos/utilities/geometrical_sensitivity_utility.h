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

#if !defined(KRATOS_GEOMETRICAL_SENSITIVITY_UTILITY_H_INCLUDED)
#define KRATOS_GEOMETRICAL_SENSITIVITY_UTILITY_H_INCLUDED

// System includes
#include <memory>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

struct ShapeParameter
{
    std::size_t NodeIndex;
    std::size_t Direction;
    class Sequence;
};

class ShapeParameter::Sequence
    {
    public:
        Sequence(std::size_t NumberOfNodes, std::size_t Dimension)
            : mNumberOfNodes(NumberOfNodes), mDimension(Dimension)
        {
        }

        operator bool() const
        {
            return (mShapeParameter.NodeIndex < mNumberOfNodes);
        }

        ShapeParameter& CurrentValue()
        {
            return mShapeParameter;
        }

        const ShapeParameter& CurrentValue() const
        {
            return mShapeParameter;
        }

        Sequence& operator++()
        {
            KRATOS_ERROR_IF_NOT(*this)
                << "Increment is out of sequence's range.\n";
            mShapeParameter.Direction = (mShapeParameter.Direction + 1) % mDimension;
            if (mShapeParameter.Direction == 0)
                ++mShapeParameter.NodeIndex;
            return *this;
        }

    private:
        const std::size_t mNumberOfNodes;
        const std::size_t mDimension;
        ShapeParameter mShapeParameter = {0, 0};
    };

class KRATOS_API(KRATOS_CORE) GeometricalSensitivityUtility
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GeometricalSensitivityUtility);

    // The interface types follow the (backend-selected) Kratos::Matrix so the
    // utility works directly on the matrices the geometries and elements
    // produce; the minor extraction is done by explicit gathering (see the
    // implementation) instead of uBLAS indirect proxies.
    typedef Matrix MatrixType;

    typedef MatrixType JacobianType;

    typedef MatrixType ShapeFunctionsLocalGradientType;

    typedef MatrixType ShapeFunctionsGradientType;

    typedef unsigned IndexType;

    template <class T>
    using matrix_row = boost::numeric::ublas::matrix_row<T>;

    ///@}
    ///@name Life Cycle
    ///@{

    GeometricalSensitivityUtility(const JacobianType& rJ, const ShapeFunctionsLocalGradientType& rDN_De);

    ///@}
    ///@name Operations
    ///@{

    void CalculateSensitivity(ShapeParameter Deriv, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const;

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

    void Initialize();

    double CalculateDeterminantOfJacobianSensitivity(ShapeParameter Deriv) const;

    MatrixType CalculateCofactorOfJacobianSensitivity(ShapeParameter Deriv) const;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/

#endif /* KRATOS_GEOMETRICAL_SENSITIVITY_UTILITY_H_INCLUDED defined */
