//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_ELEMENT_DATA_DERIVATIVE_H)
#define KRATOS_ELEMENT_DATA_DERIVATIVE_H

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim>
class ElementDataDerivative
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    ///}
    ///@name Life Cycle
    ///@{

    ElementDataDerivative(
        const IndexType NodeIndex,
        const IndexType DirectionIndex,
        const GeometryType& rGeometry,
        const double W,
        const Vector& rN,
        const Matrix& rdNdX,
        const double WDerivative,
        const double DetJDerivative,
        const Matrix& rdNdXDerivative);

    ///@}

protected:
    ///@name Private Members
    ///@{

    const IndexType mNodeIndex;
    const IndexType mDirectionIndex;
    const GeometryType& mrGeometry;
    const double mW;
    const Vector& mrN;
    const Matrix& mrdNdX;
    const double mWDerivative;
    const double mDetJDerivative;
    const Matrix& mrdNdXDerivative;

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_ELEMENT_DATA_DERIVATIVE_H