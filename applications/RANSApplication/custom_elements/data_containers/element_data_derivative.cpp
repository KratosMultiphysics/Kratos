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

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

// Application includes

// Include base h
#include "element_data_derivative.h"

namespace Kratos
{
template <unsigned int TDim>
ElementDataDerivative<TDim>::ElementDataDerivative(
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const GeometryType& rGeometry,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : mNodeIndex(NodeIndex),
      mDirectionIndex(DirectionIndex),
      mrGeometry(rGeometry),
      mW(W),
      mrN(rN),
      mrdNdX(rdNdX),
      mWDerivative(WDerivative),
      mDetJDerivative(DetJDerivative),
      mrdNdXDerivative(rdNdXDerivative)
{
    KRATOS_TRY

    const IndexType number_of_nodes = rGeometry.PointsNumber();

    KRATOS_DEBUG_ERROR_IF(mrN.size() != number_of_nodes)
        << "mrN vector is not initialized properly. [ mrN.size() != "
           "number_of_nodes, mrN.size() = "
        << mrN.size() << ", number_of_nodes = " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdX.size1() != number_of_nodes)
        << "mrdNdX matrix is not initialized properly. [ mrdNdX.size1() != "
           "number_of_nodes, mrdNdX.size1() = "
        << mrdNdX.size1() << ", number_of_nodes = " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdX.size2() != TDim)
        << "mrdNdX matrix is not initialized properly. [ mrdNdX.size2() != "
           "TDim, mrdNdX.size2() = "
        << mrdNdX.size2() << ", TDim = " << TDim << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdXDerivative.size1() != number_of_nodes)
        << "mrdNdXDerivative matrix is not initialized properly. [ "
           "mrdNdXDerivative.size1() != "
           "number_of_nodes, mrdNdXDerivative.size1() = "
        << mrdNdXDerivative.size1() << ", number_of_nodes = " << number_of_nodes
        << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mrdNdXDerivative.size2() != TDim)
        << "mrdNdXDerivative matrix is not initialized properly. [ "
           "mrdNdXDerivative.size2() != "
           "TDim, mrdNdXDerivative.size2() = "
        << mrdNdXDerivative.size2() << ", TDim = " << TDim << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mNodeIndex >= number_of_nodes)
        << "Derivative node index is not valid. [ mNodeIndex >= "
           "number_of_nodes, mNodeIndex = "
        << mNodeIndex << ", number_of_nodes = " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(mDirectionIndex >= TDim)
        << "Derivative direction index is not valid. [ mDirectionIndex >= "
           "TDim, mDirectionIndex = "
        << mDirectionIndex << ", TDim = " << TDim << " ].\n";

    KRATOS_CATCH("");
}
// template instantiations
template class ElementDataDerivative<2>;
template class ElementDataDerivative<3>;

} // namespace Kratos
