//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Jordi Cotela
//

#include "vorticity_utilities.h"

namespace Kratos {

template<std::size_t TDim>
VorticityUtilities<TDim>::~VorticityUtilities() {}

template<std::size_t TDim>
double VorticityUtilities<TDim>::CalculateQValue(
    const Geometry<Node<3>>& rGeometry,
    const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients)
{
    return 0.0;
}

template<std::size_t TDim>
double VorticityUtilities<TDim>::CalculateVorticityMagnitude(
    const Geometry<Node<3>>& rGeometry,
    const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients)
{
    return 0.0;
}

template<std::size_t TDim>
void VorticityUtilities<TDim>::CalculateVorticityVector(
    const Geometry<Node<3>>& rGeometry,
    const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients,
    array_1d<double,3>& rVorticity)
{

}

template class VorticityUtilities<2>;
template class VorticityUtilities<3>;

}