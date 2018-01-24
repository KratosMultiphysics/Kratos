//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Jordi Cotela
//

#include "fluid_element_utilities.h"

namespace Kratos {

template< std::size_t TNumNodes >
FluidElementUtilities<TNumNodes>::~FluidElementUtilities() {}


template < std::size_t TNumNodes  >
void FluidElementUtilities<TNumNodes>::GetStrainMatrix(
    const ShapeDerivatives2DType& rDNDX,
    boost::numeric::ublas::bounded_matrix<double, VoigtVector2DSize, 3*TNumNodes>& rStrainMatrix) {

    rStrainMatrix.clear();
    for (std::size_t i = 0; i < TNumNodes; i++) {
        const unsigned int col = 3*i;
        rStrainMatrix(0,col)   = rDNDX(i,0);
        rStrainMatrix(1,col+1) = rDNDX(i,1);
        rStrainMatrix(2,col)   = rDNDX(i,1);
        rStrainMatrix(2,col+1) = rDNDX(i,0);
    }
}

template < std::size_t TNumNodes >
void FluidElementUtilities<TNumNodes>::GetStrainMatrix(
    const ShapeDerivatives3DType& rDNDX,
    boost::numeric::ublas::bounded_matrix<double, VoigtVector3DSize, 4*TNumNodes>& rStrainMatrix) {

    rStrainMatrix.clear();
    for (std::size_t i = 0; i < TNumNodes; i++) {
        const unsigned int col = 4*i;
        rStrainMatrix(0,col)   = rDNDX(i,0);
        rStrainMatrix(1,col+1) = rDNDX(i,1);
        rStrainMatrix(2,col+2) = rDNDX(i,2);
        rStrainMatrix(3,col)   = rDNDX(i,1);
        rStrainMatrix(3,col+1) = rDNDX(i,0);
        rStrainMatrix(4,col+1) = rDNDX(i,2);
        rStrainMatrix(4,col+2) = rDNDX(i,1);
        rStrainMatrix(5,col)   = rDNDX(i,2);
        rStrainMatrix(5,col+2) = rDNDX(i,0);
    }
}

template< std::size_t TNumNodes >
void FluidElementUtilities<TNumNodes>::VoigtTransformForProduct(
    const array_1d<double,3>& rVector,
    boost::numeric::ublas::bounded_matrix<double, 2, VoigtVector2DSize>& rVoigtMatrix) {

    rVoigtMatrix.clear();

    rVoigtMatrix(0,0) = rVector(0);
    rVoigtMatrix(0,2) = rVector(1);
    rVoigtMatrix(1,1) = rVector(1);
    rVoigtMatrix(1,2) = rVector(0);
}

template< std::size_t TNumNodes >
void FluidElementUtilities<TNumNodes>::VoigtTransformForProduct(
    const array_1d<double,3>& rVector,
    boost::numeric::ublas::bounded_matrix<double, 3, VoigtVector3DSize>& rVoigtMatrix) {

    rVoigtMatrix.clear();

    rVoigtMatrix(0,0) = rVector(0);
    rVoigtMatrix(0,3) = rVector(1);
    rVoigtMatrix(0,5) = rVector(2);
    rVoigtMatrix(1,1) = rVector(1);
    rVoigtMatrix(1,3) = rVector(0);
    rVoigtMatrix(1,4) = rVector(2);
    rVoigtMatrix(2,2) = rVector(2);
    rVoigtMatrix(2,4) = rVector(1);
    rVoigtMatrix(2,5) = rVector(0);
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class FluidElementUtilities<3>;
template class FluidElementUtilities<4>;

///////////////////////////////////////////////////////////////////////////////////////////////////

}