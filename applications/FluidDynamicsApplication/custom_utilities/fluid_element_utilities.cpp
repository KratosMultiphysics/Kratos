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
    BoundedMatrix<double, VoigtVector2DSize, 3*TNumNodes>& rStrainMatrix) {

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
    BoundedMatrix<double, VoigtVector3DSize, 4*TNumNodes>& rStrainMatrix) {

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

template < std::size_t TNumNodes >
void FluidElementUtilities<TNumNodes>::GetNewtonianConstitutiveMatrix(
    const double DynamicViscosity,
    BoundedMatrix<double, VoigtVector2DSize, VoigtVector2DSize>& rConstitutiveMatrix) {

    constexpr double two_thirds = 2./3.;
    constexpr double four_thirds = 4./3.;

    rConstitutiveMatrix(0,0) = DynamicViscosity * four_thirds;
    rConstitutiveMatrix(0,1) = -DynamicViscosity * two_thirds;
    rConstitutiveMatrix(0,2) = 0.0;
    rConstitutiveMatrix(1,0) = -DynamicViscosity * two_thirds;
    rConstitutiveMatrix(1,1) = DynamicViscosity * four_thirds;
    rConstitutiveMatrix(1,2) = 0.0;
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
    rConstitutiveMatrix(2,2) = DynamicViscosity;
}

template < std::size_t TNumNodes >
void FluidElementUtilities<TNumNodes>::GetNewtonianConstitutiveMatrix(
    const double DynamicViscosity,
    BoundedMatrix<double, VoigtVector3DSize, VoigtVector3DSize>& rConstitutiveMatrix) {

    rConstitutiveMatrix.clear();

    constexpr double two_thirds = 2./3.;
    constexpr double four_thirds = 4./3.;

    rConstitutiveMatrix(0,0) = DynamicViscosity * four_thirds;
    rConstitutiveMatrix(0,1) = -DynamicViscosity * two_thirds;
    rConstitutiveMatrix(0,2) = -DynamicViscosity * two_thirds;

    rConstitutiveMatrix(1,0) = -DynamicViscosity * two_thirds;
    rConstitutiveMatrix(1,1) = DynamicViscosity * four_thirds;
    rConstitutiveMatrix(1,2) = -DynamicViscosity * two_thirds;

    rConstitutiveMatrix(2,0) = -DynamicViscosity * two_thirds;
    rConstitutiveMatrix(2,1) = -DynamicViscosity * two_thirds;
    rConstitutiveMatrix(2,2) = DynamicViscosity * four_thirds;

    rConstitutiveMatrix(3,3) = DynamicViscosity;
    rConstitutiveMatrix(4,4) = DynamicViscosity;
    rConstitutiveMatrix(5,5) = DynamicViscosity;
}


template< std::size_t TNumNodes >
void FluidElementUtilities<TNumNodes>::VoigtTransformForProduct(
    const array_1d<double,3>& rVector,
    BoundedMatrix<double, 2, VoigtVector2DSize>& rVoigtMatrix) {

    rVoigtMatrix.clear();

    rVoigtMatrix(0,0) = rVector(0);
    rVoigtMatrix(0,2) = rVector(1);
    rVoigtMatrix(1,1) = rVector(1);
    rVoigtMatrix(1,2) = rVector(0);
}

template< std::size_t TNumNodes >
void FluidElementUtilities<TNumNodes>::VoigtTransformForProduct(
    const array_1d<double,3>& rVector,
    BoundedMatrix<double, 3, VoigtVector3DSize>& rVoigtMatrix) {

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

template < std::size_t TNumNodes>
void FluidElementUtilities<TNumNodes>::SetNormalProjectionMatrix(
    const array_1d<double, 3>& rUnitNormal,
    BoundedMatrix<double, 2, 2>& rNormalProjMatrix) {

    rNormalProjMatrix.clear();
    rNormalProjMatrix(0,0) = rUnitNormal(0)*rUnitNormal(0);
    rNormalProjMatrix(0,1) = rUnitNormal(0)*rUnitNormal(1);
    rNormalProjMatrix(1,0) = rUnitNormal(1)*rUnitNormal(0);
    rNormalProjMatrix(1,1) = rUnitNormal(1)*rUnitNormal(1);
}

template < std::size_t TNumNodes>
void FluidElementUtilities<TNumNodes>::SetNormalProjectionMatrix(
    const array_1d<double, 3>& rUnitNormal,
    BoundedMatrix<double, 3, 3>& rNormalProjMatrix) {

    noalias(rNormalProjMatrix) = outer_prod(rUnitNormal, rUnitNormal);
}

template < std::size_t TNumNodes>
void FluidElementUtilities<TNumNodes>::SetTangentialProjectionMatrix(
    const array_1d<double, 3>& rUnitNormal,
    BoundedMatrix<double, 2, 2>& rTangProjMatrix) {

    rTangProjMatrix(0,0) = 1.0 - rUnitNormal(0)*rUnitNormal(0);
    rTangProjMatrix(0,1) = - rUnitNormal(0)*rUnitNormal(1);
    rTangProjMatrix(1,0) = - rUnitNormal(1)*rUnitNormal(0);
    rTangProjMatrix(1,1) = 1.0 - rUnitNormal(1)*rUnitNormal(1);
}

template < std::size_t TNumNodes>
void FluidElementUtilities<TNumNodes>::SetTangentialProjectionMatrix(
    const array_1d<double, 3>& rUnitNormal,
    BoundedMatrix<double, 3, 3>& rTangProjMatrix) {

    #ifdef KRATOS_USE_AMATRIX
    BoundedMatrix<double,3,3> id_matrix = IdentityMatrix(3);
    #else
    BoundedMatrix<double,3,3> id_matrix = IdentityMatrix(3,3);
    #endif
    noalias(rTangProjMatrix) = id_matrix - outer_prod(rUnitNormal, rUnitNormal);
}


template < std::size_t TNumNodes>
void FluidElementUtilities<TNumNodes>::DenseSystemSolve(
    const BoundedMatrix<double,2,2> &rA,
    const array_1d<double,2> &rB,
    array_1d<double,2> &rX)
{
    BoundedMatrix<double,2,2> inverse;

    inverse(0,0) =  rA(1,1);
    inverse(0,1) = -rA(0,1);
    inverse(1,0) = -rA(1,0);
    inverse(1,1) =  rA(0,0);

    double det = rA(0,0)*rA(1,1)-rA(0,1)*rA(1,0);
    inverse /= det;

    noalias(rX) = prod(inverse,rB);
}

template < std::size_t TNumNodes>
void FluidElementUtilities<TNumNodes>::DenseSystemSolve(
    const BoundedMatrix<double,3,3> &rA,
    const array_1d<double,3> &rB,
    array_1d<double,3> &rX)
{
    BoundedMatrix<double,3,3> inverse;

    // First column
    inverse(0,0) =  rA(1,1)*rA(2,2) - rA(1,2)*rA(2,1);
    inverse(1,0) = -rA(1,0)*rA(2,2) + rA(1,2)*rA(2,0);
    inverse(2,0) =  rA(1,0)*rA(2,1) - rA(1,1)*rA(2,0);

    // Second column
    inverse(0,1) = -rA(0,1)*rA(2,2) + rA(0,2)*rA(2,1);
    inverse(1,1) =  rA(0,0)*rA(2,2) - rA(0,2)*rA(2,0);
    inverse(2,1) = -rA(0,0)*rA(2,1) + rA(0,1)*rA(2,0);

    // Third column
    inverse(0,2) =  rA(0,1)*rA(1,2) - rA(0,2)*rA(1,1);
    inverse(1,2) = -rA(0,0)*rA(1,2) + rA(0,2)*rA(1,0);
    inverse(2,2) =  rA(0,0)*rA(1,1) - rA(0,1)*rA(1,0);

    double det = rA(0,0)*inverse(0,0) + rA(0,1)*inverse(1,0) + rA(0,2)*inverse(2,0);
    inverse /= det;

    noalias(rX) = prod(inverse,rB);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class FluidElementUtilities<3>; // triangle3
template class FluidElementUtilities<4>; // tetrahedra4, quadrilateral4
template class FluidElementUtilities<8>; // hexahedra8

///////////////////////////////////////////////////////////////////////////////////////////////////

}