// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// Application includes
#include "custom_elements/U_Pw_small_strain_axisymmetric_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainAxisymmetricElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwSmallStrainAxisymmetricElement( NewId,
                                                                    this->GetGeometry().Create( ThisNodes ),
                                                                    pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainAxisymmetricElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           GeometryType::Pointer pGeom,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwSmallStrainAxisymmetricElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainAxisymmetricElement<TDim,TNumNodes>::
    CalculateBMatrix(Matrix& rB,
                     const Matrix& GradNpT,
                     const Vector& Np)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainAxisymmetricElement::CalculateBMatrix()") << std::endl;

    const double radius = GeoElementUtilities::CalculateRadius(Np, this->GetGeometry());

    for ( unsigned int i = 0; i < TNumNodes; ++i )
    {
        const unsigned int index = TDim * i;

        rB( INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X ) = GradNpT( i, INDEX_X );
        rB( INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y ) = GradNpT( i, INDEX_Y );
        rB( INDEX_2D_PLANE_STRAIN_ZZ, index + INDEX_X ) = Np[i]/radius;
        rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X ) = GradNpT( i, INDEX_Y );
        rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y ) = GradNpT( i, INDEX_X );
    }

    // KRATOS_INFO("1-UPwSmallStrainAxisymmetricElement::CalculateBMatrix()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwSmallStrainAxisymmetricElement<TDim,TNumNodes>::
     CalculateIntegrationCoefficient( const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                                      const IndexType& PointNumber,
                                      const double& detJ)

{
    Vector N;
    N = this->GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double radius = GeoElementUtilities::CalculateRadius(N, this->GetGeometry());
    const double radiusWeight = 2.0 * Globals::Pi * radius;

    return IntegrationPoints[PointNumber].Weight() * detJ * radiusWeight;
}

//----------------------------------------------------------------------------------------------------

template class UPwSmallStrainAxisymmetricElement<2,3>;
template class UPwSmallStrainAxisymmetricElement<2,4>;

template class UPwSmallStrainAxisymmetricElement<2,6>;
template class UPwSmallStrainAxisymmetricElement<2,8>;
template class UPwSmallStrainAxisymmetricElement<2,9>;

} // Namespace Kratos
