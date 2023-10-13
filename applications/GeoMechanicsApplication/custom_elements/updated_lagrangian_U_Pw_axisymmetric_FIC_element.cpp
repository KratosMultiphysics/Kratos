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
#include "custom_elements/updated_lagrangian_U_Pw_axisymmetric_FIC_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwUpdatedLagrangianAxisymmetricFICElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwUpdatedLagrangianAxisymmetricFICElement( NewId,
                                                                       this->GetGeometry().Create( ThisNodes ),
                                                                       pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwUpdatedLagrangianAxisymmetricFICElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           GeometryType::Pointer pGeom,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwUpdatedLagrangianAxisymmetricFICElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianAxisymmetricFICElement<TDim,TNumNodes>::
    CalculateBMatrix(Matrix& rB,
                     const Matrix& GradNpT,
                     const Vector& Np)
{
    KRATOS_TRY

    const double radius = GeoElementUtilities::CalculateRadius(Np, this->GetGeometry());

    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const IndexType index = TDim * i;

        rB( INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X ) = GradNpT( i, INDEX_X );
        rB( INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y ) = GradNpT( i, INDEX_Y );
        rB( INDEX_2D_PLANE_STRAIN_ZZ, index + INDEX_X ) = Np[i]/radius;
        rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X ) = GradNpT( i, INDEX_Y );
        rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y ) = GradNpT( i, INDEX_X );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwUpdatedLagrangianAxisymmetricFICElement<TDim,TNumNodes>::
     CalculateIntegrationCoefficient( const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                                      unsigned int PointNumber,
                                      double detJ)

{
    Vector N;
    N = this->GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double radiusWeight = GeoElementUtilities::CalculateAxisymmetricCircumference(N, this->GetGeometry());

    return IntegrationPoints[PointNumber].Weight() * detJ * radiusWeight;
}

//----------------------------------------------------------------------------------------------------

template class UPwUpdatedLagrangianAxisymmetricFICElement<2,3>;
template class UPwUpdatedLagrangianAxisymmetricFICElement<2,4>;

} // Namespace Kratos
