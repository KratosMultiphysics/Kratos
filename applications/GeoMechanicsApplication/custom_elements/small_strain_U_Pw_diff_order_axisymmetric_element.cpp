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
#include "custom_elements/small_strain_U_Pw_diff_order_axisymmetric_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderAxisymmetricElement::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new SmallStrainUPwDiffOrderAxisymmetricElement( NewId,
                                                                             this->GetGeometry().Create( ThisNodes ),
                                                                             pProperties ) );
}

//----------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderAxisymmetricElement::
    Create(IndexType NewId,
           GeometryType::Pointer pGeom,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new SmallStrainUPwDiffOrderAxisymmetricElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderAxisymmetricElement::
    CalculateBMatrix(Matrix& rB,
                     const Matrix& GradNpT,
                     const Vector& Np)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderAxisymmetricElement::CalculateBMatrix()") << std::endl;

    const double radius = GeoElementUtilities::CalculateRadius(Np, this->GetGeometry());

    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().size();

    for ( unsigned int i = 0; i < NumNodes; ++i ) {
        const unsigned int index = Dim * i;

        rB( INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X ) = GradNpT( i, INDEX_X );
        rB( INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y ) = GradNpT( i, INDEX_Y );
        rB( INDEX_2D_PLANE_STRAIN_ZZ, index + INDEX_X ) = Np[i]/radius;
        rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X ) = GradNpT( i, INDEX_Y );
        rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y ) = GradNpT( i, INDEX_X );
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderAxisymmetricElement::CalculateBMatrix()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
double SmallStrainUPwDiffOrderAxisymmetricElement::
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

} // Namespace Kratos
