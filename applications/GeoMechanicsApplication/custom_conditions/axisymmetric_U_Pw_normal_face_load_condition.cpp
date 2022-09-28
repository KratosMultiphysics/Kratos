// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Vahid Galavi
//


// Application includes
#include "custom_conditions/axisymmetric_U_Pw_normal_face_load_condition.hpp"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer AxisymmetricUPwNormalFaceLoadCondition<TDim,TNumNodes>::
    Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new AxisymmetricUPwNormalFaceLoadCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double AxisymmetricUPwNormalFaceLoadCondition<TDim,TNumNodes>::
    CalculateIntegrationCoefficient( const IndexType PointNumber,
                                     const GeometryType::IntegrationPointsArrayType& IntegrationPoints ) const
{
    Vector N;
    N = this->GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double radiusWeight = GeoElementUtilities::CalculateAxisymmetricCircumference(N, this->GetGeometry());

    return IntegrationPoints[PointNumber].Weight() * radiusWeight;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
GeometryData::IntegrationMethod 
    AxisymmetricUPwNormalFaceLoadCondition<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class AxisymmetricUPwNormalFaceLoadCondition<2,2>;

} // Namespace Kratos.
