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


// Project includes
#include "custom_conditions/axisymmetric_line_normal_fluid_flux_2D_diff_order_condition.hpp"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

// Default Constructor
AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::
    AxisymmetricLineNormalFluidFlux2DDiffOrderCondition() : LineNormalFluidFlux2DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------
//Constructor 1
AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::
    AxisymmetricLineNormalFluidFlux2DDiffOrderCondition(IndexType NewId,
                                            GeometryType::Pointer pGeometry) :
                                            LineNormalFluidFlux2DDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------
//Constructor 2
AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::
    AxisymmetricLineNormalFluidFlux2DDiffOrderCondition(IndexType NewId,
                                            GeometryType::Pointer pGeometry,
                                            PropertiesType::Pointer pProperties) :
                                            LineNormalFluidFlux2DDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------
//Destructor
AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::~AxisymmetricLineNormalFluidFlux2DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Condition::Pointer AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new AxisymmetricLineNormalFluidFlux2DDiffOrderCondition(NewId,
                                                                          GetGeometry().Create(ThisNodes),
                                                                          pProperties) );
}

//----------------------------------------------------------------------------------------
double AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::
    CalculateIntegrationCoefficient(const IndexType PointNumber,
                                    const GeometryType::JacobiansType& JContainer,
                                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    KRATOS_TRY

    const double dx_dxi = JContainer[PointNumber](0,0);
    const double dy_dxi = JContainer[PointNumber](1,0);

    const double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi);

    Vector N;
    N = this->GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double radiusWeight = GeoElementUtilities::CalculateAxisymmetricCircumference(N, this->GetGeometry());

    return ds * IntegrationPoints[PointNumber].Weight() * radiusWeight;

    KRATOS_CATCH( "" )
}


} // Namespace Kratos.
