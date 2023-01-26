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
#include "custom_conditions/axisymmetric_line_normal_load_2D_diff_order_condition.hpp"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

// Default Constructor
AxisymmetricLineNormalLoad2DDiffOrderCondition::
    AxisymmetricLineNormalLoad2DDiffOrderCondition() : LineNormalLoad2DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------
//Constructor 1
AxisymmetricLineNormalLoad2DDiffOrderCondition::
    AxisymmetricLineNormalLoad2DDiffOrderCondition( IndexType NewId,
                                                    GeometryType::Pointer pGeometry ) :
                                                    LineNormalLoad2DDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------
//Constructor 2
AxisymmetricLineNormalLoad2DDiffOrderCondition::
    AxisymmetricLineNormalLoad2DDiffOrderCondition( IndexType NewId,
                                                    GeometryType::Pointer pGeometry,
                                                    PropertiesType::Pointer pProperties ) :
                                                    LineNormalLoad2DDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------
//Destructor
AxisymmetricLineNormalLoad2DDiffOrderCondition::~AxisymmetricLineNormalLoad2DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Condition::Pointer AxisymmetricLineNormalLoad2DDiffOrderCondition::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new AxisymmetricLineNormalLoad2DDiffOrderCondition(NewId,
                                                                                 GetGeometry().Create(ThisNodes),
                                                                                 pProperties));
}

//----------------------------------------------------------------------------------------
double AxisymmetricLineNormalLoad2DDiffOrderCondition::
    CalculateIntegrationCoefficient(const IndexType PointNumber,
                                    const GeometryType::JacobiansType& JContainer,
                                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    KRATOS_TRY

    Vector N;
    N = this->GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double radiusWeight = GeoElementUtilities::CalculateAxisymmetricCircumference(N, this->GetGeometry());

    return IntegrationPoints[PointNumber].Weight() * radiusWeight;

    KRATOS_CATCH( "" )
}


} // Namespace Kratos.
