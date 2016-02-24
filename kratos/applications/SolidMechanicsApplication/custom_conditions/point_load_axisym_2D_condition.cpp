//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/point_load_axisym_2D_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


//***********************************************************************************
//***********************************************************************************
PointLoadAxisym2DCondition::PointLoadAxisym2DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PointLoad2DCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
PointLoadAxisym2DCondition::PointLoadAxisym2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : PointLoad2DCondition(NewId, pGeometry, pProperties)
{

    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PointLoadAxisym2DCondition::PointLoadAxisym2DCondition( PointLoadAxisym2DCondition const& rOther )
    : PointLoad2DCondition(rOther)     
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer PointLoadAxisym2DCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointLoadAxisym2DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer PointLoadAxisym2DCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return (this->Create( NewId, rThisNodes, pGetProperties() ) );
}


//***********************************************************************************
//***********************************************************************************
PointLoadAxisym2DCondition::~PointLoadAxisym2DCondition()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void PointLoadAxisym2DCondition::CalculateKinematics(GeneralVariables& rVariables,
					       const double& rPointNumber)
{
    KRATOS_TRY

    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius);
    
    rVariables.Jacobian = 1.0;

    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void PointLoadAxisym2DCondition::CalculateRadius(double & rCurrentRadius, double & rReferenceRadius)
{

    KRATOS_TRY

    rCurrentRadius=0;
    rReferenceRadius=0;

    //Displacement from the reference to the current configuration
    array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1);
    array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
    array_1d<double, 3 > & CurrentPosition      = GetGeometry()[0].Coordinates();
    array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;

    rCurrentRadius   = CurrentPosition[0];
    rReferenceRadius = ReferencePosition[0];

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void PointLoadAxisym2DCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

    if( GetProperties()[THICKNESS] > 0 )
      IntegrationWeight /=  GetProperties()[THICKNESS];


    //contributions to stiffness matrix calculated on the reference config

    ForceLoadCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

}


//************************************************************************************
//************************************************************************************

void PointLoadAxisym2DCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
  double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

    if( GetProperties()[THICKNESS] > 0 )
      IntegrationWeight /=  GetProperties()[THICKNESS];

    //contribution to external forces

    ForceLoadCondition::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight );

}

//***********************************************************************************
//***********************************************************************************


int PointLoadAxisym2DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************


void PointLoadAxisym2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointLoad2DCondition )
}

void PointLoadAxisym2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointLoad2DCondition )
}


} // Namespace Kratos.
