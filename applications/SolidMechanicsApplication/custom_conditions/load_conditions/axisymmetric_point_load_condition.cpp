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
#include "custom_conditions/load_conditions/axisymmetric_point_load_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricPointLoadCondition::AxisymmetricPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PointLoadCondition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricPointLoadCondition::AxisymmetricPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : PointLoadCondition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  AxisymmetricPointLoadCondition::AxisymmetricPointLoadCondition( AxisymmetricPointLoadCondition const& rOther )
    : PointLoadCondition(rOther)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer AxisymmetricPointLoadCondition::Create(IndexType NewId,
							    NodesArrayType const& ThisNodes,
							    PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<AxisymmetricPointLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer AxisymmetricPointLoadCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    AxisymmetricPointLoadCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    return Kratos::make_shared<AxisymmetricPointLoadCondition>(NewCondition);
  }


  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricPointLoadCondition::~AxisymmetricPointLoadCondition()
  {
  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void AxisymmetricPointLoadCondition::CalculateKinematics(ConditionVariables& rVariables,
							   const double& rPointNumber)
  {
    KRATOS_TRY

    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius);

    rVariables.Jacobian = 1.0;

    this->CalculateExternalLoad(rVariables);

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void AxisymmetricPointLoadCondition::CalculateRadius(double & rCurrentRadius, double & rReferenceRadius)
  {

    KRATOS_TRY

    // rCurrentRadius=0;
    // rReferenceRadius=0;

    // //Displacement from the reference to the current configuration
    // array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    // array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1);
    // array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
    // array_1d<double, 3 > & CurrentPosition      = GetGeometry()[0].Coordinates();
    // array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;

    // rCurrentRadius   = CurrentPosition[0];
    // rReferenceRadius = ReferencePosition[0];

    rCurrentRadius   = GetGeometry()[0].X();
    rReferenceRadius = GetGeometry()[0].X() + GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X,1) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void AxisymmetricPointLoadCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * Globals::Pi * rVariables.CurrentRadius;

    //contributions to stiffness matrix calculated on the reference config

    LoadCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

  }


  //************************************************************************************
  //************************************************************************************

  void AxisymmetricPointLoadCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * Globals::Pi * rVariables.CurrentRadius;

    //contribution to external forces

    LoadCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

  }

  //***********************************************************************************
  //***********************************************************************************


  int AxisymmetricPointLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = PointLoadCondition::Check(rCurrentProcessInfo);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************


  void AxisymmetricPointLoadCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointLoadCondition )
  }

  void AxisymmetricPointLoadCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointLoadCondition )
  }


} // Namespace Kratos.
