//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/elastic_conditions/axisymmetric_point_elastic_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricPointElasticCondition::AxisymmetricPointElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PointElasticCondition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricPointElasticCondition::AxisymmetricPointElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : PointElasticCondition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  AxisymmetricPointElasticCondition::AxisymmetricPointElasticCondition( AxisymmetricPointElasticCondition const& rOther )
    : PointElasticCondition(rOther)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer AxisymmetricPointElasticCondition::Create(IndexType NewId,
							       NodesArrayType const& ThisNodes,
							       PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<AxisymmetricPointElasticCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer AxisymmetricPointElasticCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    AxisymmetricPointElasticCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    return Kratos::make_shared<AxisymmetricPointElasticCondition>(NewCondition);
  }


  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricPointElasticCondition::~AxisymmetricPointElasticCondition()
  {
  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void AxisymmetricPointElasticCondition::CalculateKinematics(ConditionVariables& rVariables,
							   const double& rPointNumber)
  {
    KRATOS_TRY

    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius);

    rVariables.Jacobian = 1.0;

    this->CalculateExternalStiffness(rVariables);

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void AxisymmetricPointElasticCondition::CalculateRadius(double & rCurrentRadius, double & rReferenceRadius)
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
    rReferenceRadius = GetGeometry()[0].X0();


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void AxisymmetricPointElasticCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * Globals::Pi * rVariables.CurrentRadius;

    //contributions to stiffness matrix calculated on the reference config

    ElasticCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

  }


  //************************************************************************************
  //************************************************************************************

  void AxisymmetricPointElasticCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * Globals::Pi * rVariables.CurrentRadius;

    //contribution to external forces

    ElasticCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

  }

  //***********************************************************************************
  //***********************************************************************************


  int AxisymmetricPointElasticCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = PointElasticCondition::Check(rCurrentProcessInfo);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************


  void AxisymmetricPointElasticCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointElasticCondition )
  }

  void AxisymmetricPointElasticCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointElasticCondition )
  }


} // Namespace Kratos.
