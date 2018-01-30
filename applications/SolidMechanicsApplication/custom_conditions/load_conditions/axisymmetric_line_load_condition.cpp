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
#include "custom_conditions/load_conditions/axisymmetric_line_load_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricLineLoadCondition::AxisymmetricLineLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : LineLoadCondition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricLineLoadCondition::AxisymmetricLineLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LineLoadCondition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
  }

  //************************************************************************************
  //************************************************************************************
  AxisymmetricLineLoadCondition::AxisymmetricLineLoadCondition( AxisymmetricLineLoadCondition const& rOther )
    : LineLoadCondition(rOther)     
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer AxisymmetricLineLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new AxisymmetricLineLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************CLONE*******************************************
  //************************************************************************************
  Condition::Pointer AxisymmetricLineLoadCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    AxisymmetricLineLoadCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    //-----------//      
    return Condition::Pointer( new AxisymmetricLineLoadCondition(NewCondition) );
  }


  //***********************************************************************************
  //***********************************************************************************
  AxisymmetricLineLoadCondition::~AxisymmetricLineLoadCondition()
  {
  }

  //************* GETTING METHODS
  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void AxisymmetricLineLoadCondition::CalculateKinematics(ConditionVariables& rVariables,
						      const double& rPointNumber)
  {
    KRATOS_TRY

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //get first vector of the plane
    rVariables.Tangent1[0] = rVariables.j[rPointNumber](0, 0); // x_1,e
    rVariables.Tangent1[1] = rVariables.j[rPointNumber](1, 0); // x_2,e

    rVariables.Normal[0] = -rVariables.j[rPointNumber](1, 0); //-x_2,e
    rVariables.Normal[1] =  rVariables.j[rPointNumber](0, 0); // x_1,e
    
    //Jacobian to the deformed configuration
    rVariables.Jacobian = norm_2(rVariables.Normal);

    //Compute the unit normal and weighted tangents
    if(rVariables.Jacobian>0){
      rVariables.Normal   /= rVariables.Jacobian;
      rVariables.Tangent1 /= rVariables.Jacobian;
    }

    //Jacobian to the last known configuration
    rVariables.Tangent2[0] = rVariables.J[rPointNumber](0, 0); // x_1,e
    rVariables.Tangent2[1] = rVariables.J[rPointNumber](1, 0); // x_2,e

    rVariables.Jacobian = norm_2(rVariables.Tangent2);

    //Set Shape Functions Values for this integration point
    rVariables.N =row( Ncontainer, rPointNumber);

    //Get geometry size
    rVariables.GeometrySize = GetGeometry().Length();

    //Get external load
    this->CalculateExternalLoad(rVariables);
    
    //Calculate radius
    CalculateRadius ( rVariables.CurrentRadius,  rVariables.ReferenceRadius, rVariables.N);
    
    KRATOS_CATCH( "" )
  }


  //*************************COMPUTE AXYSIMMETRIC RADIUS********************************
  //************************************************************************************

  void AxisymmetricLineLoadCondition::CalculateRadius(double & rCurrentRadius,
						  double & rReferenceRadius,
						  const Vector& rN)


  {

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rCurrentRadius=0;
    rReferenceRadius=0;

    if ( dimension == 2 )
      {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
	    array_1d<double, 3 > & CurrentPosition      = GetGeometry()[i].Coordinates();
	    array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;

            rCurrentRadius   += CurrentPosition[0]*rN[i];
            rReferenceRadius += ReferencePosition[0]*rN[i];
	  }
      }


    if ( dimension == 3 )
      {
        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
      }

    KRATOS_CATCH( "" )
  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  void AxisymmetricLineLoadCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

    //contributions to stiffness matrix calculated on the reference config

    LoadCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

  }


  //************************************************************************************
  //************************************************************************************

  void AxisymmetricLineLoadCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

    //contribution to external forces

    LoadCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

  }

  //***********************************************************************************
  //***********************************************************************************

  int AxisymmetricLineLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = LineLoadCondition::Check(rCurrentProcessInfo);

    return ErrorCode;
    
    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************


  void AxisymmetricLineLoadCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineLoadCondition )
  }

  void AxisymmetricLineLoadCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineLoadCondition )
  }



} // Namespace Kratos.
