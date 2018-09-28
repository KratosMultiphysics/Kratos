//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_conditions/axisym_thermal_contact_domain_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymThermalContactDomainPenalty2DCondition::AxisymThermalContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : ThermalContactDomainPenalty2DCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymThermalContactDomainPenalty2DCondition::AxisymThermalContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : ThermalContactDomainPenalty2DCondition( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymThermalContactDomainPenalty2DCondition::AxisymThermalContactDomainPenalty2DCondition( AxisymThermalContactDomainPenalty2DCondition const& rOther)
    :ThermalContactDomainPenalty2DCondition(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

AxisymThermalContactDomainPenalty2DCondition&  AxisymThermalContactDomainPenalty2DCondition::operator=(AxisymThermalContactDomainPenalty2DCondition const& rOther)
{
    ThermalContactDomainPenalty2DCondition::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer AxisymThermalContactDomainPenalty2DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared<AxisymThermalContactDomainPenalty2DCondition>(NewId, GetGeometry().Create( ThisNodes ), pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer AxisymThermalContactDomainPenalty2DCondition::Clone( IndexType NewId, NodesArrayType const& ThisNodes ) const
{
  return this->Create(NewId, ThisNodes, pGetProperties());
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************


AxisymThermalContactDomainPenalty2DCondition::~AxisymThermalContactDomainPenalty2DCondition()
{
}




//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//*********************************COMPUTE RADIUS*************************************
//************************************************************************************

void AxisymThermalContactDomainPenalty2DCondition::CalculateRadius(double & rCurrentRadius,
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
            //std::cout<<" node "<<i<<" -> DeltaDisplacement : "<<DeltaDisplacement<<std::endl;
        }
    }


    if ( dimension == 3 )
    {
      std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }

    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymThermalContactDomainPenalty2DCondition::CalculateKinematics(GeneralVariables& rVariables,
						    ProcessInfo& rCurrentProcessInfo,
						    const unsigned int& rPointNumber)
{
    KRATOS_TRY

    //reading shape functions
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    //Set Shape Functions Values for this integration point
    Vector N = row( Ncontainer, rPointNumber);

    //Calculate radius
    rVariables.CurrentRadius = 0;
    rVariables.ReferenceRadius = 0;
    CalculateRadius( rVariables.CurrentRadius, rVariables.ReferenceRadius, N );

    //Calculate Current Contact Projections
    this->CalcProjections(rVariables,rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void AxisymThermalContactDomainPenalty2DCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
{

  double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

  if ( GetProperties()[THICKNESS]>0 )
      IntegrationWeight /= GetProperties()[THICKNESS];


  //contributions to stiffness matrix calculated on the reference config
  this->CalculateAndAddThermalKm( rLeftHandSideMatrix, rVariables, IntegrationWeight );

  //KRATOS_WATCH(rLeftHandSideMatrix)
}


//************************************************************************************
//************************************************************************************

void AxisymThermalContactDomainPenalty2DCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, double& rIntegrationWeight)
{
  double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

  if ( GetProperties()[THICKNESS]>0 )
      IntegrationWeight /= GetProperties()[THICKNESS];


  //contribution to contact forces
  this->CalculateAndAddThermalContactForces(rRightHandSideVector, rVariables, IntegrationWeight);

  //KRATOS_WATCH(rRightHandSideVector)
}



//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  AxisymThermalContactDomainPenalty2DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}


void AxisymThermalContactDomainPenalty2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ThermalContactDomainPenalty2DCondition );

}

void AxisymThermalContactDomainPenalty2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ThermalContactDomainPenalty2DCondition );

}



} // Namespace Kratos
