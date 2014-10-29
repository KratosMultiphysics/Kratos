//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:           September 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/rigid_body_point_rigid_contact_condition.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{



//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{

    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : Condition(NewId, pGeometry, pProperties)
{
    mpRigidWall = pRigidWall;
    
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    
    mMasterElements = GetGeometry()[inode].GetValue(MASTER_ELEMENTS);
    
    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition( RigidBodyPointRigidContactCondition const& rOther )
    : Condition(rOther)      
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer RigidBodyPointRigidContactCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new RigidBodyPointRigidContactCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer RigidBodyPointRigidContactCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return this->Create( NewId, rThisNodes, pGetProperties() );
}


//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::~RigidBodyPointRigidContactCondition()
{
}

//************* GETTING METHODS

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetDofList(DofsVectorType& rConditionDofList,
				    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rConditionDofList.resize(0);

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetDofList(rConditionDofList, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::EquationIdVector(EquationIdVectorType& rResult,
					  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rResult.resize( 0, false );

    Element& MasterElement = mMasterElements.back();
    MasterElement.EquationIdVector(rResult, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetValuesVector(rValues, Step);


    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetFirstDerivativesVector(rValues, Step);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetSecondDerivativesVector(rValues, Step);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void RigidBodyPointRigidContactCondition::ClearNodalForces()
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	GetGeometry()[i].SetLock();
	array_1d<double, 3> & ContactForce  = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
	ContactForce.clear();
	GetGeometry()[i].UnSetLock();
      }

    //KRATOS_WATCH( " CLEAR NODAL FORCE " )

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::AddExplicitContribution(const VectorType& rRHSVector, 
						 const Variable<VectorType>& rRHSVariable, 
						 Variable<array_1d<double,3> >& rDestinationVariable, 
						 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == CONTACT_FORCES_VECTOR && rDestinationVariable == CONTACT_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = i * (dimension * (dimension-1));

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ContactForce[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }


    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = i * (dimension * (dimension-1));

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
}

//************* STARTING - ENDING  METHODS
//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    ClearNodalForces();

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void RigidBodyPointRigidContactCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
  CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
  CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;
  
  ClearNodalForces();
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
  KRATOS_TRY
    
  CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
  CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
  CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;
  
  KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************


int RigidBodyPointRigidContactCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
