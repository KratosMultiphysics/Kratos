//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/point_torque_3D_condition.hpp"

#include "solid_mechanics_application.h"


namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointTorque3DCondition::PointTorque3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PointTorque3DCondition::PointTorque3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************
PointTorque3DCondition::PointTorque3DCondition( PointTorque3DCondition const& rOther )
    : Condition(rOther)
{
}

//************************************************************************************
//************************************************************************************
Condition::Pointer PointTorque3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointTorque3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//************************************************************************************
//************************************************************************************
PointTorque3DCondition::~PointTorque3DCondition()
{
}


//************************************************************************************
//************************************************************************************
void PointTorque3DCondition::Initialize()
{
    KRATOS_TRY

    mEnergy = 0;

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void PointTorque3DCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    mEnergy = 0;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void PointTorque3DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rRightHandSideVector.size() != 3)
        rRightHandSideVector.resize(3,false);

    array_1d<double,3>& Torque = GetGeometry()[0].GetSolutionStepValue(POINT_TORQUE);
    rRightHandSideVector[0] = Torque[0];
    rRightHandSideVector[1] = Torque[1];
    rRightHandSideVector[2] = Torque[2];

    //current rotations to compute energy
    array_1d<double,3> Rotation = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) + GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY,1);
    Rotation *= 0.5 * rCurrentProcessInfo[DELTA_TIME];

    mEnergy += Torque[0] * Rotation[0] + Torque[1] * Rotation[1] + Torque[2] * Rotation[2];

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void PointTorque3DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 3)
        rLeftHandSideMatrix.resize(3,3,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);

    if(rRightHandSideVector.size() != 3)
        rRightHandSideVector.resize(3,false);

    array_1d<double,3>& Torque = GetGeometry()[0].GetSolutionStepValue(POINT_TORQUE);
    rRightHandSideVector[0] = Torque[0];
    rRightHandSideVector[1] = Torque[1];
    rRightHandSideVector[2] = Torque[2];


    //current rotations to compute energy
    array_1d<double,3> Rotation = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) + GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY,1);
    Rotation *= 0.5 * rCurrentProcessInfo[DELTA_TIME];

    //std::cout<<" Torque "<<Torque<<" Rotation "<<Rotation<<std::endl;

    mEnergy += Torque[0] * Rotation[0] + Torque[1] * Rotation[1] + Torque[2] * Rotation[2];



    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void PointTorque3DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int index = 0;
    const unsigned int dimension = 3;
    rResult.resize(number_of_nodes*dimension);

    for (int i=0; i<number_of_nodes; i++)
    {
        index = i*dimension;
        rResult[index]   = (GetGeometry()[i].GetDof(ROTATION_X).EquationId());
        rResult[index+1] = (GetGeometry()[i].GetDof(ROTATION_Y).EquationId());
        rResult[index+2] = (GetGeometry()[i].GetDof(ROTATION_Z).EquationId());
    }
}

//************************************************************************************
//************************************************************************************
void PointTorque3DCondition::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
    const unsigned int dimension = 3;
    ConditionalDofList.resize(GetGeometry().size()*dimension);
    unsigned int index = 0;
    for (unsigned int i=0; i<GetGeometry().size(); i++)
    {
        index = i*dimension;
        ConditionalDofList[index]   = (GetGeometry()[i].pGetDof(ROTATION_X));
        ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(ROTATION_Y));
        ConditionalDofList[index+2] = (GetGeometry()[i].pGetDof(ROTATION_Z));
    }
}

//***********************************************************************************
//***********************************************************************************

void PointTorque3DCondition::AddExplicitContribution(const VectorType& rRHS, 
						     const Variable<VectorType>& rRHSVariable, 
						     Variable<array_1d<double,3> >& rDestinationVariable, 
						     const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    
    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == MOMENT_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void PointTorque3DCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							  std::vector<double>& rValues,
							  const ProcessInfo& rCurrentProcessInfo )
{ 
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************

void PointTorque3DCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int integration_points_number = 1;

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number, false );

    if ( rVariable == EXTERNAL_ENERGY )
    {
      //reading integration points
      for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
	  rOutput[PointNumber] = mEnergy; //fabs(mEnergy);
	}
    }

    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************

void PointTorque3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
}

void PointTorque3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
}


} // Namespace Kratos


