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
#include "custom_conditions/point_moment_3D_condition.hpp"

#include "structural_mechanics_application_variables.h"


namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointMoment3DCondition::PointMoment3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PointMoment3DCondition::PointMoment3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************
PointMoment3DCondition::PointMoment3DCondition( PointMoment3DCondition const& rOther )
    : Condition(rOther)
{
}

//************************************************************************************
//************************************************************************************
Condition::Pointer PointMoment3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointMoment3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//************************************************************************************
//************************************************************************************
PointMoment3DCondition::~PointMoment3DCondition()
{
}


//************************************************************************************
//************************************************************************************
void PointMoment3DCondition::Initialize()
{
    KRATOS_TRY

    mEnergy = 0;

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void PointMoment3DCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    mEnergy = 0;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void PointMoment3DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rRightHandSideVector.size() != 3)
        rRightHandSideVector.resize(3,false);

    array_1d<double,3>& Moment = GetGeometry()[0].GetSolutionStepValue(POINT_MOMENT);
    rRightHandSideVector[0] = Moment[0];
    rRightHandSideVector[1] = Moment[1];
    rRightHandSideVector[2] = Moment[2];

    //current rotations to compute energy
    array_1d<double,3> Rotation = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) + GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY,1);
    Rotation *= 0.5 * rCurrentProcessInfo[DELTA_TIME];

    mEnergy += Moment[0] * Rotation[0] + Moment[1] * Rotation[1] + Moment[2] * Rotation[2];

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void PointMoment3DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 3)
        rLeftHandSideMatrix.resize(3,3,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);

    if(rRightHandSideVector.size() != 3)
        rRightHandSideVector.resize(3,false);

    array_1d<double,3>& Moment = GetGeometry()[0].GetSolutionStepValue(POINT_MOMENT);
    rRightHandSideVector[0] = Moment[0];
    rRightHandSideVector[1] = Moment[1];
    rRightHandSideVector[2] = Moment[2];


    //current rotations to compute energy
    array_1d<double,3> Rotation = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) + GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY,1);
    Rotation *= 0.5 * rCurrentProcessInfo[DELTA_TIME];

    //std::cout<<" Moment "<<Moment<<" Rotation "<<Rotation<<std::endl;

    mEnergy += Moment[0] * Rotation[0] + Moment[1] * Rotation[1] + Moment[2] * Rotation[2];



    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void PointMoment3DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
void PointMoment3DCondition::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
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

void PointMoment3DCondition::AddExplicitContribution(const VectorType& rRHS, 
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

	    array_1d<double, 3 > &MomentResidual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		MomentResidual[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void PointMoment3DCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							  std::vector<double>& rValues,
							  const ProcessInfo& rCurrentProcessInfo )
{ 
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************

void PointMoment3DCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
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

void PointMoment3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
}

void PointMoment3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
}


} // Namespace Kratos


