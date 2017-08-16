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
#include "custom_conditions/moment_condition.hpp"


namespace Kratos
{

  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition()
    : BoundaryCondition()
  {
  }
  
  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : BoundaryCondition(NewId, pGeometry)
  {
  }

  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : BoundaryCondition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition( MomentCondition const& rOther )
    : BoundaryCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************
  Condition::Pointer MomentCondition::Create(IndexType NewId,
					     NodesArrayType const& ThisNodes,
					     PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new MomentCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }

  
  //************************************CLONE*******************************************
  //************************************************************************************
  
  Condition::Pointer LoadCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    std::cout<<" Call base class MOMENT CONDITION Clone "<<std::endl;
  
    MomentCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());
  
    return Condition::Pointer( new MomentCondition(NewCondition) );
  }
  
  //************************************************************************************
  //************************************************************************************
  MomentCondition::~MomentCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************
  unsigned int MomentCondition::GetDofsSize()
  {
    KRATOS_TRY
     
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();    

    if(dimension == 2)
      {
	return number_of_nodes;
      }
    else
      {
	unsigned int size = number_of_nodes * dimension;
	return size;
      }
    
    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void LoadCondition::CalculateExternalMoment(ConditionVariables& rVariables)
  {
    KRATOS_TRY
      
    KRATOS_ERROR << "calling the base class CalculateExternalMoment method for a moment condition... " << std::endl;
    
    KRATOS_CATCH( "" )
  }
  
  //************************************************************************************
  //************************************************************************************
  void MomentCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    if(rRightHandSideVector.size() != 3)
      rRightHandSideVector.resize(3,false);

    array_1d<double,3>& Moment = GetGeometry()[0].GetSolutionStepValue(POINT_MOMENT);
    rRightHandSideVector[0] = Moment[0];
    rRightHandSideVector[1] = Moment[1];
    rRightHandSideVector[2] = Moment[2];

    //current rotations to compute energy
    if( GetGeometry()[0].SolutionStepsDataHas(ANGULAR_VELOCITY) ){
      array_1d<double,3> Rotation = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) + GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY,1);
      Rotation *= 0.5 * rCurrentProcessInfo[DELTA_TIME];
      
      mEnergy += Moment[0] * Rotation[0] + Moment[1] * Rotation[1] + Moment[2] * Rotation[2];
    }
        
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************
  void MomentCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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

    mEnergy += Moment[0] * Rotation[0] + Moment[1] * Rotation[1] + Moment[2] * Rotation[2];
    
    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************
  void MomentCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
  void MomentCondition::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
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

  void MomentCondition::AddExplicitContribution(const VectorType& rRHS, 
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

  void MomentCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							  std::vector<double>& rValues,
							  const ProcessInfo& rCurrentProcessInfo )
  { 
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
  }

  //************************************************************************************
  //************************************************************************************

  void MomentCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
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

  void MomentCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
  }
  
  void MomentCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
  }


} // Namespace Kratos


