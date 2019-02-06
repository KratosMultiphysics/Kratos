//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2013 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/beam_contact/beam_point_pressure_condition.hpp"


namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  BeamPointPressureCondition::BeamPointPressureCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : Condition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  BeamPointPressureCondition::BeamPointPressureCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : Condition(NewId, pGeometry, pProperties)
  {
  }

  //************************************************************************************
  //************************************************************************************
  BeamPointPressureCondition::BeamPointPressureCondition( BeamPointPressureCondition const& rOther )
  : Condition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer BeamPointPressureCondition::Create(IndexType NewId, NodesArrayType
								 const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<BeamPointPressureCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************************************************************
  //************************************************************************************


  BeamPointPressureCondition::~BeamPointPressureCondition()
  {

  }

  //************************************************************************************
  //************************************************************************************

  void BeamPointPressureCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {


  }


  void BeamPointPressureCondition::GetDofList(DofsVectorType& rConditionDofList,
				    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rConditionDofList.resize(0);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

	if( dimension == 3 ){
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	}

        rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));

    }


    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void BeamPointPressureCondition::EquationIdVector(EquationIdVectorType& rResult,
					  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int condition_size        = number_of_nodes * (dimension * (dimension-1));

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++) //TODO: fix this. Apparently, it would not work in 2D!! MA
    {
        int index = i * (dimension * (dimension-1));
        rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();

	if( dimension == 3){
	  rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	  rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
          rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
	}

	rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }

//    std::cout<<" Pressure condition "<<std::endl;
//     for ( unsigned int i = 0; i < number_of_nodes; i++ )
//       {
//     	unsigned int index = i * ( dimension * 2 );
//     	for ( unsigned int j = index; j <= (index+5); j ++ )
//     	  std::cout<<rResult[j]<<" ";
//        std::cout<<std::endl;
//       }

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void BeamPointPressureCondition::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size  = number_of_nodes * (dimension * (dimension-1));

    if ( rValues.size() != condition_size )
      rValues.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++) //TODO: fix this. Apparentñy, it would not work in 2D!! MA
    {
        unsigned int index = i * (dimension * (dimension-1));
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 ){
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
	    rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
	    rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
	}

	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void BeamPointPressureCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * (dimension * (dimension-1));

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ ) //TODO: fix this. Apparentñy, it would not work in 2D!! MA
    {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 ){
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
	    rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
	    rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
	}

	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );

    }

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void BeamPointPressureCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * (dimension * (dimension-1));

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ ) //TODO: fix this. Apparentñy, it would not work in 2D!! MA
    {
        unsigned int index = i * (dimension * (dimension-1));
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 ){
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
	  rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
	  rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );

	}

	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );

    }
    KRATOS_CATCH( "" )
}


  void BeamPointPressureCondition::InitializeSystemMatrices( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector) {

        rLeftHandSideMatrix.resize( 6, 6, false );
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( 6, 6 ); //resetting LHS
        rRightHandSideVector.resize( 6, false );
        rRightHandSideVector = ZeroVector( 6 ); //resetting RHS
  }

  void BeamPointPressureCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector);

        for(unsigned int j = 0; j < dimension; j++) {
	  rRightHandSideVector[j] = mForce[j]; //TODO: This force comes in GLOBAL AXES!!
          rRightHandSideVector[dimension+j] = 0.0;
	}

  }

  void BeamPointPressureCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType LeftHandSideMatrix = Matrix();
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector);

    for(unsigned int j = 0; j < dimension; j++) {
      rRightHandSideVector[j] = mForce[j]; //TODO: This force comes in GLOBAL AXES!!
      rRightHandSideVector[dimension+j] = 0.0;
    }
  }


  void BeamPointPressureCondition::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
{

}


  void BeamPointPressureCondition::AddExplicitContribution(const VectorType& rRHSVector,
						 const Variable<VectorType>& rRHSVariable,
						 Variable<array_1d<double,3> >& rDestinationVariable,
						 const ProcessInfo& rCurrentProcessInfo)
  {

  }


  //************************************************************************************
  //************************************************************************************
  void BeamPointPressureCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
  {

      array_1d<double, 3 > node_to_center_vector;
      node_to_center_vector[0] = node_to_center_vector[1] = node_to_center_vector[2] = 0.0;

      mForce[0] = mForce[1] = mForce[2] = 0.0;

      //Calculate all normals
      for(unsigned int i=0; i<mNodesList.size(); i++){

          node_to_center_vector[0] = GetGeometry()[0].X() - mNodesList[i]->X();
          node_to_center_vector[1] = GetGeometry()[0].Y() - mNodesList[i]->Y();
          node_to_center_vector[2] = GetGeometry()[0].Z() - mNodesList[i]->Z();

          double modulus = sqrt( node_to_center_vector[0]*node_to_center_vector[0] + node_to_center_vector[1]*node_to_center_vector[1] + node_to_center_vector[2]*node_to_center_vector[2] );

          array_1d<double, 3 > unitary_vector;
          noalias(unitary_vector) = node_to_center_vector / modulus;

          //Multiply normal times the pressure and the nodal area
          unitary_vector *= mNodesList[i]->FastGetSolutionStepValue(PRESSURE) * mNodesList[i]->FastGetSolutionStepValue(NODAL_AREA);

          //Add up forces coming from all perimetral nodes
          mForce += unitary_vector;

      }

  }




} // Namespace Kratos
