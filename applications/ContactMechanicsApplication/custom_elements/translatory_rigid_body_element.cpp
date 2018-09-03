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
#include "custom_elements/translatory_rigid_body_element.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodyElement::TranslatoryRigidBodyElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : RigidBodyElement(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodyElement::TranslatoryRigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : RigidBodyElement(NewId, pGeometry, pProperties)
{
    KRATOS_TRY

    //DO NOT ADD DOFS HERE!!!

    KRATOS_CATCH( "" )

}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodyElement::TranslatoryRigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, NodesContainerType::Pointer pNodes)
  : RigidBodyElement(NewId, pGeometry, pProperties, pNodes)
{
    KRATOS_TRY

    //DO NOT ADD DOFS HERE!!!

    KRATOS_CATCH( "" )
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TranslatoryRigidBodyElement::TranslatoryRigidBodyElement(TranslatoryRigidBodyElement const& rOther)
    :RigidBodyElement(rOther)
{
}

//*********************************CREATE*********************************************
//************************************************************************************

Element::Pointer TranslatoryRigidBodyElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<TranslatoryRigidBodyElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//*********************************CLONE**********************************************
//************************************************************************************

Element::Pointer TranslatoryRigidBodyElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{

  TranslatoryRigidBodyElement NewElement( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpNodes );

  NewElement.mInitialLocalQuaternion = this->mInitialLocalQuaternion;
  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());

  return Kratos::make_shared<TranslatoryRigidBodyElement>(NewElement);

}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodyElement::~TranslatoryRigidBodyElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************


void TranslatoryRigidBodyElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{

    ElementalDofList.resize(0);

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	if( dimension ==3 )
	  ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
      }

}

//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodyElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * ( dimension );

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension );
      rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
      rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
      if( dimension ==3 )
	rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }

}


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void TranslatoryRigidBodyElement::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * ( dimension );

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
      if( dimension ==3 )
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }

    KRATOS_CATCH( "" )
}

//************************************VELOCITY****************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************
void TranslatoryRigidBodyElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * ( dimension );

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
      if( dimension ==3 )
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }


    KRATOS_CATCH( "" )
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void TranslatoryRigidBodyElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * ( dimension );

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
      if( dimension ==3 )
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodyElement::Initialize()
{
    KRATOS_TRY

    Matrix InitialLocalMatrix = IdentityMatrix(3);

    mInitialLocalQuaternion  = QuaternionType::FromRotationMatrix( InitialLocalMatrix );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodyElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
						VectorType& rRightHandSideVector,
						Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * ( dimension );

    if ( rCalculationFlags.Is(TranslatoryRigidBodyElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( rCalculationFlags.Is(TranslatoryRigidBodyElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );

	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    }
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration
void TranslatoryRigidBodyElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension );

    if(rLeftHandSideMatrix.size1() != MatSize)
      rLeftHandSideMatrix.resize (MatSize, MatSize, false);

    rLeftHandSideMatrix = ZeroMatrix( MatSize, MatSize );

    //rCurrentProcessInfo must give it:
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    double AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];

    double Newmark1 = (1.0/ ( DeltaTime * DeltaTime * rCurrentProcessInfo[NEWMARK_BETA] ));

    //block m(1,1) of the mass matrix

    MatrixType m11 = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass = rVariables.RigidBody.Mass;

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    Matrix DiagonalMatrix = IdentityMatrix(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	m11 = ZeroMatrix(3,3);

    	RowIndex = i * (dimension);

    	for ( unsigned int j = 0; j < number_of_nodes; j++ )
    	  {

    	    ColIndex = j * (dimension);

    	    m11 = (1.0-AlphaM) * Newmark1 * TotalMass * DiagonalMatrix;

    	    //Building the Local Tangent Inertia Matrix
    	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m11, RowIndex, ColIndex );

    	  }
      }


    //std::cout<<" rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration
void TranslatoryRigidBodyElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension );

    if(rRightHandSideVector.size() != MatSize)
      rRightHandSideVector.resize(MatSize, false);

    rRightHandSideVector = ZeroVector( MatSize );

     double TotalMass = 0;
    TotalMass = rVariables.RigidBody.Mass;

    Vector LinearAccelerationVector          = ZeroVector(3);
    Vector CurrentLinearAccelerationVector   = ZeroVector(3);
    Vector PreviousLinearAccelerationVector  = ZeroVector(3);

    Vector CurrentValueVector = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//Current Linear Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ACCELERATION, CurrentValueVector, i );
	CurrentLinearAccelerationVector  = CurrentValueVector;

	//Previous Linear Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ACCELERATION, CurrentValueVector, i );
	PreviousLinearAccelerationVector = CurrentValueVector;

      }

    //Set step variables to local frame (current Frame is the local frame)
    CurrentLinearAccelerationVector     = MapToInitialLocalFrame( CurrentLinearAccelerationVector );
    PreviousLinearAccelerationVector    = MapToInitialLocalFrame( PreviousLinearAccelerationVector );

    double AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    LinearAccelerationVector  = (1.0-AlphaM) * CurrentLinearAccelerationVector + AlphaM * (PreviousLinearAccelerationVector);

    //-----------------
    //block m(1) of the inertial force vector

    //Compute Linear Term:

    Vector LinearInertialForceVector = ZeroVector(3);

    //this transformation is wrong
    //LinearInertialForceVector  = MapToMaterialFrame( TotalQuaternion, LinearInertialForceVector );

    LinearInertialForceVector = TotalMass * LinearAccelerationVector;


    //Initialize Local Matrices
    VectorType Fi = ZeroVector(6);
    unsigned int RowIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

    	RowIndex = i * (dimension);

    	Fi = ZeroVector(6);

    	//nodal force vector
    	Fi  = LinearInertialForceVector;

	BeamMathUtilsType::AddVector(Fi, rRightHandSideVector, RowIndex);

    	//std::cout<<" Fi "<<Fi<<std::endl;

      }

    //std::cout<<" Rigid Body: rRightHandSideVector "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodyElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension );

    if(rMassMatrix.size1() != MatSize)
        rMassMatrix.resize (MatSize, MatSize, false);

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    // Rigid Body Properties
    RigidBodyProperties RigidBody;
    this->CalculateRigidBodyProperties(RigidBody);

    //block m(1,1) of the mass matrix

    MatrixType m11 = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass = RigidBody.Mass;



    for( unsigned int i=0; i < number_of_nodes; i++ )
      {

        double temp = TotalMass;

    	int RowIndex = i * (dimension);

    	for( unsigned int k=0; k < dimension; k++ )
    	  {
    	    m11(k,k) = temp;
    	  }


    	//Building the Local Tangent Inertia Matrix
    	BeamMathUtilsType::AddMatrix( rMassMatrix, m11, RowIndex, RowIndex );

      }


    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodyElement::UpdateRigidBodyNodes(ProcessInfo& rCurrentProcessInfo)
{

     KRATOS_TRY

     Node<3>& rCenterOfGravity = this->GetGeometry()[0];

     array_1d<double, 3 >&  Displacement = rCenterOfGravity.FastGetSolutionStepValue(DISPLACEMENT);
     array_1d<double, 3 >&  Velocity     = rCenterOfGravity.FastGetSolutionStepValue(VELOCITY);
     array_1d<double, 3 >&  Acceleration = rCenterOfGravity.FastGetSolutionStepValue(ACCELERATION);

     for (NodesContainerType::iterator i = mpNodes->begin(); i != mpNodes->end(); ++i)
       {
	 (i)->FastGetSolutionStepValue(DISPLACEMENT) = Displacement;
	 (i)->FastGetSolutionStepValue(VELOCITY)     = Velocity;
	 (i)->FastGetSolutionStepValue(ACCELERATION) = Acceleration;
       }

     KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodyElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, RigidBodyElement )
}

void TranslatoryRigidBodyElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, RigidBodyElement )
}


} // Namespace Kratos
