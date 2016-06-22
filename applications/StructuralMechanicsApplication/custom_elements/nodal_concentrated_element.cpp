// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/nodal_concentrated_element.hpp"

#include "solid_mechanics_application_variables.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

NodalConcentratedElement::NodalConcentratedElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NodalConcentratedElement::NodalConcentratedElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
  
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NodalConcentratedElement::NodalConcentratedElement( NodalConcentratedElement const& rOther)
    :Element(rOther)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

NodalConcentratedElement&  NodalConcentratedElement::operator=(NodalConcentratedElement const& rOther)
{

    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer NodalConcentratedElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    //NEEDED TO CREATE AN ELEMENT   
    return Element::Pointer( new NodalConcentratedElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer NodalConcentratedElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    NodalConcentratedElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new NodalConcentratedElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NodalConcentratedElement::~NodalConcentratedElement()
{
}

//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    //NEEDED TO DEFINE THE DOFS OF THE ELEMENT 
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    
    rElementalDofList.resize( 0 );

    rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
    rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
    if( dimension == 3 )
    {
	rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Z ) );
    }
    
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if ( rResult.size() != dimension )
    {
        rResult.resize( dimension, false );
    }

    rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();
    rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();
    if( dimension == 3)
    {
	rResult[2] = GetGeometry()[0].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void NodalConcentratedElement::GetValuesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if ( rValues.size() != dimension ) 
    {
	rValues.resize( dimension, false );
    }

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );

    if ( dimension == 3 )
    {
	rValues[2] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}


//************************************VELOCITY****************************************
//************************************************************************************

void NodalConcentratedElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if ( rValues.size() != dimension ) 
    {
      rValues.resize( dimension, false );
    }

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );

    if ( dimension == 3 )
    {
	rValues[2] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Z, Step );
    }

}

//*********************************ACCELERATION***************************************
//************************************************************************************

void NodalConcentratedElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if ( rValues.size() != dimension ) 
    {
	rValues.resize( dimension, false );
    }

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );

    if ( dimension == 3 )
    {
	rValues[2] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Z, Step );
    }
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::Initialize()
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

////************************************************************************************
////************************************************************************************

void NodalConcentratedElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

////************************************************************************************
////************************************************************************************
void NodalConcentratedElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void NodalConcentratedElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void NodalConcentratedElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // Explicit case:
    this->ClearNodalForces();

    KRATOS_CATCH( "" );
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY;

    /* Calculate elemental system */

    // Compute RHS (RHS = rRightHandSideVector = Fext - Fint)
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Compute LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

//***********************************************************************************
//***********************************************************************************

void NodalConcentratedElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the RHS
    unsigned int system_size = dimension;

    if ( rRightHandSideVector.size() != system_size )
    {
        rRightHandSideVector.resize( system_size, false );
    }

    rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

    array_1d<double, 3 > Current_Displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    array_1d<double, 3 > Volume_Acceleration = ZeroVector(3);

    if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) )
    {
        Volume_Acceleration = GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    // Compute and add external forces
    double Nodal_Mass = Element::GetValue(NODAL_MASS);
    for ( unsigned int j = 0; j < dimension; j++ )
    {
        rRightHandSideVector[j]  += Volume_Acceleration[j] * Nodal_Mass;
    }

    // Compute and add internal forces
    array_1d<double, 3 > Nodal_Stiffness = Element::GetValue(NODAL_STIFFNESS);
    for ( unsigned int j = 0; j < dimension; j++ )
    {
        rRightHandSideVector[j]  -= Nodal_Stiffness[j] * Current_Displacement[j];
    }

}

//***********************************************************************************
//***********************************************************************************

void NodalConcentratedElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int system_size = dimension;

    if ( rLeftHandSideMatrix.size1() != system_size )
    {
        rLeftHandSideMatrix.resize( system_size, system_size, false );
    }

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

    array_1d<double, 3 > Nodal_Stiffness = Element::GetValue(NODAL_STIFFNESS);
    for ( unsigned int j = 0; j < dimension; j++ )
    {
        rLeftHandSideMatrix(j, j) += Nodal_Stiffness[j];
    }
}

//***********************************************************************************
//***********************************************************************************

void NodalConcentratedElement::AddExplicitContribution(const VectorType& rRHSVector, 
						    const Variable<VectorType>& rRHSVariable, 
						    Variable<array_1d<double,3> >& rDestinationVariable, 
						    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
    {
	GetGeometry()[0].SetLock();

	array_1d<double, 3 > &ExternalForce = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_FORCE);
	for(unsigned int j = 0; j < dimension; j++)
	{
	    ExternalForce[j] += rRHSVector[j];
	}

	GetGeometry()[0].UnSetLock();
    }

    if( rRHSVariable == INTERNAL_FORCES_VECTOR && rDestinationVariable == INTERNAL_FORCE )
    {
	GetGeometry()[0].SetLock();

	array_1d<double, 3 > &InternalForce = GetGeometry()[0].FastGetSolutionStepValue(INTERNAL_FORCE);
	for(unsigned int j = 0; j < dimension; j++)
	{
	    InternalForce[j] += rRHSVector[j];
	}

	GetGeometry()[0].UnSetLock();
    }

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
    {
	GetGeometry()[0].SetLock();

	array_1d<double, 3 > &ForceResidual = GetGeometry()[0].FastGetSolutionStepValue(FORCE_RESIDUAL);
	for(unsigned int j = 0; j < dimension; j++)
	{
	    ForceResidual[j] += rRHSVector[j];
	}

	GetGeometry()[0].UnSetLock();
    }

    KRATOS_CATCH( "" );

}

//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::ClearNodalForces()
{
    KRATOS_TRY;

    if( GetGeometry()[0].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[0].SolutionStepsDataHas(INTERNAL_FORCE) )
    {
      array_1d<double, 3 > & ExternalForce = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_FORCE);
      array_1d<double, 3 > & InternalForce = GetGeometry()[0].FastGetSolutionStepValue(INTERNAL_FORCE);

      GetGeometry()[0].SetLock();
      ExternalForce.clear();
      InternalForce.clear();
      GetGeometry()[0].UnSetLock();
    }
    
    KRATOS_CATCH( "" );
}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& NodalConcentratedElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY;

    //KRATOS NODAL CURRENT POSITION (X = X0 + DISPLACEMENT_X) IS ALWAYS COMPUTED
    GeometryType& geom = GetGeometry();
    unsigned int dimension = geom.WorkingSpaceDimension();

    rDeltaPosition = ZeroMatrix( 1 , dimension);

    rDeltaPosition(0, 0) = GetGeometry()[0].X() - GetGeometry()[0].X0();
    rDeltaPosition(0, 1) = GetGeometry()[0].Y() - GetGeometry()[0].Y0();
    if(dimension == 3)
    {
	rDeltaPosition(0, 2) = GetGeometry()[0].Z() - GetGeometry()[0].Z0();
    }

    return rDeltaPosition;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int system_size = dimension;

    if ( rMassMatrix.size1() != system_size )
    {
        rMassMatrix.resize( system_size, system_size, false );
    }

    rMassMatrix = ZeroMatrix( system_size, system_size );

    double &Nodal_Mass = Element::GetValue(NODAL_MASS);

    for ( unsigned int j = 0; j < dimension; j++ )
    {
	rMassMatrix( j, j ) = Nodal_Mass;
    }

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    //0.-Initialize the DampingMatrix:
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    const unsigned int system_size = dimension;

    rDampingMatrix = ZeroMatrix( system_size, system_size );

    //1.-Calculate StiffnessMatrix:

    MatrixType StiffnessMatrix     = ZeroMatrix( system_size, system_size );
    VectorType RightHandSideVector = ZeroVector( system_size ); 

    this->CalculateLocalSystem( StiffnessMatrix, RightHandSideVector, rCurrentProcessInfo );

    //2.-Calculate MassMatrix:

    MatrixType MassMatrix  = ZeroMatrix( system_size, system_size );

    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );
    
    //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
    { 
	alpha = GetProperties()[RAYLEIGH_ALPHA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
    { 
	alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
    {
	beta = GetProperties()[RAYLEIGH_BETA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
    { 
	beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    //4.-Compose the Damping Matrix:
   
    //Rayleigh Damping Matrix: alpha*M + beta*K
    MassMatrix      *= alpha;
    StiffnessMatrix *= beta;

    rDampingMatrix  = MassMatrix;
    rDampingMatrix += StiffnessMatrix;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

int NodalConcentratedElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // Verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );
    }

    if ( DISPLACEMENT.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );
    }

    if ( ACCELERATION.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );
    }

    if ( NODAL_MASS.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "NODAL_MASS has Key zero! (check if the application is correctly registered", "" );
    }
    
    if ( NODAL_STIFFNESS.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "NODAL_STIFFNESS has Key zero! (check if the application is correctly registered", "" );
    }

    if ( VOLUME_ACCELERATION.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered", "" );
    }

    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( VOLUME_ACCELERATION ) == false )
	{
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable VOLUME_ACCELERATION on node ", this->GetGeometry()[i].Id() );
	}
    }

    // Verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
	{
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() );
	}

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
	{
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() );
	}
    }

    return 0;

    KRATOS_CATCH( "" );
}


//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
}

void NodalConcentratedElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
}


} // Namespace Kratos


