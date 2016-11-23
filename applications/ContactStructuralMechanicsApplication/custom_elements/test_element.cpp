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
#include "custom_elements/test_element.hpp"

#include "solid_mechanics_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

TestElement::TestElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TestElement::TestElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
  
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TestElement::TestElement( TestElement const& rOther)
    :Element(rOther)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

TestElement&  TestElement::operator=(TestElement const& rOther)
{

    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer TestElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    //NEEDED TO CREATE AN ELEMENT   
    return Element::Pointer( new TestElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer TestElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    TestElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new TestElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TestElement::~TestElement()
{
}

//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

void TestElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
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

void TestElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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

void TestElement::GetValuesVector( Vector& rValues, int Step )
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

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void TestElement::Initialize()
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

////************************************************************************************
////************************************************************************************

void TestElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

////************************************************************************************
////************************************************************************************
void TestElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void TestElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
//     const array_1d<double, 3 > Current_Displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
//     
//     std::cout << "Solution in iteration: " << rCurrentProcessInfo[NL_ITERATION_NUMBER] << std::endl;
//     std::cout << "x: " << Current_Displacement[0] << std::endl;
//     std::cout << "y: " << Current_Displacement[1] << std::endl;
    
}

////************************************************************************************
////************************************************************************************

void TestElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

//     const array_1d<double, 3 > Current_Displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
//     
//     std::cout << "Final solution: " << std::endl;
//     std::cout << "x: " << Current_Displacement[0] << std::endl;
//     std::cout << "y: " << Current_Displacement[1] << std::endl;
    
    KRATOS_CATCH( "" );
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void TestElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY;

    /* Calculate elemental system */

    // Compute RHS (RHS = rRightHandSideVector)
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Compute LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

//***********************************************************************************
//***********************************************************************************

void TestElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the RHS
    unsigned int system_size = dimension;

    if ( rRightHandSideVector.size() != system_size )
    {
        rRightHandSideVector.resize( system_size, false );
    }

    rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

    const array_1d<double, 3 > Current_Displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);

    if (dimension == 2)
    {
        rRightHandSideVector[0] += 1.0 - Current_Displacement[0] * Current_Displacement[0] - Current_Displacement[1] * Current_Displacement[1];
        rRightHandSideVector[1] += Current_Displacement[0] - Current_Displacement[1];
    }
    else
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "This test element  has just 2 DOF implemented. Dimension = ", dimension);
    }
}

//***********************************************************************************
//***********************************************************************************

void TestElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int system_size = dimension;

    if ( rLeftHandSideMatrix.size1() != system_size )
    {
        rLeftHandSideMatrix.resize( system_size, system_size, false );
    }

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

    const array_1d<double, 3 > Current_Displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    
    if (dimension == 2)
    {
        rLeftHandSideMatrix(0, 0) +=  2.0 * Current_Displacement[0];
        rLeftHandSideMatrix(0, 1) +=  2.0 * Current_Displacement[1];
        rLeftHandSideMatrix(1, 0) += -1.0;
        rLeftHandSideMatrix(1, 1) +=  1.0;
    }
    else
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "This test element  has just 2 DOF implemented. Dimension = ", dimension);
    }
}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& TestElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
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

int TestElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // Verify that the variables are correctly initialized

    if ( DISPLACEMENT.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );
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

void TestElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
}

void TestElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
}


} // Namespace Kratos
