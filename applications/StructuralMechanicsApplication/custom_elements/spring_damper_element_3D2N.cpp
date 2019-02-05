// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/spring_damper_element_3D2N.hpp"

#include "structural_mechanics_application_variables.h"

#define OPT_NUM_NODES 2
#define OPT_NUM_DOFS 12
#define OPT_NUM_DIMS 3

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

SpringDamperElement3D2N::SpringDamperElement3D2N( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SpringDamperElement3D2N::SpringDamperElement3D2N( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SpringDamperElement3D2N::SpringDamperElement3D2N( SpringDamperElement3D2N const& rOther)
    :Element(rOther)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SpringDamperElement3D2N&  SpringDamperElement3D2N::operator=(SpringDamperElement3D2N const& rOther)
{

    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SpringDamperElement3D2N::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_shared<SpringDamperElement3D2N>( NewId, GetGeometry().Create( rThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

Element::Pointer SpringDamperElement3D2N::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_shared<SpringDamperElement3D2N>( NewId, pGeom, pProperties );
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SpringDamperElement3D2N::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    SpringDamperElement3D2N NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Kratos::make_shared<SpringDamperElement3D2N>(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SpringDamperElement3D2N::~SpringDamperElement3D2N()
{
}

//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    //NEEDED TO DEFINE THE DOFS OF THE ELEMENT

    rElementalDofList.resize( 0 );

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( ROTATION_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( ROTATION_Y ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( ROTATION_Z ) );
    }
}

//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    if ( rResult.size() != OPT_NUM_DOFS )
    {
        rResult.resize( OPT_NUM_DOFS, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const unsigned int index = i * 6;
        rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        rResult[index+3] = GetGeometry()[i].GetDof( ROTATION_X ).EquationId();
        rResult[index+4] = GetGeometry()[i].GetDof( ROTATION_Y ).EquationId();
        rResult[index+5] = GetGeometry()[i].GetDof( ROTATION_Z ).EquationId();
    }
}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void SpringDamperElement3D2N::GetValuesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    if ( rValues.size() != OPT_NUM_DOFS )
    {
	    rValues.resize( OPT_NUM_DOFS, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
        const array_1d<double, 3>& rot  = GetGeometry()[i].FastGetSolutionStepValue( ROTATION, Step );

        const unsigned int index = i * 6;
        rValues[index]   = disp[0];
        rValues[index+1] = disp[1];
        rValues[index+2] = disp[2];
        rValues[index+3] = rot[0];
        rValues[index+4] = rot[1];
        rValues[index+5] = rot[2];
    }
}


//************************************VELOCITY****************************************
//************************************************************************************

void SpringDamperElement3D2N::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
    if ( rValues.size() != OPT_NUM_DOFS )
    {
	    rValues.resize( OPT_NUM_DOFS, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY, Step );
        // const array_1d<double, 3>& avel = GetGeometry()[i].FastGetSolutionStepValue( ANGULAR_VELOCITY, Step );
        unsigned int index = i * 6;
        rValues[index]   = vel[0];
        rValues[index+1] = vel[1];
        rValues[index+2] = vel[2];
        // rValues[index+3] = avel[0];
        // rValues[index+4] = avel[1];
        // rValues[index+5] = avel[2];
        rValues[index+3] = 0.0;
        rValues[index+4] = 0.0;
        rValues[index+5] = 0.0;
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void SpringDamperElement3D2N::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
    if ( rValues.size() != OPT_NUM_DOFS )
    {
	    rValues.resize( OPT_NUM_DOFS, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION, Step );
        // const array_1d<double, 3>& aacc = GetGeometry()[i].FastGetSolutionStepValue( ANGULAR_ACCELERATION, Step );
        unsigned int index = i * 6;
        rValues[index]   = acc[0];
        rValues[index+1] = acc[1];
        rValues[index+2] = acc[2];
        // rValues[index+3] = aacc[0];
        // rValues[index+4] = aacc[1];
        // rValues[index+5] = aacc[2];
        rValues[index+3] = 0.0;
        rValues[index+4] = 0.0;
        rValues[index+5] = 0.0;
    }
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::Initialize()
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

////************************************************************************************
////************************************************************************************

void SpringDamperElement3D2N::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

////************************************************************************************
////************************************************************************************
void SpringDamperElement3D2N::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void SpringDamperElement3D2N::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void SpringDamperElement3D2N::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // Explicit case:
    this->ClearNodalForces();

    KRATOS_CATCH( "" );
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
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

void SpringDamperElement3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    if ( rRightHandSideVector.size() != OPT_NUM_DOFS )
    {
        rRightHandSideVector.resize( OPT_NUM_DOFS, false );
    }

    rRightHandSideVector = ZeroVector( OPT_NUM_DOFS ); //resetting RHS

    array_1d<double, OPT_NUM_DOFS > current_displacement = ZeroVector( OPT_NUM_DOFS );
    array_1d<double, 2*OPT_NUM_DIMS > elemental_stiffness = ZeroVector( 2*OPT_NUM_DIMS );
    const array_1d<double, 3> NodalStiffness = Element::GetValue( NODAL_STIFFNESS );
    elemental_stiffness[0] = NodalStiffness[0];
    elemental_stiffness[1] = NodalStiffness[1];
    elemental_stiffness[2] = NodalStiffness[2];
    // TODO: Rotational stiffness not yet implemented..
//     const array_1d<double, 3> NodalRotStiffness = Element::GetValue( NODAL_ROTATIONAL_STIFFNESS );
//     elemental_stiffness[3] = NodalRotStiffness[0];
//     elemental_stiffness[4] = NodalRotStiffness[1];
//     elemental_stiffness[5] = NodalRotStiffness[2];

    const array_1d<double, 3> ddisp =this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT)
                                   - this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);

    current_displacement[0] = -ddisp[0];
    current_displacement[1] = -ddisp[1];
    current_displacement[2] = -ddisp[2];
    current_displacement[6] = ddisp[0];
    current_displacement[7] = ddisp[1];
    current_displacement[8] = ddisp[2];

    for ( std::size_t i = 0; i < OPT_NUM_DOFS; ++i )
    {
        rRightHandSideVector[i]  -= elemental_stiffness[i % 6] * current_displacement[i];
    }

}

//***********************************************************************************
//***********************************************************************************

void SpringDamperElement3D2N::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    // Resizing the LHS
    unsigned int system_size = OPT_NUM_DOFS;

    if ( rLeftHandSideMatrix.size1() != system_size )
    {
        rLeftHandSideMatrix.resize( system_size, system_size, false );
    }

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

    // elemental_stiffness: kx, ky, kz, cpx, cpy, cpz
    array_1d<double, 2*OPT_NUM_DIMS > elemental_stiffness = ZeroVector( 2*OPT_NUM_DIMS );
    const array_1d<double, 3> NodalStiffness = Element::GetValue( NODAL_STIFFNESS );
    elemental_stiffness[0] = NodalStiffness[0];
    elemental_stiffness[1] = NodalStiffness[1];
    elemental_stiffness[2] = NodalStiffness[2];
    // TODO: Rotational stiffness not yet implemented..
//     const array_1d<double, 3> NodalRotStiffness = Element::GetValue( NODAL_ROTATIONAL_STIFFNESS );
//     elemental_stiffness[3] = NodalRotStiffness[0];
//     elemental_stiffness[4] = NodalRotStiffness[1];
//     elemental_stiffness[5] = NodalRotStiffness[2];

    for ( std::size_t i = 0; i < 2*OPT_NUM_DIMS; ++i )
    {
        rLeftHandSideMatrix(i    ,i    ) += elemental_stiffness[i];
        rLeftHandSideMatrix(i + 6,i + 6) += elemental_stiffness[i];
        rLeftHandSideMatrix(i    ,i + 6) -= elemental_stiffness[i];
        rLeftHandSideMatrix(i + 6,i    ) -= elemental_stiffness[i];
    }

}

//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::ClearNodalForces()
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

//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //this is a massless element
    unsigned int system_size = OPT_NUM_DOFS;

    if ( rMassMatrix.size1() != system_size )
    {
        rMassMatrix.resize( system_size, system_size, false );
    }

    rMassMatrix = ZeroMatrix( system_size, system_size );

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    const unsigned int system_size = OPT_NUM_DOFS;

    rDampingMatrix = ZeroMatrix( system_size, system_size );

    if ( Element::Has( NODAL_DAMPING_RATIO ) )
    {
        array_1d<double, 2*OPT_NUM_DIMS> elemental_damping_ratio = ZeroVector( 2*OPT_NUM_DIMS );
        const array_1d<double, 3> NodalDamping = Element::GetValue( NODAL_DAMPING_RATIO );
        elemental_damping_ratio[0] = NodalDamping[0];
        elemental_damping_ratio[1] = NodalDamping[1];
        elemental_damping_ratio[2] = NodalDamping[2];

        for ( std::size_t i = 0; i < 2*OPT_NUM_DIMS; ++i )
        {
            rDampingMatrix(i    , i   ) += elemental_damping_ratio[i];
            rDampingMatrix(i + 6,i + 6) += elemental_damping_ratio[i];
            rDampingMatrix(i    ,i + 6) -= elemental_damping_ratio[i];
            rDampingMatrix(i + 6,i    ) -= elemental_damping_ratio[i];
        }
    }

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

int SpringDamperElement3D2N::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
    {
        KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered)" << std::endl;
    }

    if ( DISPLACEMENT.Key() == 0 )
    {
        KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered)" << std::endl;
    }

    if ( ACCELERATION.Key() == 0 )
    {
        KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered)" << std::endl;
    }

    if ( NODAL_MASS.Key() == 0 )
    {
        KRATOS_ERROR << "NODAL_MASS has Key zero! (check if the application is correctly registered)" << std::endl;
    }

    if ( NODAL_STIFFNESS.Key() == 0 )
    {
        KRATOS_ERROR << "NODAL_STIFFNESS has Key zero! (check if the application is correctly registered)" << std::endl;
    }

    if ( VOLUME_ACCELERATION.Key() == 0 )
    {
        KRATOS_ERROR << "VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered)" << std::endl;
    }

    for ( std::size_t i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( VOLUME_ACCELERATION ) == false )
        {
            KRATOS_ERROR << "Missing variable VOLUME_ACCELERATION on node " << this->GetGeometry()[i].Id() << std::endl;
        }
    }

    // Verify that the dofs exist
    for ( std::size_t i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
        {
            KRATOS_ERROR << "Missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
        }

        if ( this->GetGeometry()[i].SolutionStepsDataHas( ROTATION ) == false )
        {
            KRATOS_ERROR << "Missing variable ROTATION on node " << this->GetGeometry()[i].Id() << std::endl;
        }

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
        {
            KRATOS_ERROR << "Missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
        }

        if ( this->GetGeometry()[i].HasDofFor( ROTATION_X ) == false || this->GetGeometry()[i].HasDofFor( ROTATION_Y ) == false || this->GetGeometry()[i].HasDofFor( ROTATION_Z ) == false )
        {
            KRATOS_ERROR << "Missing one of the dofs for the variable ROTATION on node " << GetGeometry()[i].Id() << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH( "Problem in the Check in the SpringDamperElement3D2N" )
}


//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
}

void SpringDamperElement3D2N::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
}


} // Namespace Kratos


