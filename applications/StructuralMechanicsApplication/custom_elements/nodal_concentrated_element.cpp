// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/nodal_concentrated_element.hpp"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

NodalConcentratedElement::NodalConcentratedElement( 
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    const bool UseRayleighDamping
    )
    : Element( NewId, pGeometry )
    , mUseRayleighDamping( UseRayleighDamping )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NodalConcentratedElement::NodalConcentratedElement( 
    IndexType NewId, GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties, 
    const bool UseRayleighDamping
    )
    : Element( NewId, pGeometry, pProperties )
    , mUseRayleighDamping( UseRayleighDamping )
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NodalConcentratedElement::NodalConcentratedElement( NodalConcentratedElement const& rOther)
    :Element(rOther)
    ,mUseRayleighDamping(rOther.mUseRayleighDamping)
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

Element::Pointer NodalConcentratedElement::Create( 
    IndexType NewId, 
    NodesArrayType const& rThisNodes, 
    PropertiesType::Pointer pProperties 
    ) const
{
    //NEEDED TO CREATE AN ELEMENT   
    return boost::make_shared<NodalConcentratedElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mUseRayleighDamping );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer NodalConcentratedElement::Clone( 
    IndexType NewId, 
    NodesArrayType const& rThisNodes 
    ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    NodalConcentratedElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties(), mUseRayleighDamping );

    return boost::make_shared<NodalConcentratedElement>(new_element);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NodalConcentratedElement::~NodalConcentratedElement()
{
}

//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::GetDofList( 
    DofsVectorType& rElementalDofList, 
    ProcessInfo& rCurrentProcessInfo
    )
{
    //NEEDED TO DEFINE THE DOFS OF THE ELEMENT 
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    rElementalDofList.resize( 0 );

    rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
    rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
    if( dimension == 3 )
        rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Z ) );    
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::EquationIdVector( 
    EquationIdVectorType& rResult, 
    ProcessInfo& rCurrentProcessInfo
    )
{
    //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( rResult.size() != dimension )
        rResult.resize( dimension, false );

    rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();
    rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();
    if( dimension == 3)
        rResult[2] = GetGeometry()[0].GetDof( DISPLACEMENT_Z ).EquationId();
}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void NodalConcentratedElement::GetValuesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( rValues.size() != dimension ) 
        rValues.resize( dimension, false );

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );

    if ( dimension == 3 )
        rValues[2] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z, Step );
}


//************************************VELOCITY****************************************
//************************************************************************************

void NodalConcentratedElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( rValues.size() != dimension ) 
        rValues.resize( dimension, false );

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );

    if ( dimension == 3 )
        rValues[2] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Z, Step );
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void NodalConcentratedElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( rValues.size() != dimension ) 
        rValues.resize( dimension, false );

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );

    if ( dimension == 3 )
        rValues[2] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Z, Step );
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
    const unsigned int system_size = dimension;

    if ( rRightHandSideVector.size() != system_size )
        rRightHandSideVector.resize( system_size, false );

    rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

    const array_1d<double, 3 >& current_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    array_1d<double, 3 > volume_acceleration = ZeroVector(3);

    if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) )
        volume_acceleration = GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);

    // Compute and add external forces
    if (this->Has(NODAL_MASS)){
        const double Nodal_Mass = this->GetValue(NODAL_MASS);
        for ( unsigned int j = 0; j < dimension; ++j )
            rRightHandSideVector[j]  += volume_acceleration[j] * Nodal_Mass;
    }
#ifdef KRATOS_DEBUG
    else
        std::cout << "WARNING :: YOUR NODAL ELEMENT ID: << this->Id() << LACKS OF NODAL MASS" << std::endl;
#endif

    // Compute and add internal forces
    if (this->Has(NODAL_STIFFNESS)){
        const array_1d<double, 3 >& nodal_stiffness = this->GetValue(NODAL_STIFFNESS);
        for ( unsigned int j = 0; j < dimension; ++j )
            rRightHandSideVector[j]  -= nodal_stiffness[j] * current_displacement[j];
    }
#ifdef KRATOS_DEBUG
    else
        std::cout << "WARNING :: YOUR NODAL ELEMENT ID: << this->Id() << LACKS OF NODAL STIFFNESS" << std::endl;
#endif

}

//***********************************************************************************
//***********************************************************************************

void NodalConcentratedElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int system_size = dimension;

    if ( rLeftHandSideMatrix.size1() != system_size )
        rLeftHandSideMatrix.resize( system_size, system_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

    if (this->Has(NODAL_STIFFNESS)){
        const array_1d<double, 3 >& nodal_stiffness = this->GetValue(NODAL_STIFFNESS);
        for ( unsigned int j = 0; j < dimension; ++j )
            rLeftHandSideMatrix(j, j) += nodal_stiffness[j];
    }
#ifdef KRATOS_DEBUG
    else
        std::cout << "WARNING :: YOUR NODAL ELEMENT ID: << this->Id() << LACKS OF NODAL STIFFNESS" << std::endl;
#endif
    
}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& NodalConcentratedElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY;

    //KRATOS NODAL CURRENT POSITION (X = X0 + DISPLACEMENT_X) IS ALWAYS COMPUTED
    GeometryType& geom = GetGeometry();
    const unsigned int dimension = geom.WorkingSpaceDimension();

    rDeltaPosition = ZeroMatrix( 1 , dimension);

    rDeltaPosition(0, 0) = GetGeometry()[0].X() - GetGeometry()[0].X0();
    rDeltaPosition(0, 1) = GetGeometry()[0].Y() - GetGeometry()[0].Y0();
    if(dimension == 3)
        rDeltaPosition(0, 2) = GetGeometry()[0].Z() - GetGeometry()[0].Z0();

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
        rMassMatrix.resize( system_size, system_size, false );

    rMassMatrix = ZeroMatrix( system_size, system_size );

    const double Nodal_Mass = this->GetValue(NODAL_MASS);

    for ( unsigned int j = 0; j < dimension; ++j )
        rMassMatrix( j, j ) = Nodal_Mass;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::CalculateDampingMatrix( 
    MatrixType& rDampingMatrix, 
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    //0.-Initialize the DampingMatrix:
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int system_size = dimension;

    rDampingMatrix = ZeroMatrix( system_size, system_size );

    //Check, if Rayleigh damping is available; use nodal damping, if not
    if( mUseRayleighDamping ) {
        //1.-Calculate StiffnessMatrix:

        MatrixType stiffness_matrix     = ZeroMatrix( system_size, system_size );
        VectorType right_hand_side_vector = ZeroVector( system_size ); 

        this->CalculateLocalSystem( stiffness_matrix, right_hand_side_vector, rCurrentProcessInfo );

        //2.-Calculate MassMatrix:

        MatrixType mass_matrix  = ZeroMatrix( system_size, system_size );

        this->CalculateMassMatrix ( mass_matrix, rCurrentProcessInfo );
        
        //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        double alpha = 0.0;
        if( GetProperties().Has(RAYLEIGH_ALPHA) )
            alpha = GetProperties()[RAYLEIGH_ALPHA];
        else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
            alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

        double beta  = 0.0;
        if( GetProperties().Has(RAYLEIGH_BETA) )
            beta = GetProperties()[RAYLEIGH_BETA];
        else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
            beta = rCurrentProcessInfo[RAYLEIGH_BETA];

        //4.-Compose the Damping Matrix:
    
        //Rayleigh Damping Matrix: alpha*M + beta*K
        mass_matrix      *= alpha;
        stiffness_matrix *= beta;

        rDampingMatrix  = mass_matrix;
        rDampingMatrix += stiffness_matrix;
    } else {
        const array_1d<double, 3 >& nodal_damping_ratio = this->GetValue(NODAL_DAMPING_RATIO);
        for ( unsigned int j = 0; j < dimension; ++j )
            rDampingMatrix(j, j) += nodal_damping_ratio[j];
    }

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

int NodalConcentratedElement::Check( const ProcessInfo& rCurrentProcessInfo )
{    
    KRATOS_TRY
    
    // Base class 
    Element::Check(rCurrentProcessInfo);
    
    // Verify that the variables are correctly initialized
    KRATOS_ERROR_IF(VELOCITY.Key() == 0) << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF(DISPLACEMENT.Key() == 0 ) << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF(ACCELERATION.Key() == 0) << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF(NODAL_MASS.Key() == 0) << "NODAL_MASS has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF(NODAL_STIFFNESS.Key() == 0) << "NODAL_STIFFNESS has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF(NODAL_DAMPING_RATIO.Key() == 0) << "NODAL_DAMPING_RATIO has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF(VOLUME_ACCELERATION.Key() == 0) << "VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;

    for ( std::size_t i = 0; i < this->GetGeometry().size(); ++i ) {
        KRATOS_ERROR_IF(this->GetGeometry()[i].SolutionStepsDataHas( VOLUME_ACCELERATION ) == false) << "Missing variable VOLUME_ACCELERATION on node " << this->GetGeometry()[i].Id() << std::endl;
        // Verify that the dofs exist
        KRATOS_ERROR_IF(this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false) << "Missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
        KRATOS_ERROR_IF(( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )) << "Missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
    }
    
    return 0;

    KRATOS_CATCH( "Problem in the Check in the NodalConcentratedElement" )
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


