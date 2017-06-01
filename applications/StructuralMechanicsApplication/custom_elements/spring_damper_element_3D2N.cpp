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

#include "solid_mechanics_application_variables.h"
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
    return Element::Pointer( new SpringDamperElement3D2N( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SpringDamperElement3D2N::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    SpringDamperElement3D2N NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new SpringDamperElement3D2N(NewElement) );
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

    for ( size_t i = 0; i < GetGeometry().size(); ++i)
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

    for ( size_t i = 0; i < GetGeometry().size(); ++i)
    {
        unsigned int index = i * 6;
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

    for ( size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
        const array_1d<double, 3>& rot  = GetGeometry()[i].FastGetSolutionStepValue( ROTATION, Step );

        unsigned int index = i * 6;
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

    for ( size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY, Step );
        std::cout << "TODO: angular velocity??" << std::endl;
        unsigned int index = i * 6;
        rValues[index]   = vel[0];
        rValues[index+1] = vel[1];
        rValues[index+2] = vel[2];
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

    for ( size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION, Step );
        std::cout << "TODO: angular acceleration??" << std::endl;
        unsigned int index = i * 6;
        rValues[index]   = acc[0];
        rValues[index+1] = acc[1];
        rValues[index+2] = acc[2];
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
    //sowas in der art
// KRATOS_TRY
// 		this->mArea = this->GetProperties()[CROSS_AREA];
// 		this->mYoungsModulus = this->GetProperties()[YOUNG_MODULUS];
// 		this->mLength = this->CalculateReferenceLength();
// 		this->mCurrentLength = this->CalculateCurrentLength();
// 		this->mDensity = this->GetProperties()[DENSITY];

// 		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == false) {
// 			this->mPreStress = 0.00;
// 		}
// 		else this->mPreStress = this->GetProperties()[TRUSS_PRESTRESS_PK2];

// 		if (this->GetProperties().Has(TRUSS_IS_CABLE) == false) {
// 			this->mIsCable = false;
// 		}
// 		else this->mIsCable = this->GetProperties()[TRUSS_IS_CABLE];

// 		if (this->mLength == 0.00) {
// 			KRATOS_ERROR << ("Zero length found in element #", this->Id()) <<
// 				std::endl;
// 		}

// KRATOS_CATCH("")
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

void SpringDamperElement3D2N::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY;

    /* Calculate elemental system */
    std::cout << "=====================================================" << std::endl;
    std::cout << "DISCRETE SPRING/DAMPER ELEMENT" << std::endl;
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
    std::cout << "CalculateRightHandSide for id=" << this->Id() << std::endl;

    if ( rRightHandSideVector.size() != OPT_NUM_DOFS )
    {
        rRightHandSideVector.resize( OPT_NUM_DOFS, false );
    }

    rRightHandSideVector = ZeroVector( OPT_NUM_DOFS ); //resetting RHS

    // array_1d<double, OPT_NUM_DOFS > current_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    array_1d<double, OPT_NUM_DOFS > current_displacement = ZeroVector( OPT_NUM_DOFS );
    array_1d<double, OPT_NUM_DOFS > external_forces = ZeroVector( OPT_NUM_DOFS );
    array_1d<double, OPT_NUM_NODES > disp = ZeroVector( OPT_NUM_NODES );
    // array_1d<double, OPT_NUM_NODES> nodal_mass = ZeroVector( OPT_NUM_NODES );
    for ( size_t i = 0; i < OPT_NUM_NODES; ++i )
    {
        const unsigned int index = i * 6;
        const NodeType& iNode = GetGeometry()[i];

        current_displacement[index]   = iNode.FastGetSolutionStepValue( DISPLACEMENT_X );
        current_displacement[index+1] = iNode.FastGetSolutionStepValue( DISPLACEMENT_Y );
        current_displacement[index+2] = iNode.FastGetSolutionStepValue( DISPLACEMENT_Z );

        // if ( iNode.SolutionStepDataHas( NODAL_STIFFNESS ) )
        // {
        //     elemental_stiffness[index] = iNode.FastGetSolutionStepValue( NODAL_STIFFNESS_X );
        //     elemental_stiffness[index+1] = iNode.FastGetSolutionStepValue( NODAL_STIFFNESS_Y );
        //     elemental_stiffness[index+2] = iNode.FastGetSolutionStepValue( NODAL_STIFFNESS_Z );
        // }

        if ( iNode.SolutionStepsDataHas( VOLUME_ACCELERATION ) && iNode.Has( NODAL_MASS ) )
        {
            // external_forces[index]   = iNode.FastGetSolutionStepValue( VOLUME_ACCELERATION_X ) * iNode.GetValue( NODAL_MASS );
            // external_forces[index+1] = iNode.FastGetSolutionStepValue( VOLUME_ACCELERATION_Y ) * iNode.GetValue( NODAL_MASS );
            // external_forces[index+2] = iNode.FastGetSolutionStepValue( VOLUME_ACCELERATION_Z ) * iNode.GetValue( NODAL_MASS );
        }

        // if ( iNode.Has( NODAL_MASS ) )
        // {
        //     nodal_mass[i] = iNode.GetValue( NODAL_MASS );
        // }
        // disp[i] = this->CalculateCurrentLength();
    }

    KRATOS_WATCH(external_forces);
    KRATOS_WATCH(current_displacement);
    // KRATOS_WATCH(nodal_mass);

    // Compute and add external forces
    // double Nodal_Mass = Element::GetValue(NODAL_MASS);
    // for ( size_t i = 0; i < OPT_NUM_NODES; ++i )
    // {
    //     for ( size_t j = 0; j < )
    //     const unsigned int index = i * 6;
    //     rRightHandSideVector[i] += external_forces[i] * nodal_mass[i/OPT_NUM_DIMS];
    //     // rRightHandSideVector[j]  += external_forces[j] * Nodal_Mass;
    // }
    rRightHandSideVector += external_forces;
    // noalias( rRightHandSideVector ) += external_forces;

    // Compute and add internal forces
    // Calculate current length of the element
    // const NodeType& node_0 = GetGeometry()[0];
    // const NodeType& node_1 = GetGeometry()[1];
    // const double lx = node_0.X() - node_1.X();
    // const double ly = node_0.Y() - node_1.Y();
    // const double lz = node_0.Z() - node_1.Z();
    // const double current_length = sqrt( lx*lx + ly*ly + lz*lz );

    // array_1d<double, 3 > elemental_stiffness = Element::GetValue(NODAL_STIFFNESS);

    array_1d<double, 2*OPT_NUM_DIMS > elemental_stiffness = ZeroVector( 2*OPT_NUM_DIMS ); 
    elemental_stiffness[0] = Element::GetValue( NODAL_STIFFNESS_X );
    elemental_stiffness[1] = Element::GetValue( NODAL_STIFFNESS_Y );
    elemental_stiffness[2] = Element::GetValue( NODAL_STIFFNESS_Z );
    std::cout << "ROTATIONAL STIFFNESS MISSING" << std::endl;

    const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
        - this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
        - this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
        - this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);

    std::cout << "du/dv/dw=" << du <<"/" << dv << "/" << dw << std::endl;

    current_displacement[0] = -du;
    current_displacement[1] = -dv;
    current_displacement[2] = -dw;
    current_displacement[6] = du;
    current_displacement[7] = dv;
    current_displacement[8] = dw;

    for ( size_t i = 0; i < OPT_NUM_DOFS; ++i )
    {
        rRightHandSideVector[i]  -= elemental_stiffness[i % 6] * current_displacement[i];
    }
    KRATOS_WATCH(rRightHandSideVector);

}

//***********************************************************************************
//***********************************************************************************

void SpringDamperElement3D2N::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    std::cout << "CalculateLeftHandSide" << this->Id() << std::endl;
    // const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int system_size = OPT_NUM_DOFS;

    if ( rLeftHandSideMatrix.size1() != system_size )
    {
        rLeftHandSideMatrix.resize( system_size, system_size, false );
    }

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

    // array_1d<double, 3 > Nodal_Stiffness = Element::GetValue(NODAL_STIFFNESS);
    // for ( unsigned int j = 0; j < dimension; j++ )
    // {
    //     rLeftHandSideMatrix(j, j) += Nodal_Stiffness[j];
    // }

    //elemental_stiffness: kx, ky, kz, cpx, cpy, cpz
    array_1d<double, 2*OPT_NUM_DIMS > elemental_stiffness = ZeroVector( 2*OPT_NUM_DIMS ); 
    elemental_stiffness[0] = Element::GetValue( NODAL_STIFFNESS_X );
    elemental_stiffness[1] = Element::GetValue( NODAL_STIFFNESS_Y );
    elemental_stiffness[2] = Element::GetValue( NODAL_STIFFNESS_Z );
    std::cout << "ROTATIONAL STIFFNESS MISSING" << std::endl;

    // for ( size_t i = 0; i < OPT_NUM_DOFS; ++i )
    // {
    //     rLeftHandSideMatrix( i, i )  += elemental_stiffness[i % 6];
    // }
    for ( size_t i = 0; i < 2*OPT_NUM_DIMS; ++i )
    {
        rLeftHandSideMatrix( i, i) += elemental_stiffness[i];
        rLeftHandSideMatrix( i+6, i+6) += elemental_stiffness[i];
        rLeftHandSideMatrix( i, i+6) -= elemental_stiffness[i];
        rLeftHandSideMatrix( i+6, i) -= elemental_stiffness[i];
    }

    KRATOS_WATCH(rLeftHandSideMatrix);
}

//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::ClearNodalForces()
{
    KRATOS_TRY;
    std::cout << "ClearNodalForces" << std::endl;
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


Matrix& SpringDamperElement3D2N::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY;
    std::cout << "CalculateDeltaPosition" << std::endl;
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

void SpringDamperElement3D2N::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    std::cout << "CalculateMassMatrix" << std::endl;
    //lumped
    // unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int system_size = OPT_NUM_DOFS;

    if ( rMassMatrix.size1() != system_size )
    {
        rMassMatrix.resize( system_size, system_size, false );
    }

    rMassMatrix = ZeroMatrix( system_size, system_size );

    // KRATOS_WATCH(Element::GetValue(NODAL_MASS));
    // double &Nodal_Mass = Element::GetValue(NODAL_MASS);

    // array_1d<double, OPT_NUM_DOFS> nodal_mass = ZeroVector( OPT_NUM_DOFS );

    for ( size_t i = 0; i < OPT_NUM_NODES; ++i )
    {
        const unsigned int index = i * 6;
        // const NodeType& iNode = GetGeometry()[i];
        const double& nodal_mass = GetGeometry()[i].GetValue( NODAL_MASS );

        rMassMatrix( index, index ) = nodal_mass;
        rMassMatrix( index+1, index+1 ) = nodal_mass;
        rMassMatrix( index+2, index+2 ) = nodal_mass;
        // for ( size_t j = i*6; j < i*6+dimension; ++j)
        // {
        //     rMassMatrix( j, j ) = i_nodal_mass;
        // }
        // nodal_mass[index]   = i_nodal_mass;
        // nodal_mass[index+1] = i_nodal_mass;
        // nodal_mass[index+2] = i_nodal_mass;
        
    }

    // for ( unsigned int j = 0; j < dimension; j++ )
    // {
	//     rMassMatrix( j, j ) = Nodal_Mass;
    // }
    KRATOS_WATCH(rMassMatrix);

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void SpringDamperElement3D2N::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    std::cout << "CalculateDampingMatrix" << std::endl;

    //0.-Initialize the DampingMatrix:
    // const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    const unsigned int system_size = OPT_NUM_DOFS;

    rDampingMatrix = ZeroMatrix( system_size, system_size );

    //1.-Calculate StiffnessMatrix:

    // MatrixType StiffnessMatrix     = ZeroMatrix( system_size, system_size );
    // VectorType RightHandSideVector = ZeroVector( system_size ); 

    // this->CalculateLocalSystem( StiffnessMatrix, RightHandSideVector, rCurrentProcessInfo );

    // //2.-Calculate MassMatrix:

    // MatrixType MassMatrix  = ZeroMatrix( system_size, system_size );

    // this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );
    
    // //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    // double alpha = 0;
    // if( GetProperties().Has(RAYLEIGH_ALPHA) )
    // { 
	// alpha = GetProperties()[RAYLEIGH_ALPHA];
    // }
    // else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
    // { 
	// alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    // }

    // double beta  = 0;
    // if( GetProperties().Has(RAYLEIGH_BETA) )
    // {
	// beta = GetProperties()[RAYLEIGH_BETA];
    // }
    // else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
    // { 
	// beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    // }

    // //4.-Compose the Damping Matrix:
   
    // //Rayleigh Damping Matrix: alpha*M + beta*K
    // MassMatrix      *= alpha;
    // StiffnessMatrix *= beta;

    // rDampingMatrix  = MassMatrix;
    // rDampingMatrix += StiffnessMatrix;

    std::cout << "CalculateDampingMatrix end" << std::endl;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

// double SpringDamperElement3D2N::CalculateCurrentLength(){

// 		KRATOS_TRY
// 		const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
// 			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
// 		const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
// 			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
// 		const double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
// 			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
// 		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
// 		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
// 		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
// 		const double l = sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy) 
// 			+ (dw + dz)*(dw + dz));
// 		return l;
// 		KRATOS_CATCH("")
// }

int SpringDamperElement3D2N::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    std::cout << "Check" << std::endl;

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
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
	    }

        if ( this->GetGeometry()[i].SolutionStepsDataHas( ROTATION ) == false )
	    {
            KRATOS_ERROR << "missing variable ROTATION on node " << this->GetGeometry()[i].Id() << std::endl;
	    }

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
	    {
            KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
	    }
    }

    return 0;

    KRATOS_CATCH( "" );
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


