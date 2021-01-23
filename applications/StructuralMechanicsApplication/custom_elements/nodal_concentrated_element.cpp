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
#include "includes/checks.h"
#include "custom_elements/nodal_concentrated_element.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"

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
    return Kratos::make_intrusive<NodalConcentratedElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mUseRayleighDamping );
}

//************************************************************************************
//************************************************************************************

Element::Pointer NodalConcentratedElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_intrusive<NodalConcentratedElement>( NewId, pGeom, pProperties, mUseRayleighDamping );
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

    return Kratos::make_intrusive<NodalConcentratedElement>(new_element);
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
    const ProcessInfo& rCurrentProcessInfo
    ) const
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
    const ProcessInfo& rCurrentProcessInfo
    ) const
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

void NodalConcentratedElement::GetValuesVector( Vector& rValues, int Step ) const
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

void NodalConcentratedElement::GetFirstDerivativesVector( Vector& rValues, int Step ) const
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

void NodalConcentratedElement::GetSecondDerivativesVector( Vector& rValues, int Step ) const
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


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

void NodalConcentratedElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

    // We get the reference
    const auto& rconst_this = *this;

    // Compute and add external forces
    const double nodal_mass = rconst_this.GetValue(NODAL_MASS);

    for ( unsigned int j = 0; j < dimension; ++j )
        rRightHandSideVector[j]  += volume_acceleration[j] * nodal_mass;

    // Compute and add internal forces
    const array_1d<double, 3 >& nodal_stiffness = rconst_this.GetValue(NODAL_DISPLACEMENT_STIFFNESS);
    for ( unsigned int j = 0; j < dimension; ++j )
        rRightHandSideVector[j]  -= nodal_stiffness[j] * current_displacement[j];
}

//***********************************************************************************
//***********************************************************************************

void NodalConcentratedElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int system_size = dimension;

    if ( rLeftHandSideMatrix.size1() != system_size )
        rLeftHandSideMatrix.resize( system_size, system_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

    // We get the reference
    const auto& rconst_this = *this;

    // We add the nodal stiffness
    const array_1d<double, 3 >& nodal_stiffness = rconst_this.GetValue(NODAL_DISPLACEMENT_STIFFNESS);
    for ( unsigned int j = 0; j < dimension; ++j )
        rLeftHandSideMatrix(j, j) += nodal_stiffness[j];
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

void NodalConcentratedElement::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int system_size = dimension;

    if ( rMassMatrix.size1() != system_size )
        rMassMatrix.resize( system_size, system_size, false );

    rMassMatrix = ZeroMatrix( system_size, system_size );

    // We get the reference
    const auto& rconst_this = *this;

    // Get the nodal mass
    const double nodal_mass = rconst_this.GetValue(NODAL_MASS);

    for ( unsigned int j = 0; j < dimension; ++j )
        rMassMatrix( j, j ) = nodal_mass;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
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
        StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
            *this,
            rDampingMatrix,
            rCurrentProcessInfo,
            system_size);
    } else {
        const auto& rconst_this = *this;
        const array_1d<double, 3 >& nodal_damping_ratio = rconst_this.GetValue(NODAL_DAMPING_RATIO);
        for ( unsigned int j = 0; j < dimension; ++j ) {
            rDampingMatrix(j, j) += nodal_damping_ratio[j];
        }
    }

    KRATOS_CATCH( "" );
}


void NodalConcentratedElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& rconst_this = *this;

    if (rDestinationVariable == NODAL_MASS) {
        double& r_nodal_mass = GetGeometry()[0].GetValue(NODAL_MASS);
        #pragma omp atomic
        r_nodal_mass += rconst_this.GetValue(NODAL_MASS);
    }

    KRATOS_CATCH("")
}

void NodalConcentratedElement::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType local_size = dimension*number_of_nodes;
    const auto& rconst_this = *this;

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {

        Vector damping_residual_contribution = ZeroVector(local_size);
        Vector current_nodal_velocities = ZeroVector(local_size);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix;
        ProcessInfo temp_process_information; // cant pass const ProcessInfo
        CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);

        for (SizeType i = 0; i < number_of_nodes; ++i) {
            SizeType index = dimension * i;
            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < dimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    } else if (rDestinationVariable == NODAL_INERTIA) {
        double& r_nodal_mass = GetGeometry()[0].GetValue(NODAL_MASS);
        #pragma omp atomic
        r_nodal_mass += rconst_this.GetValue(NODAL_MASS);
        // no contribution for GetGeometry()[0].GetValue(NODAL_INERTIA);
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

int NodalConcentratedElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_MASS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_DISPLACEMENT_STIFFNESS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_DAMPING_RATIO)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( std::size_t i = 0; i < this->GetGeometry().size(); ++i ) {
        const Node<3>& rnode = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rnode)
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
