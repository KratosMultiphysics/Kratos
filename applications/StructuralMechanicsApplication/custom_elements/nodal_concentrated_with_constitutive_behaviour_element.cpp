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
#include "custom_elements/nodal_concentrated_with_constitutive_behaviour_element.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
/***********************************************************************************/

NodalConcentratedWithConstitutiveBehaviourElement::NodalConcentratedWithConstitutiveBehaviourElement(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    ConstitutiveLaw::Pointer pConstitutiveLaw,
    const bool UseRayleighDamping,
    const bool ComputeActiveNodeFlag
    )
    : NodalConcentratedElement( NewId, pGeometry, UseRayleighDamping, ComputeActiveNodeFlag)
    , mpConstitutiveLaw(pConstitutiveLaw)
{

}

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

NodalConcentratedWithConstitutiveBehaviourElement::NodalConcentratedWithConstitutiveBehaviourElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    const bool UseRayleighDamping,
    const bool ComputeActiveNodeFlag
    )
    : NodalConcentratedElement( NewId, pGeometry, UseRayleighDamping, ComputeActiveNodeFlag)
{
}

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

NodalConcentratedWithConstitutiveBehaviourElement::NodalConcentratedWithConstitutiveBehaviourElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties, 
    const bool UseRayleighDamping,
    const bool ComputeActiveNodeFlag
    )
    : NodalConcentratedElement( NewId, pGeometry, pProperties, UseRayleighDamping, ComputeActiveNodeFlag)
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

NodalConcentratedWithConstitutiveBehaviourElement::NodalConcentratedWithConstitutiveBehaviourElement( NodalConcentratedWithConstitutiveBehaviourElement const& rOther)
    :BaseType(rOther)
    ,mpConstitutiveLaw(rOther.mpConstitutiveLaw)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
/***********************************************************************************/

NodalConcentratedWithConstitutiveBehaviourElement&  NodalConcentratedWithConstitutiveBehaviourElement::operator=(NodalConcentratedWithConstitutiveBehaviourElement const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    BaseType::operator=(rOther);
    mpConstitutiveLaw = rOther.mpConstitutiveLaw;

    return *this;
}

//*********************************OPERATIONS*****************************************
/***********************************************************************************/

Element::Pointer NodalConcentratedWithConstitutiveBehaviourElement::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_shared<NodalConcentratedWithConstitutiveBehaviourElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_RAYLEIGH_DAMPING) );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer NodalConcentratedWithConstitutiveBehaviourElement::Create(
    IndexType NewId, 
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties 
    ) const
{
    //NEEDED TO CREATE AN ELEMENT   
    return Kratos::make_shared<NodalConcentratedWithConstitutiveBehaviourElement>( NewId, pGeom, pProperties, mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_RAYLEIGH_DAMPING) );
}

//************************************CLONE*******************************************
/***********************************************************************************/

Element::Pointer NodalConcentratedWithConstitutiveBehaviourElement::Clone(
    IndexType NewId, 
    NodesArrayType const& rThisNodes 
    ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    return Kratos::make_shared<NodalConcentratedWithConstitutiveBehaviourElement>(NewId, GetGeometry().Create( rThisNodes ), pGetProperties(), mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_RAYLEIGH_DAMPING) );
}


//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

NodalConcentratedWithConstitutiveBehaviourElement::~NodalConcentratedWithConstitutiveBehaviourElement()
{
}
//*********************************DISPLACEMENT***************************************
/***********************************************************************************/

//************* STARTING - ENDING  METHODS
/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::Initialize()
{
    KRATOS_TRY;

    // We get the constitutive law from properties if not directly defined
    if (mpConstitutiveLaw == nullptr) {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(CONSTITUTIVE_LAW)) << "NodalConcentratedWithConstitutiveBehaviourElement requires the definition of a constitutive law in the properties" << std::endl;
        if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
                mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        } else {
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
        }
    }

    // Defining auxiliar zero array
    const array_1d<double, 3> zero_array(3, 0.0);

    // We check the nodal stiffness
    if (mpConstitutiveLaw->Has(NODAL_STIFFNESS)) {
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_DISPLACEMENT_STIFFNESS, true);
        this->SetValue(INITIAL_DISPLACEMENT, zero_array);
    } else
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_DISPLACEMENT_STIFFNESS, false);

    // We check the nodal rotational stiffness
    if (GetGeometry()[0].SolutionStepsDataHas(ROTATION_X) &&
        mpConstitutiveLaw->Has(NODAL_ROTATIONAL_STIFFNESS)) {
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_ROTATIONAL_STIFFNESS, true);
        this->SetValue(INITIAL_ROTATION, zero_array);
    } else
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_ROTATIONAL_STIFFNESS, false);

    // We check the nodal mass
    if (mpConstitutiveLaw->Has(NODAL_MASS))
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_MASS, true);
    else
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_MASS, false);

    // We check the nodal inertia
    if (GetGeometry()[0].SolutionStepsDataHas(ROTATION_X) &&
        mpConstitutiveLaw->Has(NODAL_INERTIA))
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_INERTIA, true);
    else
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_INERTIA, false);

    // We check the nodal damping
    if (mpConstitutiveLaw->Has(NODAL_DAMPING_RATIO))
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_DAMPING_RATIO, true);
    else
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_DAMPING_RATIO, false);

    // We check the nodal rtotational damping
    if (GetGeometry()[0].SolutionStepsDataHas(ROTATION_X) &&
        mpConstitutiveLaw->Has(NODAL_ROTATIONAL_DAMPING_RATIO))
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_ROTATIONAL_DAMPING_RATIO, true);
    else
        mELementalFlags.Set(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_ROTATIONAL_DAMPING_RATIO, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

//************* COMPUTING  METHODS
/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::CalculateLocalSystem(
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

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // The domain size
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters constitutive_law_parameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Resizing as needed the RHS
    const SizeType system_size = ComputeSizeOfSystem();

    if ( rRightHandSideVector.size() != system_size )
        rRightHandSideVector.resize( system_size, false );

    rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

    // Volume acceleration
    if (mELementalFlags.Is(NodalConcentratedElement::COMPUTE_NODAL_MASS) && GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION)) {

        const array_1d<double, 3 >& volume_acceleration = GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        // Compute and add external forces
        double nodal_mass = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_MASS, nodal_mass);
        for ( IndexType j = 0; j < dimension; ++j )
            rRightHandSideVector[j]  += volume_acceleration[j] * nodal_mass;
    }

    // Auxiliar index
    IndexType aux_index = 0;

    // The displacement terms
    if( mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_DISPLACEMENT_STIFFNESS) ||
        mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_MASS)) {

        // Compute and add internal forces
        const array_1d<double, 3 >& current_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 >& initial_displacement = this->GetValue(INITIAL_DISPLACEMENT);
        array_1d<double, 3> nodal_stiffness = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_STIFFNESS, nodal_stiffness) ;

        for ( IndexType j = 0; j < dimension; ++j )
            rRightHandSideVector[j]  -= nodal_stiffness[j] * (current_displacement[j] - initial_displacement[j]);

        aux_index += dimension;
    }

    // The rotational terms
    if( mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_ROTATIONAL_STIFFNESS) ||
        mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_INERTIA)) {

        // Compute and add internal forces
        const array_1d<double, 3 >& current_rotation = GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
        const array_1d<double, 3 >& initial_rotation = this->GetValue(INITIAL_ROTATION);
        array_1d<double, 3 > nodal_rotational_stiffness = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_ROTATIONAL_STIFFNESS, nodal_rotational_stiffness);

        for ( IndexType j = 0; j < dimension; ++j )
            rRightHandSideVector[aux_index + j]  -= nodal_rotational_stiffness[j] * (current_rotation[j] - initial_rotation[j]);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // The domain size
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters constitutive_law_parameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Resizing as needed the LHS
    const SizeType system_size = ComputeSizeOfSystem();

    if ( rLeftHandSideMatrix.size1() != system_size )
        rLeftHandSideMatrix.resize( system_size, system_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS
    
    // Auxiliar index
    IndexType aux_index = 0;

    // The displacement terms
    if( mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_DISPLACEMENT_STIFFNESS) ||
        mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_MASS)) {

        array_1d<double, 3 > nodal_stiffness = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_STIFFNESS, nodal_stiffness);

        for ( IndexType j = 0; j < dimension; ++j )
            rLeftHandSideMatrix( j, j ) += nodal_stiffness[j];

        aux_index += dimension;
    }

    // The rotational terms
    if( mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_ROTATIONAL_STIFFNESS) ||
        mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_INERTIA)) {

        array_1d<double, 3 > nodal_rotational_stiffness = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_ROTATIONAL_STIFFNESS, nodal_rotational_stiffness);

        for ( IndexType j = 0; j < dimension; ++j )
            rLeftHandSideMatrix( aux_index + j, aux_index + j ) += nodal_rotational_stiffness[j];
    }
}

//*************************COMPUTE DELTA POSITION*************************************
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters constitutive_law_parameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Lumped (by definition)
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = ComputeSizeOfSystem();

    if ( rMassMatrix.size1() != system_size )
        rMassMatrix.resize( system_size, system_size, false );

    rMassMatrix = ZeroMatrix( system_size, system_size );
    
    // Auxiliar index
    IndexType aux_index = 0;

    // The displacement terms
    if( mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_DISPLACEMENT_STIFFNESS) ||
        mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_MASS)) {

        double nodal_mass = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_MASS, nodal_mass);

        for ( IndexType j = 0; j < dimension; ++j )
            rMassMatrix( j, j ) = nodal_mass;

        aux_index += dimension;
    }

    // The rotational terms
    if( mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_ROTATIONAL_STIFFNESS) ||
        mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_NODAL_INERTIA)) {

        array_1d<double, 3> nodal_inertia = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_INERTIA, nodal_inertia);

        for ( IndexType j = 0; j < dimension; ++j )
            rMassMatrix( aux_index + j, aux_index + j ) = nodal_inertia[j];
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, 
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // The domain size
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType system_size = ComputeSizeOfSystem();

    // 0.-Initialize the DampingMatrix:
    rDampingMatrix = ZeroMatrix( system_size, system_size );

    //Check, if Rayleigh damping is available; use nodal damping, if not
    if( mELementalFlags.Is(NodalConcentratedWithConstitutiveBehaviourElement::COMPUTE_RAYLEIGH_DAMPING) ) {
        //1.-Calculate StiffnessMatrix:

        MatrixType stiffness_matrix       = ZeroMatrix( system_size, system_size );
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
        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters constitutive_law_parameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        IndexType aux_index = 0;
        if( mELementalFlags.Is(NodalConcentratedElement::COMPUTE_DAMPING_RATIO) ) {
            array_1d<double, 3> nodal_damping_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_DAMPING_RATIO, nodal_damping_ratio);
            for ( IndexType j = 0; j < dimension; ++j )
                rDampingMatrix(j, j) += nodal_damping_ratio[j];

            aux_index += dimension;
        }
        if( mELementalFlags.Is(NodalConcentratedElement::COMPUTE_ROTATIONAL_DAMPING_RATIO) ) {
            array_1d<double, 3> nodal_rotational_damping_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_parameters, NODAL_ROTATIONAL_DAMPING_RATIO, nodal_rotational_damping_ratio);
            for ( IndexType j = 0; j < dimension; ++j )
                rDampingMatrix(aux_index + j, aux_index+ j) += nodal_rotational_damping_ratio[j];
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void NodalConcentratedWithConstitutiveBehaviourElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    rSerializer.save("ConstitutiveLaw",mpConstitutiveLaw);
}

void NodalConcentratedWithConstitutiveBehaviourElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    rSerializer.load("ConstitutiveLaw",mpConstitutiveLaw);
}

} // Namespace Kratos


