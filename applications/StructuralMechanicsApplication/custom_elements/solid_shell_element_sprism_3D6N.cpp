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
#include "includes/checks.h"
#include "input_output/logger.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_elements/solid_shell_element_sprism_3D6N.h"

namespace Kratos
{
/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, EAS_IMPLICIT_EXPLICIT,              4 ); // True means implicit // TODO: change this using templates!!!
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, TOTAL_UPDATED_LAGRANGIAN,           5 ); // True means total lagrangian // TODO: change this using templates!!!
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, QUADRATIC_ELEMENT,                  6 ); // True means quadratic in-plane behaviour // TODO: Idem

// ------------------------------------------------------------------------- //
// ------------------------------ PUBLIC ----------------------------------- //
// ------------------------------------------------------------------------- //

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SolidShellElementSprism3D6N::SolidShellElementSprism3D6N( )
    : BaseType( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SolidShellElementSprism3D6N::SolidShellElementSprism3D6N(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SolidShellElementSprism3D6N::SolidShellElementSprism3D6N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

/*********************************** COPY CONSTRUCTOR ******************************/
/***********************************************************************************/

SolidShellElementSprism3D6N::SolidShellElementSprism3D6N( SolidShellElementSprism3D6N const& rOther)
    :BaseType(rOther)
    ,mFinalizedStep(rOther.mFinalizedStep)
    ,mAuxContainer(rOther.mAuxContainer)
{
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

SolidShellElementSprism3D6N::~SolidShellElementSprism3D6N()
{
}

/********************************** ASSIGMENT OPERATOR *****************************/
/***********************************************************************************/

SolidShellElementSprism3D6N&  SolidShellElementSprism3D6N::operator=(SolidShellElementSprism3D6N const& rOther)
{
    BaseType::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mAuxContainer.clear();
    mAuxContainer.resize( rOther.mAuxContainer.size());

    for(IndexType i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
        mAuxContainer[i]=rOther.mAuxContainer[i];
    }

    return *this;
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Element::Pointer SolidShellElementSprism3D6N::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<SolidShellElementSprism3D6N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//************************************************************************************
//************************************************************************************

Element::Pointer SolidShellElementSprism3D6N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<SolidShellElementSprism3D6N>(NewId, pGeom, pProperties);
}

/*********************************** CLONE ******************************************/
/************************************************************************************/

Element::Pointer SolidShellElementSprism3D6N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    SolidShellElementSprism3D6N new_element( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    new_element.mThisIntegrationMethod = this->GetIntegrationMethod();

    const SizeType integration_point_number = mConstitutiveLawVector.size();

    if ( new_element.mConstitutiveLawVector.size() != integration_point_number)
        new_element.mConstitutiveLawVector.resize(integration_point_number);

    KRATOS_ERROR_IF( new_element.mConstitutiveLawVector.size() != new_element.GetGeometry().IntegrationPointsNumber() ) << "Constitutive law not has the correct size " << new_element.mConstitutiveLawVector.size() << std::endl;

    for(IndexType i = 0; i < integration_point_number; ++i)
        new_element.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();

    if ( new_element.mAuxContainer.size() != mAuxContainer.size() )
        new_element.mAuxContainer.resize(mAuxContainer.size());

    for(IndexType i = 0; i < mAuxContainer.size(); ++i)
        new_element.mAuxContainer[i] = mAuxContainer[i];

    return Kratos::make_shared<SolidShellElementSprism3D6N>(new_element);
}

//******************************* GETTING METHODS *********************************//
/***********************************************************************************/

void SolidShellElementSprism3D6N::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    const IndexType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);
    const IndexType dim = number_of_nodes * 3;

    if (rResult.size() != dim)
        rResult.resize(dim, false);

    // Nodes of the central element
    IndexType index = 0;
    for (IndexType i = 0; i < 6; ++i) {
        rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        index += 3;
    }

    // Adding the ids of the neighbouring nodes
    for (IndexType i = 0; i < 6; ++i) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            rResult[index]     = p_neighbour_nodes[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = p_neighbour_nodes[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = p_neighbour_nodes[i].GetDof(DISPLACEMENT_Z).EquationId();
            index += 3;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    rElementalDofList.resize(0);

    // Nodes of the central element
    for (IndexType i = 0; i < GetGeometry().size(); ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }

    // Adding the dofs of the neighbouring nodes
    for (IndexType i = 0; i < 6; ++i) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            rElementalDofList.push_back(p_neighbour_nodes[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(p_neighbour_nodes[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(p_neighbour_nodes[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("");
}

/******************************** DISPLACEMENT **************************************/
/************************************************************************************/

void SolidShellElementSprism3D6N::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    const SizeType mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    IndexType index = 0;

    // Nodes of the central element
    for (IndexType i = 0; i < 6; ++i) {
        const array_1d<double, 3 > & disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        for (IndexType j = 0; j < 3; ++j)
            rValues[index + j] = disp[j];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; ++i) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            const array_1d<double, 3 > & disp = p_neighbour_nodes[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            for (IndexType j = 0; j < 3; ++j)
                rValues[index + j] = disp[j];
            index += 3;
        }
    }
}

/********************************** VELOCITY ****************************************/
/************************************************************************************/

void SolidShellElementSprism3D6N::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    const SizeType mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    IndexType index = 0;

    // Nodes of the central element
    for (IndexType i = 0; i < 6; ++i) {
        const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        for (IndexType j = 0; j < 3; ++j)
            rValues[index + j] = vel[j];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; ++i) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            const array_1d<double, 3 > & vel = p_neighbour_nodes[i].FastGetSolutionStepValue(VELOCITY, Step);
            for (IndexType j = 0; j < 3; ++j)
                rValues[index + j] = vel[j];
            index += 3;
        }
    }
}

/******************************** ACCELERATION **************************************/
/************************************************************************************/

void SolidShellElementSprism3D6N::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    const SizeType mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    IndexType index = 0;

    // Nodes of the central element
    for (IndexType i = 0; i < 6; ++i) {
        const array_1d<double, 3 > & acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        for (IndexType j = 0; j < 3; ++j)
            rValues[index + j] = acc[j];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; ++i) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            const array_1d<double, 3 > & acc = p_neighbour_nodes[i].FastGetSolutionStepValue(ACCELERATION, Step);
            for (IndexType j = 0; j < 3; ++j)
                rValues[index + j] = acc[j];
            index += 3;
        }
    }
}

//****************************** COMPUTING METHODS ********************************//
/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR);

    MatrixType left_hand_side_matrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( left_hand_side_matrix, rRightHandSideVector, local_system.CalculationFlags );

    //Set general_variables to Local system components
    local_system.SetLeftHandSideMatrix(left_hand_side_matrix);
    local_system.SetRightHandSideVector(rRightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateRightHandSide(
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR);
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType left_hand_side_matrix = Matrix();

    /* Initialize sizes for the system components: */
    if( rRHSVariables.size() != rRightHandSideVectors.size() ) {
        rRightHandSideVectors.resize(rRHSVariables.size());
    }

    for( IndexType i = 0; i < rRightHandSideVectors.size(); ++i ) {
        this->InitializeSystemMatrices( left_hand_side_matrix, rRightHandSideVectors[i], local_system.CalculationFlags );
    }

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrix(left_hand_side_matrix);
    local_system.SetRightHandSideVectors(rRightHandSideVectors);

    local_system.SetRightHandSideVariables(rRHSVariables);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX, true);
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR, false);

    VectorType right_hand_side_vector = Vector();

    /* Initialize sizes for the system components: */
    this->InitializeSystemMatrices( rLeftHandSideMatrix, right_hand_side_vector, local_system.CalculationFlags );

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    local_system.SetRightHandSideVector(right_hand_side_vector);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX, true);
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR, true);

    /* Initialize sizes for the system components: */
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, local_system.CalculationFlags );

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    local_system.SetRightHandSideVector(rRightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateLocalSystem(
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags*/
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    /* Initialize sizes for the system components: */
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() ) {
        rLeftHandSideMatrices.resize(rLHSVariables.size());
    }

    if( rRHSVariables.size() != rRightHandSideVectors.size() ) {
        rRightHandSideVectors.resize(rRHSVariables.size());
    }

    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX);
    for( IndexType i = 0; i < rLeftHandSideMatrices.size(); ++i ) {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], local_system.CalculationFlags );
    }

    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR, true);
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX, false);

    for( IndexType i = 0; i < rRightHandSideVectors.size(); ++i ) {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], local_system.CalculationFlags );
    }

    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX, true);

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    local_system.SetRightHandSideVectors(rRightHandSideVectors);

    local_system.SetLeftHandSideVariables(rLHSVariables);
    local_system.SetRightHandSideVariables(rRHSVariables);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);
    const SizeType mat_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != mat_size) {
        rMassMatrix.resize(mat_size, mat_size, false);
    }

    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

    Matrix aux_matrix;
    const SizeType aux_mat_size = GetGeometry().size() * 3;
    BaseType::CalculateMassMatrix(aux_matrix, rCurrentProcessInfo);
    noalias(subrange(rMassMatrix, 0, aux_mat_size, 0, aux_mat_size)) = aux_matrix;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const IndexType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    // Resizing as needed the LHS
    const IndexType mat_size = number_of_nodes * 3;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Calculate StiffnessMatrix:

    MatrixType stiffness_matrix  = Matrix();

    this->CalculateLeftHandSide( stiffness_matrix, rCurrentProcessInfo );

    // 2.-Calculate mass matrix:

    MatrixType mass_matrix  = Matrix();

    this->CalculateMassMatrix ( mass_matrix, rCurrentProcessInfo );

    // 3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) ) {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ) {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) ) {
        beta = GetProperties()[RAYLEIGH_BETA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ) {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    // 4.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * mass_matrix;
    noalias( rDampingMatrix ) += beta  * stiffness_matrix;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const MatrixType& rStiffnessMatrix,
    const MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * 3;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) ) {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ) {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) ) {
        beta = GetProperties()[RAYLEIGH_BETA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ) {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    // 2.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * rMassMatrix;
    noalias( rDampingMatrix ) += beta  * rStiffnessMatrix;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const IndexType integration_point_number = GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number )
        rOutput.resize( integration_point_number, false );

    if ( rVariable == VON_MISES_STRESS ) {
        /* Create and initialize element variables: */
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags &ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy (Values);

            const Matrix& stress_tensor = MathUtils<double>::StressVectorToTensor(general_variables.StressVector); //reduced dimension stress tensor


            // In general coordinates:
            double sigma_equivalent =  (0.5)*((stress_tensor(0,0)-stress_tensor(1,1))*((stress_tensor(0,0)-stress_tensor(1,1)))+
                                              (stress_tensor(1,1)-stress_tensor(2,2))*((stress_tensor(1,1)-stress_tensor(2,2)))+
                                              (stress_tensor(2,2)-stress_tensor(0,0))*((stress_tensor(2,2)-stress_tensor(0,0)))+
                                            6*(stress_tensor(0,1)*stress_tensor(1,0)+stress_tensor(1,2)*stress_tensor(2,1)+stress_tensor(2,0)*stress_tensor(0,2)));

            if( sigma_equivalent < 0 )
                sigma_equivalent = 0;

            sigma_equivalent = std::sqrt(sigma_equivalent);

            rOutput[point_number] =  sigma_equivalent;
        }
    } else if ( rVariable == NORM_ISOCHORIC_STRESS ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points,point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy (Values);

            const Matrix& stress_tensor  = MathUtils<double>::StressVectorToTensor(general_variables.StressVector); //reduced dimension stress tensor

            double stress_norm =  ((stress_tensor(0,0)*stress_tensor(0,0))+(stress_tensor(1,1)*stress_tensor(1,1))+(stress_tensor(2,2)*stress_tensor(2,2))+
                                   (stress_tensor(0,1)*stress_tensor(0,1))+(stress_tensor(0,2)*stress_tensor(0,2))+(stress_tensor(1,2)*stress_tensor(1,2))+
                                   (stress_tensor(1,0)*stress_tensor(1,0))+(stress_tensor(2,0)*stress_tensor(2,0))+(stress_tensor(2,1)*stress_tensor(2,1)));

            stress_norm = std::sqrt(stress_norm);

            rOutput[point_number] = stress_norm;
        }
    } else if ( rVariable == STRAIN_ENERGY ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double ZetaGauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, ZetaGauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            double strain_energy = 0.0;

            // Compute stresses and constitutive parameters
            if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN))
                mConstitutiveLawVector[point_number]->CalculateMaterialResponseKirchhoff(Values);
            else
                mConstitutiveLawVector[point_number]->CalculateMaterialResponsePK2(Values);

            mConstitutiveLawVector[point_number]->GetValue(STRAIN_ENERGY, strain_energy);

            rOutput[point_number] = general_variables.detJ * integration_points[point_number].Weight() * strain_energy;  // 1/2 * sigma * epsilon
        }
    } else if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType ii = 0; ii < integration_point_number; ++ii )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    } else {
        /* Create and initialize element variables: */
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags &ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            rOutput[point_number] = mConstitutiveLawVector[point_number]->CalculateValue( Values, rVariable, rOutput[point_number] );
        }
    }

    if ( rOutput.size() != 6 ) {
        std::vector<double> r_output_aux;
        r_output_aux = rOutput;

        rOutput.resize( 6, false );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (IndexType iii = 0; iii < 6; ++iii) {
            rOutput[iii] = 0.0;

            for (IndexType i_gp = 0; i_gp < integration_point_number; i_gp++)
                rOutput[iii] += interpol(i_gp, iii) * r_output_aux[i_gp];
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const IndexType integration_point_number = GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number )
        rOutput.resize( integration_point_number );

    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &constitutive_laws_options=Values.GetOptions();

        constitutive_laws_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        constitutive_laws_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables, Values, point_number);

            // Call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR)
                general_variables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
            else
                general_variables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

            mConstitutiveLawVector[point_number]->CalculateMaterialResponse(Values, general_variables.StressMeasure);

            if (rOutput[point_number].size() != general_variables.StressVector.size())
                rOutput[point_number].resize( general_variables.StressVector.size(), false);
            rOutput[point_number] = general_variables.StressVector;
        }
    } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR || rVariable == HENCKY_STRAIN_VECTOR) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags &ConstitutiveLawOptions=values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        values.SetStrainVector(general_variables.StrainVector);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables, values, point_number);

            // Compute Green-Lagrange Strain
            if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ) {
                mConstitutiveLawVector[point_number]->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
            } else if( rVariable == ALMANSI_STRAIN_VECTOR ) {
                mConstitutiveLawVector[point_number]->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_Cauchy);
            } else if( rVariable == HENCKY_STRAIN_VECTOR ) {
                mConstitutiveLawVector[point_number]->CalculateValue(values, HENCKY_STRAIN_VECTOR, general_variables.StrainVector);
            }

            if (rOutput[point_number].size() != general_variables.StrainVector.size())
                rOutput[point_number].resize( general_variables.StrainVector.size(), false );

            rOutput[point_number] = general_variables.StrainVector;
        }
    } else if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType ii = 0; ii < integration_point_number; ++ii )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    } else {
        /* Create and initialize element variables: */
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags &ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            rOutput[point_number] = mConstitutiveLawVector[point_number]->CalculateValue( Values, rVariable, rOutput[point_number] );
        }
    }

    if ( rOutput.size() != 6 ) {
        std::vector<Vector> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (IndexType iii = 0; iii < 6; ++iii) {
            rOutput[iii] = ZeroVector(rOutput[0].size());

            for (IndexType i_gp = 0; i_gp < integration_point_number; i_gp++)
                rOutput[iii] += interpol(i_gp, iii) * rOutput_aux[i_gp];
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateOnIntegrationPoints(
    const Variable<Matrix >& rVariable,
    std::vector< Matrix >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const IndexType integration_point_number = GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number )
        rOutput.resize( integration_point_number );

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
        std::vector<Vector> stress_vector;
        if( rVariable == CAUCHY_STRESS_TENSOR )
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
        else
            this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

        // Loop integration points
        if ( rOutput.size() != stress_vector.size() )
            rOutput.resize( stress_vector.size() );

        for ( IndexType point_number = 0; point_number < rOutput.size(); ++point_number ) {
            if (rOutput[point_number].size2() != 3)
                rOutput[point_number].resize(3, 3, false);
            rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR || rVariable == HENCKY_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        else if ( rVariable == ALMANSI_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        else if ( rVariable == HENCKY_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( HENCKY_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        // Loop integration points
        if ( rOutput.size() != StrainVector.size() )
            rOutput.resize( StrainVector.size() );

        for ( IndexType point_number = 0; point_number < rOutput.size(); ++point_number ) {
            if (rOutput[point_number].size2() != 3)
                rOutput[point_number].resize(3, 3, false);

            rOutput[point_number] = MathUtils<double>::StrainVectorToTensor(StrainVector[point_number]);
        }
    } else if ( rVariable == CONSTITUTIVE_MATRIX ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &constitutive_laws_options=Values.GetOptions();
        constitutive_laws_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy(Values);

            if( rOutput[point_number].size2() != general_variables.ConstitutiveMatrix.size2() ) {
                rOutput[point_number].resize( general_variables.ConstitutiveMatrix.size1() , general_variables.ConstitutiveMatrix.size2() , false );
            }
            rOutput[point_number] = general_variables.ConstitutiveMatrix;
        }
    } else if ( rVariable == DEFORMATION_GRADIENT ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            if( rOutput[point_number].size2() != general_variables.F.size2() )
                rOutput[point_number].resize( general_variables.F.size1() , general_variables.F.size2() , false );
            rOutput[point_number] = general_variables.F;
        }
    } else if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType ii = 0; ii < integration_point_number; ++ii )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    } else {
        /* Create and initialize element variables: */
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags &ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            rOutput[point_number] = mConstitutiveLawVector[point_number]->CalculateValue( Values, rVariable, rOutput[point_number] );
        }
    }

    if ( rOutput.size() != 6 ) {
        std::vector<Matrix> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (IndexType iii = 0; iii < 6; ++iii) {
            rOutput[iii] = ZeroMatrix(rOutput[0].size1(), rOutput[0].size2());

            for (IndexType Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
        }
    }

    KRATOS_CATCH( "" );
}

//**************************** ON INTEGRATION POINTS ******************************//
/******************************** SET DOUBLE VALUE *********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    BaseType::SetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/******************************** SET VECTOR VALUE *********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/******************************** SET MATRIX VALUE *********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/****************************** SET CONSTITUTIVE VALUE *****************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/******************************** GET DOUBLE VALUE *********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/********************************** GET VECTOR VALUE *******************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/*********************************** GET MATRIX VALUE ******************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/******************************** GET CONSTITUTIVE VALUE ***************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//********************************* CHECK VALUES **********************************//
/***********************************************************************************/
/***********************************************************************************/

int  SolidShellElementSprism3D6N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    /* Check the neighbours have been calculated */
    // Neighbour elements
    WeakPointerVector< Element >& p_neighbour_elements = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(p_neighbour_elements.size() == 0) << "The neighbour elements are not calculated" << std::endl;

    // Neighbour nodes
    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    KRATOS_ERROR_IF(p_neighbour_nodes.size() == 0) << "The neighbour nodes are not calculated" << std::endl;

    const int check = BaseType::Check(rCurrentProcessInfo);

    // Verify that nodal variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(VON_MISES_STRESS)
    KRATOS_CHECK_VARIABLE_KEY(NORM_ISOCHORIC_STRESS)
    KRATOS_CHECK_VARIABLE_KEY(CAUCHY_STRESS_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(CAUCHY_STRESS_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(PK2_STRESS_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(PK2_STRESS_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(GREEN_LAGRANGE_STRAIN_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(GREEN_LAGRANGE_STRAIN_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(ALMANSI_STRAIN_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(ALMANSI_STRAIN_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(HENCKY_STRAIN_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(HENCKY_STRAIN_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_MATRIX)
    KRATOS_CHECK_VARIABLE_KEY(DEFORMATION_GRADIENT)

    /* Verify compatibility with the constitutive law */
    ConstitutiveLaw::Features law_features;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(law_features);

    // Check strain measure
    bool correct_strain_measure = false;
    for(IndexType i = 0; i < law_features.mStrainMeasures.size(); ++i)
        if (law_features.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
            correct_strain_measure = true;

    KRATOS_ERROR_IF_NOT(correct_strain_measure) << "Constitutive law is not compatible with the element type SolidShellElementSprism3D6N" << std::endl;

    return check;

    KRATOS_CATCH( "" );
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    mFinalizedStep = false;
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Create and initialize element variables:
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Get constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    double& alpha_eas = this->GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives this_cartesian_derivatives;
    this->CalculateCartesianDerivatives(this_cartesian_derivatives);

    /* Calculate common components (B, C) */
    CommonComponents common_components;
    common_components.clear();
    this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

    // Reading integration points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components,integration_points, point_number, alpha_eas, zeta_gauss);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(general_variables,Values,point_number);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(Values, general_variables.StressMeasure);

        // Call the constitutive law to finalize the solution step
        mConstitutiveLawVector[point_number]->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( this->GetIntegrationMethod() ), point_number ),
                rCurrentProcessInfo );

        // Call the element internal variables update
        this->FinalizeStepVariables(general_variables, point_number);
    }

    mFinalizedStep = true;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add something if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    /* Create and initialize element variables: */
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags &ConstitutiveLawOptions=values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    /* Getting the alpha parameter of the EAS improvement */
    double& alpha_eas = this->GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives this_cartesian_derivatives;
    this->CalculateCartesianDerivatives(this_cartesian_derivatives);

    /* Calculate common components (B, C) */
    CommonComponents common_components;
    common_components.clear();
    this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

    /* Reset the EAS integrated components */
    EASComponents EAS;
    EAS.clear();

    // Reading integration points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(general_variables.B, common_components, zeta_gauss, alpha_eas);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(general_variables, values, point_number);

        // Compute stresses and constitutive parameters
        mConstitutiveLawVector[point_number]->CalculateMaterialResponse(values, general_variables.StressMeasure);

        // Calculating weights for integration on the "reference configuration"
        const double integration_weight = integration_points[point_number].Weight() * general_variables.detJ;

        /* Integrate in Zeta */
        // EAS components
        IntegrateEASInZeta(general_variables, EAS, zeta_gauss, integration_weight);
    }

    /* Getting the increase of displacements */
    BoundedMatrix<double, 36, 1 > delta_disp;

    delta_disp = GetVectorCurrentPosition() - GetVectorPreviousPosition(); // Calculates the increase of displacements

    /* Update alpha EAS */
    if (EAS.mStiffAlpha > std::numeric_limits<double>::epsilon()) // Avoid division by zero
        alpha_eas -= prod(EAS.mHEAS, delta_disp)(0, 0) / EAS.mStiffAlpha;
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::Initialize()
{
    KRATOS_TRY;

    mFinalizedStep = true; // the creation is out of the time step, it must be true

    if( GetProperties().Has(INTEGRATION_ORDER) ) {
        const SizeType integration_order = GetProperties()[INTEGRATION_ORDER];
        switch ( integration_order )
        {
        case 1:
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_1;
            break;
        case 2:
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_2;
            break;
        case 3:
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_3;
            break;
        case 4:
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_4;
            break;
        case 5:
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
            break;
        default:
            KRATOS_WARNING("SolidShellElementSprism3D6N") << "Integration order " << integration_order << " is not available, using default integration order for the geometry" << std::endl;
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_1;
        }
    } else {
        mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_1;
    }

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    /* Constitutive Law initialisation */
    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    /* Implicit or explicit EAS update */
    if( GetProperties().Has(CONSIDER_IMPLICIT_EAS_SPRISM_ELEMENT) )
        mELementalFlags.Set(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT, GetProperties()[CONSIDER_IMPLICIT_EAS_SPRISM_ELEMENT]);
    else
        mELementalFlags.Set(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT, true);

    /* Total or updated lagrangian */
    if( GetProperties().Has(CONSIDER_TOTAL_LAGRANGIAN_SPRISM_ELEMENT) )
        mELementalFlags.Set(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN, GetProperties()[CONSIDER_TOTAL_LAGRANGIAN_SPRISM_ELEMENT]);
    else
        mELementalFlags.Set(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN, true);

    /* Quadratic or linear element */
    if( GetProperties().Has(CONSIDER_QUADRATIC_SPRISM_ELEMENT) )
        mELementalFlags.Set(SolidShellElementSprism3D6N::QUADRATIC_ELEMENT, GetProperties()[CONSIDER_QUADRATIC_SPRISM_ELEMENT]);
    else
        mELementalFlags.Set(SolidShellElementSprism3D6N::QUADRATIC_ELEMENT, true);

    // Resizing the containers
    mAuxContainer.resize( integration_points.size() );

    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) { // Jacobian inverses
        // Compute jacobian inverses and set the domain initial size:
        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian(J0, this->GetIntegrationMethod());

        /* Calculating the inverse J0 */
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            // Calculating and storing inverse of the jacobian and the parameters needed
            double aux_detJ;
            MathUtils<double>::InvertMatrix( J0[point_number], mAuxContainer[point_number], aux_detJ );
        }
    } else { // Historic deformation gradient
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            mAuxContainer[point_number] = IdentityMatrix(3);
        }
    }

    /* Initialize AlphaEAS */
    this->SetValue(ALPHA_EAS, 0.0);

    /* Material initialisation */
    InitializeMaterial();

    KRATOS_CATCH("");
}

/***********************************PROTECTED***************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateElementalSystem(
    LocalSystemComponents& rLocalSystem,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create and initialize element variables: */
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags &ConstitutiveLawOptions=values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    /* Getting the alpha parameter of the EAS improvement */
    double& alpha_eas = this->GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives this_cartesian_derivatives;
    this->CalculateCartesianDerivatives(this_cartesian_derivatives);

    /* Calculate common components (B, C) */
    CommonComponents common_components;
    common_components.clear();
    this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

    /* Reset the integrated stress components */
    StressIntegratedComponents rIntegratedStress;
    rIntegratedStress.clear();

    /* Reset the EAS integrated components */
    EASComponents EAS;
    EAS.clear();

    /* Auxiliary terms: Allocating the VolumeForce*/
    Vector volume_force = ZeroVector(3);

    // Reading integration points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(general_variables.B, common_components, zeta_gauss, alpha_eas);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components, integration_points, point_number, alpha_eas, zeta_gauss);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(general_variables, values, point_number);

        // Compute stresses and constitutive parameters
        mConstitutiveLawVector[point_number]->CalculateMaterialResponse(values, general_variables.StressMeasure);

        // Calculating weights for integration on the "reference configuration"
        const double integration_weight = integration_points[point_number].Weight() * general_variables.detJ;

        /* Integrate in Zeta */
        // Stresses
        IntegrateStressesInZeta(general_variables, rIntegratedStress, alpha_eas, zeta_gauss, integration_weight);
        // EAS components
        IntegrateEASInZeta(general_variables, EAS, zeta_gauss, integration_weight);

        if ( rLocalSystem.CalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR) ) { // Calculation of the vector is required
            /* Volume forces */
            this->CalculateVolumeForce( volume_force, general_variables, integration_weight );
        }
    }

    /* Calculate the RHS */
    if ( rLocalSystem.CalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR) ) { // Calculation of the vector is required
        /* Contribution to external and internal forces */
        this->CalculateAndAddRHS ( rLocalSystem, general_variables, volume_force, rIntegratedStress, common_components, EAS, alpha_eas );
    }

    if ( rLocalSystem.CalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX) ) { // Calculation of the matrix is required
        /* Contribution to the tangent stiffness matrix */
        this->CalculateAndAddLHS( rLocalSystem, general_variables, values, rIntegratedStress, common_components, this_cartesian_derivatives, EAS, alpha_eas );
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::PrintElementCalculation(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables
    )
{
    KRATOS_TRY;

    KRATOS_INFO("SolidShellElementSprism3D6N") << " Element: " << this->Id() << std::endl;

    WeakPointerVectorNodesType& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    const IndexType number_of_neighbours = NumberOfActiveNeighbours(NeighbourNodes);

    for ( IndexType i = 0; i < 6; ++i ) {
        const array_1d<double, 3> &current_position  = GetGeometry()[i].Coordinates();
        const array_1d<double, 3 > & current_displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & previous_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        const array_1d<double, 3> previous_position  = current_position - (current_displacement-previous_displacement);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous  Position  node[" << GetGeometry()[i].Id() << "]: "<<previous_position << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3> &current_position  = NeighbourNodes[i].Coordinates();
        const array_1d<double, 3 > & current_displacement  = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & previous_displacement = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        const array_1d<double, 3> previous_position  = current_position - (current_displacement-previous_displacement);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous  Position  neighbour node[" << NeighbourNodes[i].Id() << "]: "<<previous_position << std::endl;
    }

    for ( IndexType i = 0; i < 6; ++i ) {
        const array_1d<double, 3> & current_position  = GetGeometry()[i].Coordinates();
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Position  node[" << GetGeometry()[i].Id()<<"]: " << current_position << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3> & current_position  = NeighbourNodes[i].Coordinates();
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Position neighbour node[" << NeighbourNodes[i].Id()<<"]: " << current_position << std::endl;
    }

    for ( IndexType i = 0; i < 6; ++i ) {
        const array_1d<double, 3 > & previous_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous Displacement node[" << GetGeometry()[i].Id() << "]: " << previous_displacement << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3 > & previous_displacement = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous Displacement neighbour node[" << NeighbourNodes[i].Id() << "]: " << previous_displacement << std::endl;
    }

    for ( IndexType i = 0; i < 6; ++i ) {
        const array_1d<double, 3 > & current_displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Displacement  node[" << GetGeometry()[i].Id() << "]: " << current_displacement << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3 > & current_displacement  = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Displacement  node[" << NeighbourNodes[i].Id() << "]: " << current_displacement << std::endl;
    }

    KRATOS_INFO("SolidShellElementSprism3D6N") << " Stress " << rVariables.StressVector << std::endl;
    KRATOS_INFO("SolidShellElementSprism3D6N") << " Strain " << rVariables.StrainVector << std::endl;
    KRATOS_INFO("SolidShellElementSprism3D6N") << " F  " << rVariables.F<<std::endl;
    KRATOS_INFO("SolidShellElementSprism3D6N") << " ConstitutiveMatrix " <<rVariables.ConstitutiveMatrix << std::endl;
    KRATOS_INFO("SolidShellElementSprism3D6N") << " K " << rLocalSystem.GetLeftHandSideMatrix() << std::endl;
    KRATOS_INFO("SolidShellElementSprism3D6N") << " f " << rLocalSystem.GetRightHandSideVector() << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

bool SolidShellElementSprism3D6N::HasNeighbour(
    const IndexType Index,
    const NodeType& NeighbourNode
    )
{
    if (NeighbourNode.Id() == GetGeometry()[Index].Id()) {
        return false;
    } else {
        if ( mELementalFlags.Is(SolidShellElementSprism3D6N::QUADRATIC_ELEMENT))
            return true;
        else
            return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SolidShellElementSprism3D6N::NumberOfActiveNeighbours(WeakPointerVectorNodesType& pNeighbourNodes)
{
    std::size_t active_neighbours = 0;
    for (IndexType i = 0; i < pNeighbourNodes.size(); ++i) {
        if (HasNeighbour(i, pNeighbourNodes[i]))
           ++active_neighbours;
    }
    return active_neighbours;
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetNodalCoordinates(
    BoundedMatrix<double, 12, 3 > & NodesCoord,
    WeakPointerVectorNodesType& NeighbourNodes,
    const Configuration ThisConfiguration
    )
{
     NodesCoord = ZeroMatrix(12, 3);
     const IndexType number_of_neighbours = NumberOfActiveNeighbours(NeighbourNodes);

     if (ThisConfiguration == Configuration::INITIAL) {
         /* Fill the aux matrix of coordinates */
         for (IndexType i = 0; i < 6; ++i) {
             const array_1d<double, 3> &initial_position = GetGeometry()[i].GetInitialPosition().Coordinates();
             for (IndexType j = 0; j < 3; ++j)
                 NodesCoord(i, j) = initial_position[j];
         }

         if (number_of_neighbours == 6) { // All the possible neighbours
             for (IndexType i = 0; i < 6; ++i) {
                 const array_1d<double, 3> &initial_position = NeighbourNodes[i].GetInitialPosition().Coordinates();
                 for (IndexType j = 0; j < 3; ++j)
                    NodesCoord(i + 6, j) = initial_position[j];
             }
         } else {
             for (IndexType i = 0; i < 6; ++i) {
                 if (HasNeighbour(i, NeighbourNodes[i])) {
                     const array_1d<double, 3> &initial_position = NeighbourNodes[i].GetInitialPosition().Coordinates();

                     for (IndexType j = 0; j < 3; ++j)
                        NodesCoord(i + 6, j) = initial_position[j];

                 } else {
                     for (IndexType j = 0; j < 3; ++j)
                        NodesCoord(i + 6, j) = 0.0;
                 }
             }
         }
     } else if (ThisConfiguration == Configuration::CURRENT) {
         /* Fill the aux matrix of coordinates */
         for (IndexType i = 0; i < 6; ++i) {
             const array_1d<double, 3>& current_position  = GetGeometry()[i].Coordinates();
             for (IndexType j = 0; j < 3; ++j)
                NodesCoord(i, j) = current_position[j];
         }

         if (number_of_neighbours == 6) { // All the possible neighours
             for (IndexType i = 0; i < 6; ++i) {
                 const array_1d<double, 3>& current_position  = NeighbourNodes[i].Coordinates();
                 for (IndexType j = 0; j < 3; ++j)
                    NodesCoord(i + 6, j) = current_position[j];
             }
         } else {
             for (IndexType i = 0; i < 6; ++i) {
                 if (HasNeighbour(i, NeighbourNodes[i])) {
                     const array_1d<double, 3>& current_position  = NeighbourNodes[i].Coordinates();
                     for (IndexType j = 0; j < 3; ++j)
                        NodesCoord(i + 6, j) = current_position[j];
                 } else {
                     for (IndexType j = 0; j < 3; ++j)
                        NodesCoord(i + 6, j) = 0.0;
                 }
             }
         }
     } else {
         const std::string& config = (ThisConfiguration == Configuration::INITIAL) ? "Initial" : "Current";
         KRATOS_ERROR << " The configuration is not possible, the posibilities are Current and Initial: " << config << std::endl;
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerivatives(CartesianDerivatives& rCartesianDerivatives)
{
    BoundedMatrix<double, 12, 3 > nodes_coord; // Coordinates of the nodes
    WeakPointerVectorNodesType& neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        this->GetNodalCoordinates(nodes_coord, neighbour_nodes, Configuration::INITIAL);
    } else {
        this->GetNodalCoordinates(nodes_coord, neighbour_nodes, Configuration::CURRENT);
    }

    /* Calculate local system of coordinates of the element */
    const double ang_rot = GetProperties().Has(ANG_ROT) ? GetProperties()[ANG_ROT] : 0.0; // TODO: Change to consider multiple plies
    OrthogonalBase this_orthogonal_base;
    this->CalculateLocalCoordinateSystem(this_orthogonal_base, OrthogonalBaseApproach::Z, ang_rot);

    //******************************** CENTRAL POINT ******************************
    // Calculate cartesian derivatives
    BoundedMatrix<double, 2, 4 > cartesian_derivatives_center_lower,  cartesian_derivatives_center_upper;

    // Lower face
    CalculateCartesianDerivativesOnCenterPlane(cartesian_derivatives_center_lower, this_orthogonal_base, GeometricLevel::LOWER);
    // Upperr face
    CalculateCartesianDerivativesOnCenterPlane(cartesian_derivatives_center_upper, this_orthogonal_base, GeometricLevel::UPPER );

    /* Transversal derivative */
    CalculateCartesianDerOnCenterTrans(rCartesianDerivatives, nodes_coord, this_orthogonal_base, GeometricLevel::CENTER); // Center
    CalculateCartesianDerOnCenterTrans(rCartesianDerivatives, nodes_coord, this_orthogonal_base, GeometricLevel::LOWER);  // Lower part
    CalculateCartesianDerOnCenterTrans(rCartesianDerivatives, nodes_coord, this_orthogonal_base, GeometricLevel::UPPER);  // Upper part

    //******************************** GAUSS POINTS *******************************

    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 0.5;
    local_coordinates[1] = 0.5;
    local_coordinates[2] = -1.0;

    /* Transversal derivatives */
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[0], nodes_coord, this_orthogonal_base, local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[3], nodes_coord, this_orthogonal_base, local_coordinates);
    local_coordinates[0] = 0.0;
    local_coordinates[2] = -1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[1], nodes_coord, this_orthogonal_base, local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[4], nodes_coord, this_orthogonal_base, local_coordinates);
    local_coordinates[0] = 0.5;
    local_coordinates[1] = 0.0;
    local_coordinates[2] = -1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[2], nodes_coord, this_orthogonal_base, local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[5], nodes_coord, this_orthogonal_base, local_coordinates);

    /* In-plane derivative */
    for (IndexType i = 0; i < 3 ;++i) {
        if (HasNeighbour(i, neighbour_nodes[i])) { // Assuming that if the upper element has neighbours the lower has too
            CalculateCartesianDerOnGaussPlane(rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i    ], nodes_coord, this_orthogonal_base, i, GeometricLevel::LOWER);
            CalculateCartesianDerOnGaussPlane(rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i + 3], nodes_coord, this_orthogonal_base, i, GeometricLevel::UPPER);
        } else {
            noalias(rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i    ]) = cartesian_derivatives_center_lower;
            noalias(rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i + 3]) = cartesian_derivatives_center_upper;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCommonComponents(
    CommonComponents& rCommonComponents,
    const CartesianDerivatives& rCartesianDerivatives
    )
{
    KRATOS_TRY;

    BoundedMatrix<double, 12, 3 > NodesCoord; // Coordinates of the nodes
    WeakPointerVectorNodesType& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    this->GetNodalCoordinates(NodesCoord, NeighbourNodes, Configuration::CURRENT);

    /* Declare deformation Gradient F components */
    // In plane components
    BoundedMatrix<double, 3, 2 > in_plane_gradient_F_gauss;
    // Transversal components
    TransverseGradient transverse_gradient;

    //*****************************************************************************

    /* COMPUTATION OF B TANGENTS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Calculating the membrane strain-displacement matrix */
    // Lower face
    for (IndexType i = 0; i < 3; ++i) {
        CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i], NodesCoord, i, GeometricLevel::LOWER);
        CalculateAndAddBMembrane(rCommonComponents.BMembraneLower, rCommonComponents.CMembraneLower, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i], in_plane_gradient_F_gauss, i);
    }

    rCommonComponents.BMembraneLower *= 1.0/3.0;
    rCommonComponents.CMembraneLower *= 1.0/3.0;

    // Upper face
    for (IndexType i = 0; i < 3; ++i) {
        CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i + 3], NodesCoord, i, GeometricLevel::UPPER);
        CalculateAndAddBMembrane(rCommonComponents.BMembraneUpper, rCommonComponents.CMembraneUpper, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i + 3], in_plane_gradient_F_gauss, i);
    }

    rCommonComponents.BMembraneUpper *= 1.0/3.0;
    rCommonComponents.CMembraneUpper *= 1.0/3.0;

    /* SHEAR CONTRIBUTION */
    /* Calculating the shear strain-displacement matrix */

    // Declaring the isoparametric transverse gradient variables
    TransverseGradientIsoParametric transverse_gradient_isoparametric;

    // Lower face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(transverse_gradient_isoparametric, NodesCoord, GeometricLevel::LOWER);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(transverse_gradient.F0, rCartesianDerivatives.TransversalCartesianDerivativesGauss[0], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F1, rCartesianDerivatives.TransversalCartesianDerivativesGauss[1], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F2, rCartesianDerivatives.TransversalCartesianDerivativesGauss[2], NodesCoord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAddBShear(rCommonComponents.BShearLower, rCommonComponents.CShearLower, rCartesianDerivatives, transverse_gradient, transverse_gradient_isoparametric, GeometricLevel::LOWER);

    // Upper face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(transverse_gradient_isoparametric, NodesCoord, GeometricLevel::UPPER);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(transverse_gradient.F0, rCartesianDerivatives.TransversalCartesianDerivativesGauss[3], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F1, rCartesianDerivatives.TransversalCartesianDerivativesGauss[4], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F2, rCartesianDerivatives.TransversalCartesianDerivativesGauss[5], NodesCoord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAddBShear(rCommonComponents.BShearUpper, rCommonComponents.CShearUpper, rCartesianDerivatives, transverse_gradient, transverse_gradient_isoparametric, GeometricLevel::UPPER);

    /* NORMAL TRANSVERSE */
    /* Calculate f normal components */
    array_1d<double, 3 > F3;
    CalculateTransverseGradientF(F3, rCartesianDerivatives.TransversalCartesianDerivativesCenter, NodesCoord);

    /* Calculating the normal transverse strain-displacement matrix */
    CalculateAndAddBNormal(rCommonComponents.BNormal, rCommonComponents.CNormal, rCartesianDerivatives.TransversalCartesianDerivativesCenter, F3);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateLocalCoordinateSystem(
    OrthogonalBase& ThisOrthogonalBase,
    const OrthogonalBaseApproach ThisOrthogonalBaseApproach,
    const double ThisAngle
    )
{
    KRATOS_TRY;

    /* Mid-surface vectors */
    double norm; // TODO: Use the geometry normal when avalaible
    array_1d<double, 3 > vxe, vye;
    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        vxe[0] = 0.5 * ((GetGeometry()[2].X0() + GetGeometry()[5].X0()) - (GetGeometry()[1].X0() + GetGeometry()[4].X0()));
        vxe[1] = 0.5 * ((GetGeometry()[2].Y0() + GetGeometry()[5].Y0()) - (GetGeometry()[1].Y0() + GetGeometry()[4].Y0()));
        vxe[2] = 0.5 * ((GetGeometry()[2].Z0() + GetGeometry()[5].Z0()) - (GetGeometry()[1].Z0() + GetGeometry()[4].Z0()));

        vye[0] = 0.5 * ((GetGeometry()[0].X0() + GetGeometry()[3].X0()) - (GetGeometry()[2].X0() + GetGeometry()[5].X0()));
        vye[1] = 0.5 * ((GetGeometry()[0].Y0() + GetGeometry()[3].Y0()) - (GetGeometry()[2].Y0() + GetGeometry()[5].Y0()));
        vye[2] = 0.5 * ((GetGeometry()[0].Z0() + GetGeometry()[3].Z0()) - (GetGeometry()[2].Z0() + GetGeometry()[5].Z0()));
    } else {
        vxe[0] = 0.5 * ((GetGeometry()[2].X() + GetGeometry()[5].X()) - (GetGeometry()[1].X() + GetGeometry()[4].X()));
        vxe[1] = 0.5 * ((GetGeometry()[2].Y() + GetGeometry()[5].Y()) - (GetGeometry()[1].Y() + GetGeometry()[4].Y()));
        vxe[2] = 0.5 * ((GetGeometry()[2].Z() + GetGeometry()[5].Z()) - (GetGeometry()[1].Z() + GetGeometry()[4].Z()));

        vye[0] = 0.5 * ((GetGeometry()[0].X() + GetGeometry()[3].X()) - (GetGeometry()[2].X() + GetGeometry()[5].X()));
        vye[1] = 0.5 * ((GetGeometry()[0].Y() + GetGeometry()[3].Y()) - (GetGeometry()[2].Y() + GetGeometry()[5].Y()));
        vye[2] = 0.5 * ((GetGeometry()[0].Z() + GetGeometry()[3].Z()) - (GetGeometry()[2].Z() + GetGeometry()[5].Z()));
    }

    MathUtils<double>::CrossProduct(ThisOrthogonalBase.Vzeta, vxe, vye);
    norm = norm_2(ThisOrthogonalBase.Vzeta);
    ThisOrthogonalBase.Vzeta /= norm;

    const double threshold = std::numeric_limits<double>::epsilon();
    double ortho_comp;

    /* Performing the calculation */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    if (ThisOrthogonalBaseApproach == OrthogonalBaseApproach::X) {
        ortho_comp = ThisOrthogonalBase.Vzeta[1] * ThisOrthogonalBase.Vzeta[1] + ThisOrthogonalBase.Vzeta[2] * ThisOrthogonalBase.Vzeta[2]; // Component in th Y-Z plane
        if (ortho_comp < threshold) { // If ThisOrthogonalBase.Vzeta is almost orthogonal to  Y-Z plane
            ThisOrthogonalBase.Veta[0] = - ThisOrthogonalBase.Vzeta[2]; // Choose ThisOrthogonalBase.Vxi orthogonal to global Y direction
            ThisOrthogonalBase.Veta[1] = 0.0;
            ThisOrthogonalBase.Veta[2] = ThisOrthogonalBase.Vzeta[0];

            norm = norm_2(ThisOrthogonalBase.Vxi);
            ThisOrthogonalBase.Vxi /= norm;
            MathUtils<double>::CrossProduct(ThisOrthogonalBase.Vxi, ThisOrthogonalBase.Veta, ThisOrthogonalBase.Vzeta);
        } else { // SELECT local y=ThisOrthogonalBase.Vxi in the global YZ plane
            ThisOrthogonalBase.Vxi[0] = 0.0;
            ThisOrthogonalBase.Vxi[1] = ThisOrthogonalBase.Vzeta[2];
            ThisOrthogonalBase.Vxi[2] = - ThisOrthogonalBase.Vzeta[1];

            norm = norm_2(ThisOrthogonalBase.Vxi);
            ThisOrthogonalBase.Vxi /= norm;

            ThisOrthogonalBase.Veta[0] = ortho_comp; // Choose ThisOrthogonalBase.Vxi orthogonal to global X direction
            ThisOrthogonalBase.Veta[1] = - ThisOrthogonalBase.Vzeta[0] * ThisOrthogonalBase.Vzeta[1];
            ThisOrthogonalBase.Veta[2] = - ThisOrthogonalBase.Vzeta[0] * ThisOrthogonalBase.Vzeta[2];

            norm = norm_2(ThisOrthogonalBase.Veta);
            ThisOrthogonalBase.Veta /= norm;
        }
    } else if (ThisOrthogonalBaseApproach == OrthogonalBaseApproach::Y) {
        ortho_comp = ThisOrthogonalBase.Vzeta[0] * ThisOrthogonalBase.Vzeta[0] + ThisOrthogonalBase.Vzeta[2] * ThisOrthogonalBase.Vzeta[2]; // Component in th Z-X plane
        if (ortho_comp < threshold) { // If vze is almost orthogonal to  Z-X plane
            ThisOrthogonalBase.Veta[0] =       0.0; // Choose ThisOrthogonalBase.Vxi orthogonal to global X direction
            ThisOrthogonalBase.Veta[1] =   ThisOrthogonalBase.Vzeta[2];
            ThisOrthogonalBase.Veta[2] = - ThisOrthogonalBase.Vzeta[1];

            norm = norm_2(ThisOrthogonalBase.Veta);
            ThisOrthogonalBase.Veta /= norm;
            MathUtils<double>::CrossProduct(ThisOrthogonalBase.Vxi, ThisOrthogonalBase.Veta, ThisOrthogonalBase.Vzeta);
        } else { // SELECT local z=ThisOrthogonalBase.Vxi in the global ZX plane
            ThisOrthogonalBase.Vxi[0] = - ThisOrthogonalBase.Vzeta[2]; // Choose ThisOrthogonalBase.Vxi orthogonal to global Y direction
            ThisOrthogonalBase.Vxi[1] = 0.0;
            ThisOrthogonalBase.Vxi[2] = - ThisOrthogonalBase.Vzeta[0];

            norm = norm_2(ThisOrthogonalBase.Vxi);
            ThisOrthogonalBase.Vxi /= norm;

            ThisOrthogonalBase.Veta[0] = - ThisOrthogonalBase.Vzeta[0] * ThisOrthogonalBase.Vzeta[1];
            ThisOrthogonalBase.Veta[1] = ortho_comp;
            ThisOrthogonalBase.Veta[2] = - ThisOrthogonalBase.Vzeta[2] * ThisOrthogonalBase.Vzeta[1];

            norm = norm_2(ThisOrthogonalBase.Veta);
            ThisOrthogonalBase.Veta /= norm;
        }
    } else if (ThisOrthogonalBaseApproach == OrthogonalBaseApproach::Z) {
        ortho_comp = ThisOrthogonalBase.Vzeta[0] * ThisOrthogonalBase.Vzeta[0] + ThisOrthogonalBase.Vzeta[1] * ThisOrthogonalBase.Vzeta[1]; // Component in th X-Y plane
        if (ortho_comp < threshold) { // If vze is almost orthogonal to  X-Y plane
            ThisOrthogonalBase.Veta[0] = 0.0; // Choose ThisOrthogonalBase.Vxi orthogonal to global X direction
            ThisOrthogonalBase.Veta[1] = ThisOrthogonalBase.Vzeta[2];
            ThisOrthogonalBase.Veta[2] = - ThisOrthogonalBase.Vzeta[1];

            norm = norm_2(ThisOrthogonalBase.Veta);
            ThisOrthogonalBase.Veta /= norm;
            MathUtils<double>::CrossProduct(ThisOrthogonalBase.Vxi, ThisOrthogonalBase.Veta, ThisOrthogonalBase.Vzeta);
        } else { // SELECT local x=ThisOrthogonalBase.Vxi in the global XY plane
            ThisOrthogonalBase.Vxi[0] = - ThisOrthogonalBase.Vzeta[1];
            ThisOrthogonalBase.Vxi[1] = ThisOrthogonalBase.Vzeta[0];
            ThisOrthogonalBase.Vxi[2] = 0.0;

            norm = norm_2(ThisOrthogonalBase.Vxi);
            ThisOrthogonalBase.Vxi /= norm;

            ThisOrthogonalBase.Veta[0] = - ThisOrthogonalBase.Vzeta[0] * ThisOrthogonalBase.Vzeta[2]; // Choose ThisOrthogonalBase.Vxi orthogonal to global Z direction
            ThisOrthogonalBase.Veta[1] = - ThisOrthogonalBase.Vzeta[1] * ThisOrthogonalBase.Vzeta[2];
            ThisOrthogonalBase.Veta[2] = ortho_comp;

            norm = norm_2(ThisOrthogonalBase.Veta);
            ThisOrthogonalBase.Veta /= norm;
        }
    } else {
        ThisOrthogonalBase.Vxi[0] = 1.0;
        ThisOrthogonalBase.Vxi[1] = 0.0;
        ThisOrthogonalBase.Vxi[2] = 0.0;

        ThisOrthogonalBase.Veta[0] = 0.0;
        ThisOrthogonalBase.Veta[1] = 1.0;
        ThisOrthogonalBase.Veta[2] = 0.0;
    }

    if (ThisAngle != 0.0) {
        // Compute angle between local system ThisOrthogonalBase.Vxi-ThisOrthogonalBase.Veta and L1
        const double cosa = std::cos(ThisAngle);
        const double sina = std::sin(ThisAngle);
        // Rotate local system ThisOrthogonalBase.Vxi-ThisOrthogonalBase.Veta to best fit L1-L2
        ThisOrthogonalBase.Vzeta = ThisOrthogonalBase.Vxi; // Reusing as auxiliar value
        ThisOrthogonalBase.Vxi =   cosa * ThisOrthogonalBase.Vxi + sina * ThisOrthogonalBase.Veta;
        ThisOrthogonalBase.Veta = - sina * ThisOrthogonalBase.Vzeta  + cosa * ThisOrthogonalBase.Veta;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateIdVector(array_1d<IndexType, 18 >& rIdVector)
{
    KRATOS_TRY;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    /* Compute ID vector */ // TODO: Optimize this
    IndexType index = 18;
    for (IndexType i = 0; i < 6; ++i) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            for (IndexType j = 0; j < 3; ++j) {
                rIdVector[i * 3 + j] = index;
                ++index;
            }
        } else {
            for (IndexType j = 0; j < 3; ++j) {
                rIdVector[i * 3 + j] = 36;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::ComputeLocalDerivatives(
    BoundedMatrix<double, 6, 3 > & LocalDerivativePatch,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    const double L_1 = 0.5 * (1.0 - rLocalCoordinates[2]);
    const double L_2 = 0.5 * (1.0 + rLocalCoordinates[2]);
//     const double zzeta = 1.0 - rLocalCoordinates[0] - rLocalCoordinates[1];

    /* Derivative in direction nu and xi */
    // Lower face
    LocalDerivativePatch(0, 0) = - L_1;
    LocalDerivativePatch(1, 0) =   L_1;
    LocalDerivativePatch(2, 0) =   0.0;

    LocalDerivativePatch(0, 1) = - L_1;
    LocalDerivativePatch(1, 1) =   0.0;
    LocalDerivativePatch(2, 1) =   L_1;

    // Upper face
    LocalDerivativePatch(3, 0) = - L_2;
    LocalDerivativePatch(4, 0) =   L_2;
    LocalDerivativePatch(5, 0) =   0.0;

    LocalDerivativePatch(3, 1) = - L_2;
    LocalDerivativePatch(4, 1) =   0.0;
    LocalDerivativePatch(5, 1) =   L_2;

//     /* Derivative in direction zeta */
//     LocalDerivativePatch(0, 2) = - zzeta/2.0;
//     LocalDerivativePatch(1, 2) = - rLocalCoordinates[0]/2.0;
//     LocalDerivativePatch(2, 2) = - rLocalCoordinates[1]/2.0;
//     LocalDerivativePatch(3, 2) =   zzeta/2.0;
//     LocalDerivativePatch(4, 2) =   rLocalCoordinates[0]/2.0;
//     LocalDerivativePatch(5, 2) =   rLocalCoordinates[1]/2.0;

    /* Derivative in direction zeta */
    LocalDerivativePatch(0, 2) = - 1.0 + rLocalCoordinates[1] + rLocalCoordinates[0];
    LocalDerivativePatch(1, 2) = - rLocalCoordinates[0];
    LocalDerivativePatch(2, 2) = - rLocalCoordinates[1];
    LocalDerivativePatch(3, 2) =   1.0 - rLocalCoordinates[1] - rLocalCoordinates[0];
    LocalDerivativePatch(4, 2) =   rLocalCoordinates[0];
    LocalDerivativePatch(5, 2) =   rLocalCoordinates[1];
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::ComputeLocalDerivativesQuadratic(
    BoundedMatrix<double, 4, 2 >& rLocalDerivativePatch,
    const IndexType NodeGauss
    )
{
    /* Local coordinates */
    double xi  = 0.0;
    double eta = 0.0;

    if (NodeGauss == 0) {
        xi  = 0.5;
        eta = 0.5;
    } else if (NodeGauss == 1) {
        xi  = 0.0;
        eta = 0.5;
    } else if (NodeGauss == 2) {
        xi  = 0.5;
        eta = 0.0;
    }

    /* Derivative in main nodes */
    rLocalDerivativePatch(0, 0) = - 1.0 + eta;
    rLocalDerivativePatch(0, 1) = - 1.0 + xi;
    rLocalDerivativePatch(1, 0) =   1.0 - eta;
    rLocalDerivativePatch(1, 1) =   1.0 - xi - 2.0 * eta;
    rLocalDerivativePatch(2, 0) =   1.0 - 2.0 * xi - eta;
    rLocalDerivativePatch(2, 1) =   1.0 - xi;

    /* Derivative in neighbour nodes */
    if (NodeGauss == 0) {
        rLocalDerivativePatch(3, 0) = xi + eta - 0.5;
        rLocalDerivativePatch(3, 1) = xi + eta - 0.5;
    } else if (NodeGauss == 1) {
        rLocalDerivativePatch(3, 0) = xi - 0.5;
        rLocalDerivativePatch(3, 1) = 0.0;
    } else if (NodeGauss == 2) {
        rLocalDerivativePatch(3, 0) = 0.0;
        rLocalDerivativePatch(3, 1) = eta - 0.5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateJacobianCenterGauss(
    GeometryType::JacobiansType& J,
    std::vector< Matrix >& Jinv,
    Vector& detJ,
    const IndexType rPointNumber,
    const double ZetaGauss
    )
{
    /* Fill the aux matrix of coordinates */
    BoundedMatrix<double, 3, 6 > nodes_coord;
    for (IndexType i = 0; i < 6; ++i) {
        const array_1d<double, 3> &current_position  = GetGeometry()[i].Coordinates();
        for (IndexType j = 0; j < 3; ++j)
            nodes_coord(j, i) = current_position[j];
    }

    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 1.0/3.0;
    local_coordinates[1] = 1.0/3.0;
    local_coordinates[2] = ZetaGauss;

    /* Local derivatives patch */
    BoundedMatrix<double, 6, 3 > LocalDerivativePatch;
    ComputeLocalDerivatives(LocalDerivativePatch, local_coordinates);

    /* Compute Jacobian */
    noalias(J[rPointNumber]) = prod(nodes_coord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    MathUtils<double>::InvertMatrix( J[rPointNumber], Jinv[rPointNumber], detJ[rPointNumber] );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateJacobian(
    double & detJ,
    BoundedMatrix<double, 3, 3 > & J,
    BoundedMatrix<double, 6, 3 > & LocalDerivativePatch,
    const BoundedMatrix<double, 12, 3 > & NodesCoord,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Auxiliar coordinates of the nodes */
    BoundedMatrix<double, 3, 6 > nodes_coord_aux;

    for (IndexType i = 0; i < 6; ++i)
        for (IndexType j = 0; j < 3; ++j)
            nodes_coord_aux(j, i) = NodesCoord(i, j);

    /* Local derivatives patch */
    ComputeLocalDerivatives(LocalDerivativePatch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(J) = prod(nodes_coord_aux, LocalDerivativePatch);

    /* Compute determinant */
    detJ = MathUtils<double>::Det3(J);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateJacobianAndInv(
    BoundedMatrix<double, 3, 3 >& J,
    BoundedMatrix<double, 3, 3 >& Jinv,
    BoundedMatrix<double, 6, 3 >& LocalDerivativePatch,
    const BoundedMatrix<double, 3, 6 >& NodesCoord,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Local derivatives patch */
    ComputeLocalDerivatives(LocalDerivativePatch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(J) = prod(NodesCoord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    Jinv = MathUtils<double>::InvertMatrix<3>(J, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateJacobianAndInv(
    BoundedMatrix<double, 3, 3 > & J,
    BoundedMatrix<double, 3, 3 > & Jinv,
    const BoundedMatrix<double, 3, 6 > & NodesCoord,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Local derivatives patch */
    BoundedMatrix<double, 6, 3 > local_derivative_patch;
    ComputeLocalDerivatives(local_derivative_patch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(J) = prod(NodesCoord, local_derivative_patch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    Jinv = MathUtils<double>::InvertMatrix<3>(J, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerivativesOnCenterPlane(
    BoundedMatrix<double, 2, 4 >& CartesianDerivativesCenter,
    const OrthogonalBase& ThisOrthogonalBase,
    const GeometricLevel Part
    )
{
    const IndexType index = Part == GeometricLevel::UPPER ? 3 : 0;

    double norm0, norm;
    array_1d<double, 3 > vxe, vye;
    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        vxe[0] = GetGeometry()[2 + index].X0() - GetGeometry()[1 + index].X0();
        vxe[1] = GetGeometry()[2 + index].Y0() - GetGeometry()[1 + index].Y0();
        vxe[2] = GetGeometry()[2 + index].Z0() - GetGeometry()[1 + index].Z0();

        vye[0] = GetGeometry()[0 + index].X0() - GetGeometry()[2 + index].X0();
        vye[1] = GetGeometry()[0 + index].Y0() - GetGeometry()[2 + index].Y0();
        vye[2] = GetGeometry()[0 + index].Z0() - GetGeometry()[2 + index].Z0();
    } else {
        vxe[0] = GetGeometry()[2 + index].X() - GetGeometry()[1 + index].X();
        vxe[1] = GetGeometry()[2 + index].Y() - GetGeometry()[1 + index].Y();
        vxe[2] = GetGeometry()[2 + index].Z() - GetGeometry()[1 + index].Z();

        vye[0] = GetGeometry()[0 + index].X() - GetGeometry()[2 + index].X();
        vye[1] = GetGeometry()[0 + index].Y() - GetGeometry()[2 + index].Y();
        vye[2] = GetGeometry()[0 + index].Z() - GetGeometry()[2 + index].Z();
    }

    array_1d<double, 3 > t1g, t2g, t3g;
    MathUtils<double>::CrossProduct(t3g, vxe, vye);
    norm0 = norm_2(t3g);
    t3g /= norm0;

    MathUtils<double>::CrossProduct(t2g, t3g, ThisOrthogonalBase.Vxi);
    norm = norm_2(t2g);
    t2g /= norm;

    MathUtils<double>::CrossProduct(t1g, t2g, t3g);
    norm = norm_2(t1g);
    t1g /= norm;

    array_1d<double, 3 > a, b;

    a[0] = inner_prod(vxe, t1g)/norm0;
    a[1] = inner_prod(vye, t1g)/norm0;
    a[2] = -(a[0] + a[1]);
    b[0] = inner_prod(vxe, t2g)/norm0;
    b[1] = inner_prod(vye, t2g)/norm0;
    b[2] = -(b[0] + b[1]);

    CartesianDerivativesCenter = ZeroMatrix(2, 4);
    for (IndexType i = 0; i < 3; ++i) {
       CartesianDerivativesCenter(0, i) = - b[i];
       CartesianDerivativesCenter(1, i) =   a[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerOnGaussPlane(
    BoundedMatrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3 > & NodesCoord,
    const OrthogonalBase& ThisOrthogonalBase,
    const IndexType NodeGauss,
    const GeometricLevel Part
    )
{
    const IndexType index = Part == GeometricLevel::UPPER ? 3 : 0;

    /* Local derivatives patch */
    BoundedMatrix<double, 4, 2 > local_derivative_patch;
    ComputeLocalDerivativesQuadratic(local_derivative_patch,NodeGauss);

    /* Auxiliar coordinates of the nodes */
    BoundedMatrix<double, 3, 4 > nodes_coord_aux;

    for (IndexType i = 0; i < 3; ++i) {
        for (IndexType j = 0; j < 3; ++j) {
            nodes_coord_aux(j, i) = NodesCoord(i + index, j);
        }
    }

    for (IndexType j = 0; j < 3; ++j)
        nodes_coord_aux(j, 3) = NodesCoord(NodeGauss + 6 + index, j);

    /* Compute local derivatives */
    const BoundedMatrix<double, 3, 2 > Xd = prod(nodes_coord_aux, local_derivative_patch);

    /* Split local derivatives */
    array_1d<double, 3 > Xdxi, Xdeta;
    Xdxi[0]  = Xd(0, 0);
    Xdxi[1]  = Xd(1, 0);
    Xdxi[2]  = Xd(2, 0);
    Xdeta[0] = Xd(0, 1);
    Xdeta[1] = Xd(1, 1);
    Xdeta[2] = Xd(2, 1);

    /* Compute orthonormal vectors */
    array_1d<double, 3 > t1g, t2g, t3g;
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t1g, t2g, t3g, ThisOrthogonalBase.Vxi, Xdxi, Xdeta);

    /* Compute Jacobian */
    BoundedMatrix<double, 2, 2 > jac;
    jac(0, 0) = inner_prod(Xdxi,  t1g);
    jac(0, 1) = inner_prod(Xdxi,  t2g);
    jac(1, 0) = inner_prod(Xdeta, t1g);
    jac(1, 1) = inner_prod(Xdeta, t2g);

    /* Compute the inverse of the Jacobian */
    double aux_det;
    const BoundedMatrix<double, 2, 2 > JinvPlane = MathUtils<double>::InvertMatrix<2>(jac, aux_det);

    /* Compute the Cartesian derivatives */
    noalias(InPlaneCartesianDerivativesGauss) = prod(JinvPlane, trans(local_derivative_patch));
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerOnGaussTrans(
    BoundedMatrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3 > & NodesCoord,
    const OrthogonalBase& ThisOrthogonalBase,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Compute local derivatives */
    double det;
    BoundedMatrix<double, 3, 3 > Xd;
    BoundedMatrix<double, 6, 3 > local_derivatives_patch;
    CalculateJacobian(det, Xd, local_derivatives_patch, NodesCoord, rLocalCoordinates);

    /* Split local derivatives */
    array_1d<double, 3 > Xdxi, Xdeta;
    Xdxi[0]  = Xd(0, 0);
    Xdxi[1]  = Xd(1, 0);
    Xdxi[2]  = Xd(2, 0);
    Xdeta[0] = Xd(0, 1);
    Xdeta[1] = Xd(1, 1);
    Xdeta[2] = Xd(2, 1);

    /* Compute orthonormal vectors */
    BoundedMatrix<double, 3, 3 > t;
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t, ThisOrthogonalBase.Vxi, Xdxi, Xdeta);

    /* Compute Jacobian */
    const BoundedMatrix<double, 3, 3 > jac = prod(t, Xd);

    /* Compute inverse of the Jaccobian (just third column) */
    BoundedMatrix<double, 3 ,1> JinvTrans;
    JinvTrans(0, 0) =   (jac(0, 1) * jac(1, 2) - jac(0, 2) * jac(1, 1)) / det;
    JinvTrans(1, 0) = - (jac(0, 0) * jac(1, 2) - jac(0, 2) * jac(1, 0)) / det;
    JinvTrans(2, 0) =   (jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1)) / det;

    /* Compute Cartesian derivatives */
    noalias(TransversalCartesianDerivativesGauss) = prod(local_derivatives_patch, JinvTrans);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerOnCenterTrans(
    CartesianDerivatives& rCartesianDerivatives,
    const BoundedMatrix<double, 12, 3 >& NodesCoord,
    const OrthogonalBase& ThisOrthogonalBase,
    const GeometricLevel Part
    )
{
    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 1.0/3.0;
    local_coordinates[1] = 1.0/3.0;

    if (Part == GeometricLevel::CENTER)
        local_coordinates[2] =   0.0;
    else if (Part == GeometricLevel::LOWER)
        local_coordinates[2] = - 1.0;
    else if (Part == GeometricLevel::UPPER)
        local_coordinates[2] =   1.0;

    /* Auxiliar coordinates of the nodes */
    BoundedMatrix<double, 3, 6 > nodes_coord_aux;
    for (IndexType i = 0; i < 6; ++i)
        for (IndexType j = 0; j < 3; ++j)
            nodes_coord_aux(j, i) = NodesCoord(i, j);

    /* Auxiliar components to calculate the Jacobian and his inverse */
    BoundedMatrix<double, 3, 3 > J, Jinv;

    if (Part == GeometricLevel::CENTER) {
        /* Calculate the Jacobian and his inverse */
        BoundedMatrix<double, 6, 3 > local_derivatives_patch;
        CalculateJacobianAndInv(J, Jinv, local_derivatives_patch, nodes_coord_aux, local_coordinates);

        // Compute cartesian (y3) derivatives of the shape functions necessary to compute f_3
        /* Compute Cartesian derivatives */
        const BoundedMatrix<double, 6, 3 > transverse_cartesian_derivatives_gauss_aux = prod(local_derivatives_patch, Jinv);

        for (IndexType i = 0; i < 6 ; ++i)
            rCartesianDerivatives.TransversalCartesianDerivativesCenter(i, 0) = inner_prod(ThisOrthogonalBase.Vzeta, row(transverse_cartesian_derivatives_gauss_aux, i));
     } else {
        /* Calculate the Jacobian and his inverse */
        CalculateJacobianAndInv(J, Jinv, nodes_coord_aux, local_coordinates);

         /* Split local derivatives */
         array_1d<double, 3 > Xdxi, Xdeta;
         Xdxi[0]   = Jinv(0, 0);
         Xdxi[1]   = Jinv(0, 1);
         Xdxi[2]   = Jinv(0, 2);
         Xdeta[0]  = Jinv(1, 0);
         Xdeta[1]  = Jinv(1, 1);
         Xdeta[2]  = Jinv(1, 2);

         /* Compute inverse of the Jaccobian (just in plane components)*/
         if (Part == GeometricLevel::LOWER) {
             rCartesianDerivatives.JInvPlaneLower(0, 0) = inner_prod(Xdxi,  ThisOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneLower(0, 1) = inner_prod(Xdeta, ThisOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneLower(1, 0) = inner_prod(Xdxi,  ThisOrthogonalBase.Veta);
             rCartesianDerivatives.JInvPlaneLower(1, 1) = inner_prod(Xdeta, ThisOrthogonalBase.Veta);
         } else if (Part == GeometricLevel::UPPER) {
             rCartesianDerivatives.JInvPlaneUpper(0, 0) = inner_prod(Xdxi,  ThisOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneUpper(0, 1) = inner_prod(Xdeta, ThisOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneUpper(1, 0) = inner_prod(Xdxi,  ThisOrthogonalBase.Veta);
             rCartesianDerivatives.JInvPlaneUpper(1, 1) = inner_prod(Xdeta, ThisOrthogonalBase.Veta);
         }
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateInPlaneGradientFGauss(
    BoundedMatrix<double, 3, 2>& InPlaneGradientFGauss,
    const BoundedMatrix<double, 2, 4>& InPlaneCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3>& NodesCoord,
    const IndexType NodeGauss,
    const GeometricLevel Part
    )
{
    /* Auxiliar operators */
    const IndexType index = Part == GeometricLevel::UPPER ? 3 : 0;
    BoundedMatrix<double, 3, 3 > nodes_coord_aux;
    BoundedMatrix<double, 3, 2 > in_plane_cartesian_derivatives_gauss_aux;

    for (IndexType i = 0; i < 3; ++i) {
        for (IndexType j = 0; j < 3; ++j)
            nodes_coord_aux(j, i) = NodesCoord(i + index, j);

        in_plane_cartesian_derivatives_gauss_aux(i, 0) = InPlaneCartesianDerivativesGauss(0, i);
        in_plane_cartesian_derivatives_gauss_aux(i, 1) = InPlaneCartesianDerivativesGauss(1, i);
    }

    noalias(InPlaneGradientFGauss) = prod(nodes_coord_aux, in_plane_cartesian_derivatives_gauss_aux);

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    if (HasNeighbour(NodeGauss, p_neighbour_nodes[NodeGauss])) {
        for (IndexType j = 0; j < 3 ; ++j) {
            InPlaneGradientFGauss(j, 0) += NodesCoord(NodeGauss + 6 + index, j) * InPlaneCartesianDerivativesGauss(0, 3);
            InPlaneGradientFGauss(j, 1) += NodesCoord(NodeGauss + 6 + index, j) * InPlaneCartesianDerivativesGauss(1, 3);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateTransverseGradientF(
    array_1d<double, 3 >& TransverseGradientF,
    const BoundedMatrix<double, 6, 1 >& TransversalCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3 >& NodesCoord
    )
{
    for (IndexType j = 0; j < 3; ++j) {
        TransverseGradientF[j] = 0.0;
    }

    for (IndexType i = 0; i < 6; ++i) {
        for (IndexType j = 0; j < 3; ++j) {
            TransverseGradientF[j] += TransversalCartesianDerivativesGauss(i, 0) * NodesCoord(i, j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateTransverseGradientFinP(
    TransverseGradientIsoParametric& TransverseGradientIsoParametric,
    const BoundedMatrix<double, 12, 3 > & NodesCoord,
    const GeometricLevel Part
    )
{
    const IndexType index = Part == GeometricLevel::UPPER ? 3 : 0;

    for (IndexType i = 0; i < 3; ++i) {
        TransverseGradientIsoParametric.Ft[i]   = NodesCoord(2 + index, i) - NodesCoord(1 + index, i);
        TransverseGradientIsoParametric.Fxi[i]  = NodesCoord(0 + index, i) - NodesCoord(2 + index, i);
        TransverseGradientIsoParametric.Feta[i] = NodesCoord(1 + index, i) - NodesCoord(0 + index, i);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddBMembrane(
    BoundedMatrix<double, 3, 18 >& BMembrane,
    BoundedMatrix<double, 3, 1  >& CMembrane,
    const BoundedMatrix<double, 2, 4 >& InPlaneCartesianDerivativesGauss,
    const BoundedMatrix<double, 3, 2 >& InPlaneGradientFGauss,
    const IndexType NodeGauss
    )
{
    for (IndexType i = 0; i < 4; ++i) {
        IndexType base = i * 3;
        if (i == 3) {
            base += NodeGauss * 3;
        }
        for (IndexType j = 0; j < 3; ++j) {
            BMembrane(0, base + j) += InPlaneCartesianDerivativesGauss(0, i) * InPlaneGradientFGauss(j, 0);
            BMembrane(1, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 1);
            BMembrane(2, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 0)
                                   +  InPlaneCartesianDerivativesGauss(0, i) * InPlaneGradientFGauss(j, 1);
        }
    }

    /* Calculate de componets of Cauchy tensor */
    // In plane auxiliar components
    array_1d<double, 3 > aux_deformation_gradient_F1, aux_deformation_gradient_F2;

    for (IndexType i = 0; i < 3; ++i) {
        aux_deformation_gradient_F1[i] = InPlaneGradientFGauss(i, 0);
        aux_deformation_gradient_F2[i] = InPlaneGradientFGauss(i, 1);
    }

    CMembrane(0, 0) += inner_prod(aux_deformation_gradient_F1, aux_deformation_gradient_F1);
    CMembrane(1, 0) += inner_prod(aux_deformation_gradient_F2, aux_deformation_gradient_F2);
    CMembrane(2, 0) += inner_prod(aux_deformation_gradient_F1, aux_deformation_gradient_F2);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddMembraneKgeometric(
    BoundedMatrix<double, 36, 36 > & Kgeometricmembrane,
    const CartesianDerivatives& rCartesianDerivatives,
    const array_1d<double, 3 > & SMembrane,
    const GeometricLevel Part
    )
{
    const IndexType index = static_cast<IndexType>(Part);
    const IndexType auxiliar_index = Part == GeometricLevel::UPPER ? 3 : 0;

    BoundedMatrix<double, 6, 6 > H = ZeroMatrix(6, 6);

    IndexType ii, jj;
    for (IndexType i = 0; i < 4; ++i) {
        for (IndexType j = 0; j < 4; ++j) {
            // Gauss 1
            ii = i;
            jj = j;
            H(ii, jj) += SMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, j)
                       + SMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, j)
                       + SMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, j)
                                       + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, j));

            // Gauss 2
            ii = (i ==  3) ? 4 : i;
            jj = (j ==  3) ? 4 : j;

            H(ii, jj) += SMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, j)
                       + SMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, j)
                       + SMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, j)
                                       + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, j));

            // Gauss 3
            ii = (i ==  3) ? 5 : i;
            jj = (j ==  3) ? 5 : j;

            H(ii, jj) += SMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, j)
                       + SMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, j)
                       + SMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, j)
                                       + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, j));
        }
    }

    H *= 1.0/3.0;

    // Assembling in Kgeometricmembrane
    IndexType rowindex, colindex;
    for (IndexType i = 0; i < 6; ++i) {
        rowindex = (i < 3) ? i * 3 + index : i * 3 + index + 9;
        for (IndexType j = i; j < 6; ++j) {
            colindex = (j < 3) ? j * 3 + index : j * 3 + index + 9;
            for(IndexType ii = 0; ii < 3; ++ii) {
                Kgeometricmembrane(rowindex + ii,colindex + ii) += H (i, j);
                if (rowindex != colindex) { // Skip diagonal
                    Kgeometricmembrane(colindex + ii, rowindex + ii) += H (i, j); // Symmetric part
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddBShear(
    BoundedMatrix<double, 2, 18 >& BShear,
    BoundedMatrix<double, 2, 1 >& CShear,
    const CartesianDerivatives& rCartesianDerivatives,
    const TransverseGradient& rTransverseGradient,
    const TransverseGradientIsoParametric& rTransverseGradientIsoParametric,
    const GeometricLevel Part
    )
{
    const IndexType index = static_cast<IndexType>(Part);
    const IndexType auxiliar_index = Part == GeometricLevel::UPPER ? 3 : 0;

    const BoundedMatrix<double, 2, 2 >& JInvPlane = Part == GeometricLevel::UPPER ? rCartesianDerivatives.JInvPlaneUpper : rCartesianDerivatives.JInvPlaneLower;

    // Considering the Gauss point in the middle of the element
    const double eta_p = 1.0/3.0;
    const double xi_p  = 1.0/3.0;
    BoundedMatrix<double, 2, 3 > Pa;
    Pa(0, 0) = - xi_p;
    Pa(0, 1) = - xi_p;
    Pa(0, 2) = 1.0 - xi_p;
    Pa(1, 0) = eta_p;
    Pa(1, 1) = eta_p - 1.0;
    Pa(1, 2) = eta_p;

    BoundedMatrix<double, 3, 18 > aux_b_shear = ZeroMatrix(3, 18);

    /* First contribution*/
    for (IndexType i = 0; i < 6; ++i) {
        IndexType base = i * 3;
        for (IndexType j = 0; j < 3; ++j) {
            aux_b_shear(0, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliar_index + 0](i, 0) * rTransverseGradientIsoParametric.Ft[j];
            aux_b_shear(1, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliar_index + 1](i, 0) * rTransverseGradientIsoParametric.Fxi[j];
            aux_b_shear(2, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliar_index + 2](i, 0) * rTransverseGradientIsoParametric.Feta[j];
        }
    }

    /* Second contibution */
    for (IndexType i = 0; i < 3; ++i) {
        /* First row */
        aux_b_shear(0, i + index + 3) -= rTransverseGradient.F0[i];
        aux_b_shear(0, i + index + 6) += rTransverseGradient.F0[i];

        /* Second row */
        aux_b_shear(1, i + index)     += rTransverseGradient.F1[i];
        aux_b_shear(1, i + index + 6) -= rTransverseGradient.F1[i];

        /* Third row */
        aux_b_shear(2, i + index)     -= rTransverseGradient.F2[i];
        aux_b_shear(2, i + index + 3) += rTransverseGradient.F2[i];
    }

    const BoundedMatrix<double, 2, 3 > aux_prod = prod(JInvPlane, Pa);
    noalias(BShear) = prod(aux_prod, aux_b_shear);

    // Calculating the components of C
    BoundedMatrix<double, 3, 1 > aux_c_shear;
    aux_c_shear(0, 0) = inner_prod(rTransverseGradientIsoParametric.Ft  , rTransverseGradient.F0);
    aux_c_shear(1, 0) = inner_prod(rTransverseGradientIsoParametric.Fxi , rTransverseGradient.F1);
    aux_c_shear(2, 0) = inner_prod(rTransverseGradientIsoParametric.Feta, rTransverseGradient.F2);

    noalias(CShear) = prod(aux_prod, aux_c_shear);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddShearKgeometric(
    BoundedMatrix<double, 36, 36 > & Kgeometricshear,
    const CartesianDerivatives& rCartesianDerivatives,
    const array_1d<double, 2 > & SShear,
    const GeometricLevel Part
    )
{
//     const IndexType index = static_cast<IndexType>(Part);
    const IndexType auxiliar_index = Part == GeometricLevel::UPPER ? 3 : 0;

    const BoundedMatrix<double, 2, 2 >& JInvPlane = Part == GeometricLevel::UPPER ? rCartesianDerivatives.JInvPlaneUpper : rCartesianDerivatives.JInvPlaneLower;

    const double Q1 = 1.0/3.0 * (SShear[0] * JInvPlane(0, 0) + SShear[1] * JInvPlane(0, 1));
    const double Q2 = 1.0/3.0 * (SShear[0] * JInvPlane(1, 0) + SShear[1] * JInvPlane(1, 1));

//    array_1d<double, 3 > q;
//    q[0] = -Q1 + Q2;
//    q[1] = -(Q1 + 2.0 * Q2);
//    q[2] = (2.0 * Q1 + Q2);

//    int delta;
//    if (index == 9)
//        delta = 3;
//    else
//        delta = 0;

//    for (IndexType i = 0; i < 3; ++i) { // For each DOF
//        /* First assembling */
//        Kgeometricshear(i + index + 3, i + index + 3) -= q[0] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliar_index](1 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) += q[0] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliar_index](2 + delta, 0);

//        /* Second assembling */
//        Kgeometricshear(i + index, i + index)         += q[1] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliar_index](0 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) -= q[1] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliar_index](2 + delta, 0);
//        /* Third assembling */
//        Kgeometricshear(i + index, i + index)         -= q[2] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliar_index](0 + delta, 0);
//        Kgeometricshear(i + index + 3, i + index + 3) += q[2] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliar_index](1 + delta, 0);
//    }

    array_1d<double, 3 > n1; // Side node with + contribution (previous DOF position)
    array_1d<double, 3 > n2; // Side node with - contribution (previous DOF position)

//    if (Part == GeometricLevel::LOWER) {
//        n1[0] = 6;
//        n1[1] = 3;
//        n1[2] = 0;

//        n2[0] = 6;
//        n2[1] = 3;
//        n2[2] = 0;
//    } else {
//        n1[0] = 15;
//        n1[1] = 12;
//        n1[2] = 9;

//        n2[0] = 15;
//        n2[1] = 12;
//        n2[2] = 9;
//    }

    // Note: Technically this is the correct one
    if (Part == GeometricLevel::LOWER) {
        n1[0] = 6;
        n1[1] = 0;
        n1[2] = 3;

        n2[0] = 3;
        n2[1] = 6;
        n2[2] = 0;
    } else {
        n1[0] = 15;
        n1[1] = 9;
        n1[2] = 12;

        n2[0] = 12;
        n2[1] = 15;
        n2[2] = 9;
    }

    double value = 0.0;
    for (IndexType k = 0; k < 3; k++) {
        IndexType l = 0; // Initializes DOF associated to N_3
        for (IndexType i = 0; i < 6; ++i) { //  For each node
            if (k == 0)
                value = (-Q1 + Q2) *  rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliar_index](i, 0);
            else if (k == 1)
                value = -(Q1 + 2.0 * Q2) * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliar_index](i, 0);
            else if (k == 2)
                value = (2.0 * Q1 + Q2) * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliar_index](i, 0);

            for (IndexType j = 0; j < 3; ++j) { // For each DOF (diagonal only)
                Kgeometricshear(n1[k] + j, l + j) += value;
                Kgeometricshear(l + j, n1[k] + j) += value;
            }

            for (IndexType j = 0; j < 3; ++j) { // For each DOF (diagonal only)
                Kgeometricshear(n2[k] + j, l + j) -= value;
                Kgeometricshear(l + j, n2[k] + j) -= value;
            }

            l += 3; // Increment DOF position I
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddBNormal(
    BoundedMatrix<double, 1, 18 >& BNormal,
    double& CNormal,
    const BoundedMatrix<double, 6, 1 >& TransversalCartesianDerivativesCenter,
    const array_1d<double, 3 >& TransversalDeformationGradientF
    )
{
    IndexType base;
    for (IndexType i = 0; i < 6; ++i) {
        base = i * 3;
        for (IndexType j = 0; j < 3; ++j) {
            BNormal(0, base + j) = TransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[j];
        }
    }

    CNormal = inner_prod(TransversalDeformationGradientF, TransversalDeformationGradientF);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddNormalKgeometric(
    BoundedMatrix<double, 36, 36 >& Kgeometricnormal,
    const BoundedMatrix<double, 6, 1 >& TransversalCartesianDerivativesCenter,
    const double SNormal
    )
{
    BoundedMatrix<double, 6, 6 > H = ZeroMatrix(6, 6);
    for (IndexType i = 0; i < 6; ++i) {
        const double aux = SNormal * TransversalCartesianDerivativesCenter(i, 0);
        for (IndexType j = 0; j < 6; ++j) {
            H(i, j) =  aux * TransversalCartesianDerivativesCenter(j, 0);
        }
    }

    noalias(H) = SNormal * prod(TransversalCartesianDerivativesCenter, trans(TransversalCartesianDerivativesCenter));

    IndexType rowindex, colindex;
    for (IndexType i = 0; i < 6; ++i) {
        rowindex = i * 3;
        for (IndexType j = 0; j < 6; ++j) {
            colindex = j * 3;
            for(IndexType ii = 0; ii < 3; ++ii) {
                Kgeometricnormal(rowindex + ii,colindex + ii) += H(i, j);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 36, 1 > SolidShellElementSprism3D6N::GetVectorCurrentPosition()
{
    KRATOS_TRY;

    BoundedMatrix<double, 36, 1 > vector_current_position;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (IndexType index = 0; index < 6; ++index) {
        const array_1d<double,3>& current_position = GetGeometry()[index].Coordinates();
        for (IndexType j = 0; j < 3; ++j)
            vector_current_position(index * 3 + j, 0) = current_position[j];
    }

    /* Neighbour nodes */
    const SizeType number_of_neighbours = NumberOfActiveNeighbours(p_neighbour_nodes);

    if (number_of_neighbours == 6) { // All the possible neighbours
        for (IndexType index = 0; index < 6; ++index) {
            const array_1d<double,3>& current_position = p_neighbour_nodes[index].Coordinates();
            for (IndexType j = 0; j < 3; ++j)
                vector_current_position(18 + index * 3 + j, 0) = current_position[j];
        }
    } else {
        for (IndexType index = 0; index < 6; ++index) {
            if (HasNeighbour(index, p_neighbour_nodes[index])) {
                const array_1d<double,3>& current_position = p_neighbour_nodes[index].Coordinates();
                for (IndexType j = 0; j < 3; ++j)
                    vector_current_position(18 + index * 3 + j, 0) = current_position[j];
            } else {
                for (IndexType j = 0; j < 3; ++j)
                    vector_current_position(18 + index * 3 + j, 0) = 0.0;
            }
        }
    }

    return vector_current_position;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 36, 1 > SolidShellElementSprism3D6N::GetVectorPreviousPosition()
{
    KRATOS_TRY;

    BoundedMatrix<double, 36, 1 > vector_current_position;

    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (IndexType index = 0; index < 6; ++index) {
        const array_1d<double,3>& previous_position = GetGeometry()[index].GetInitialPosition().Coordinates()
                                                    + GetGeometry()[index].FastGetSolutionStepValue(DISPLACEMENT, 1);
        for (IndexType j = 0; j < 3; ++j)
            vector_current_position(index * 3 + j, 0) = previous_position[j];
    }

    /* Neighbour nodes */
    const SizeType number_of_neighbours = NumberOfActiveNeighbours(p_neighbour_nodes);

    if (number_of_neighbours == 6) { // All the possible neighbours
        for (IndexType index = 0; index < 6; ++index) {
            const array_1d<double,3>& previous_position = p_neighbour_nodes[index].GetInitialPosition().Coordinates()
                                                 + p_neighbour_nodes[index].FastGetSolutionStepValue(DISPLACEMENT, 1);

            for (IndexType j = 0; j < 3; ++j)
                vector_current_position(18 + index * 3 + j, 0) = previous_position[j];
        }
    } else {
        for (IndexType index = 0; index < 6; ++index) {
            if (HasNeighbour(index, p_neighbour_nodes[index])) {
                const array_1d<double,3>& previous_position = p_neighbour_nodes[index].GetInitialPosition().Coordinates()
                                                     + p_neighbour_nodes[index].FastGetSolutionStepValue(DISPLACEMENT, 1);
                for (IndexType j = 0; j < 3; ++j)
                    vector_current_position(18 + index * 3 + j, 0) = previous_position[j];
            } else {
                for (IndexType j = 0; j < 3; ++j)
                    vector_current_position(18 + index * 3 + j, 0) = 0.0;
            }
        }
    }

    return vector_current_position;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::IntegrateStressesInZeta(
    GeneralVariables& rVariables,
    StressIntegratedComponents& rIntegratedStress,
    const double AlphaEAS,
    const double ZetaGauss,
    const double IntegrationWeight
    )
{
    KRATOS_TRY;

    const double L1 = 0.5 * (1.0 - ZetaGauss);
    const double L2 = 0.5 * (1.0 + ZetaGauss);

    const double factor_eas = std::exp(2.0 * AlphaEAS * ZetaGauss);

    /* INTEGRATE PK2 IN ZETA */
    // Integrate stresses in the reference configuration
    /* In plane stresses */
    // Lower
    rIntegratedStress.SMembraneLower(0) +=  L1 * IntegrationWeight * rVariables.StressVector[0]; // xx
    rIntegratedStress.SMembraneLower(1) +=  L1 * IntegrationWeight * rVariables.StressVector[1]; // yy
    rIntegratedStress.SMembraneLower(2) +=  L1 * IntegrationWeight * rVariables.StressVector[3]; // xy
    // Upper
    rIntegratedStress.SMembraneUpper(0) +=  L2 * IntegrationWeight * rVariables.StressVector[0]; // xx
    rIntegratedStress.SMembraneUpper(1) +=  L2 * IntegrationWeight * rVariables.StressVector[1]; // yy
    rIntegratedStress.SMembraneUpper(2) +=  L2 * IntegrationWeight * rVariables.StressVector[3]; // xy

    /* Transversal stresses */ // Note: Order according to the Voigt Notation in the Wiki
    // Lower face
    rIntegratedStress.SShearLower(0)    +=  L1 * IntegrationWeight * rVariables.StressVector[5]; // xz
    rIntegratedStress.SShearLower(1)    +=  L1 * IntegrationWeight * rVariables.StressVector[4]; // yz
    // Upper face
    rIntegratedStress.SShearUpper(0)    +=  L2 * IntegrationWeight * rVariables.StressVector[5]; // xz
    rIntegratedStress.SShearUpper(1)    +=  L2 * IntegrationWeight * rVariables.StressVector[4]; // yz

    /* Normal stress */
    rIntegratedStress.SNormal           +=  factor_eas * IntegrationWeight * rVariables.StressVector[2]; // zz

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::IntegrateEASInZeta(
    GeneralVariables& rVariables,
    EASComponents& rEAS,
    const double ZetaGauss,
    const double IntegrationWeight
    )
{
    KRATOS_TRY;

    /* INTEGRATE EAS IN ZETA */
    // Calculate EAS residual
    rEAS.mRHSAlpha += IntegrationWeight * ZetaGauss * rVariables.StressVector[2] * rVariables.C[2];

    // Calculate EAS stiffness
    rEAS.mStiffAlpha += IntegrationWeight * ZetaGauss * ZetaGauss * rVariables.C[2] * (rVariables.ConstitutiveMatrix(2, 2) * rVariables.C[2] + 2.0 * rVariables.StressVector[2]);

    BoundedMatrix<double, 1, 36 > B3;
    BoundedMatrix<double, 1,  6 > D3;

    for (IndexType i = 0; i < 6; ++i) {
        D3(0, i) = rVariables.ConstitutiveMatrix(2, i);
    }
    for (IndexType i = 0; i < 36; ++i) {
        B3(0, i) = rVariables.B(2, i);
    }

    // Calculate H operator
    noalias(rEAS.mHEAS) += IntegrationWeight * ZetaGauss * (rVariables.C[2] * prod(D3, rVariables.B) + 2.0 * rVariables.StressVector[2] * B3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    ConstitutiveLaw::Parameters& rValues,
    const StressIntegratedComponents& rIntegratedStress,
    const CommonComponents& rCommonComponents,
    const CartesianDerivatives& rCartesianDerivatives,
    const EASComponents& rEAS,
    double& AlphaEAS
    )
{
    /* Contributions of the stiffness matrix calculated on the reference configuration */
    if( rLocalSystem.CalculationFlags.Is( SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) ) {
        std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
        const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

        for( IndexType i = 0; i < rLeftHandSideVariables.size(); ++i ) {
            bool calculated = false;
            /* Calculate the Material Stiffness Matrix */
            if( rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX ) {
                /* Reading integration points */
                const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

                for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
                    const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

                    /* Assemble B */
                    this->CalculateDeformationMatrix(rVariables.B, rCommonComponents, zeta_gauss, AlphaEAS);

                    // Compute element kinematics C, F ...
                    this->CalculateKinematics(rVariables, rCommonComponents, integration_points, point_number, AlphaEAS, zeta_gauss);

                    // Set general variables to constitutivelaw parameters
                    this->SetGeneralVariables(rVariables, rValues, point_number);

                    // Compute stresses and constitutive parameters
                    mConstitutiveLawVector[point_number]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

                    // Calculating weights for integration on the "reference configuration"
                    const double integration_weight = integration_points[point_number].Weight() * rVariables.detJ;

                    /* Operation performed: add Km to the LefsHandSideMatrix */
                    this->CalculateAndAddKuum( rLeftHandSideMatrices[i], rVariables, integration_weight);
                }
                calculated = true;
            }

            /* Calculate the Geometric Stiffness Matrix */
            if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX ) {
                /* Operation performed: add Kg to the LefsHandSideMatrix */
                this->CalculateAndAddKuug( rLeftHandSideMatrices[i], rIntegratedStress, rCartesianDerivatives );
                calculated = true;
            }

            /* Implicit or explicit EAS update*/
            if ( mELementalFlags.Is(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT)) {
                /* Apply EAS stabilization */
                ApplyEASLHS(rLeftHandSideMatrices[i], rEAS);
            }

            KRATOS_ERROR_IF_NOT(calculated) << " ELEMENT can not supply the required local system variable: " << rLeftHandSideVariables[i] << std::endl;
        }
    } else {
        MatrixType& LeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        /* Calculate the Material Stiffness Matrix */
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            /* Assemble B */
            this->CalculateDeformationMatrix(rVariables.B, rCommonComponents, zeta_gauss, AlphaEAS);

            // Compute element kinematics C, F ...
            this->CalculateKinematics(rVariables, rCommonComponents, integration_points, point_number, AlphaEAS, zeta_gauss);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(rVariables, rValues, point_number);

            // Compute stresses and constitutive parameters
            mConstitutiveLawVector[point_number]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

            // Calculating weights for integration on the "reference configuration"
            const double integration_weight = integration_points[point_number].Weight() * rVariables.detJ;

            /* Operation performed: add Km to the LefsHandSideMatrix */
            this->CalculateAndAddKuum( LeftHandSideMatrix, rVariables, integration_weight);
        }

        /* Calculate the Geometric Stiffness Matrix */
        /* Operation performed: add Kg to the LefsHandSideMatrix */
        this->CalculateAndAddKuug( LeftHandSideMatrix, rIntegratedStress, rCartesianDerivatives );

        /* Implicit or explicit EAS update*/
        if ( mELementalFlags.Is(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT)) {
            /* Apply EAS stabilization */
            ApplyEASLHS(LeftHandSideMatrix, rEAS);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    Vector& rVolumeForce,
    const StressIntegratedComponents& rIntegratedStress,
    const CommonComponents& rCommonComponents,
    const EASComponents& rEAS,
    double& AlphaEAS
    )
{
    /* Contribution of the internal and external forces */
    if( rLocalSystem.CalculationFlags.Is( SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) ) {
        std::vector<VectorType>& RightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
        const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
        for( IndexType i = 0; i < rRightHandSideVariables.size(); ++i ) {
            bool calculated = false;
            if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR ) {
                /* Operation performed: RightHandSideVector += ExtForce */
                this->CalculateAndAddExternalForces( RightHandSideVectors[i], rVariables, rVolumeForce );
                calculated = true;
            }

            if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR ) {
                /* Operation performed: RightHandSideVector -= IntForce */
                this->CalculateAndAddInternalForces( RightHandSideVectors[i], rIntegratedStress, rCommonComponents, rEAS, AlphaEAS );
                calculated = true;
            }

            KRATOS_ERROR_IF_NOT(calculated) << " ELEMENT can not supply the required local system variable: " << rRightHandSideVariables[i] << std::endl;
        }
    } else {
        VectorType& RightHandSideVector = rLocalSystem.GetRightHandSideVector();

        /* Operation performed: RightHandSideVector += ExtForce */
        this->CalculateAndAddExternalForces( RightHandSideVector, rVariables, rVolumeForce );

        /* Operation performed: RightHandSideVector -= IntForce */
        this->CalculateAndAddInternalForces( RightHandSideVector, rIntegratedStress, rCommonComponents, rEAS, AlphaEAS );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddKuum(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double IntegrationWeight
    )
{
    KRATOS_TRY;

    /* Calculate K */
    typedef BoundedMatrix<double,  6, 36 > temp_type;
    const BoundedMatrix<double, 36, 36 > K = IntegrationWeight * prod(trans(rVariables.B), prod<temp_type>(rVariables.ConstitutiveMatrix, rVariables.B));

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    IndexType index_i, index_j;

    for (IndexType i = 0; i < 36; ++i) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (index_i < 36) {
            for (IndexType j = 0; j < 36; ++j) {
                index_j = j < 18 ? j : id_vector[j - 18];
                if (index_j < 36) {
                    rLeftHandSideMatrix(index_i, index_j) += K(i, j);
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddKuug(
    MatrixType& rLeftHandSideMatrix,
    const StressIntegratedComponents& rIntegratedStress,
    const CartesianDerivatives& rCartesianDerivatives
    )
{
    KRATOS_TRY;

    /* The stress is already integrated, we just calculate it once */

    /* Auxiliar stiffness matrix */
    BoundedMatrix<double, 36, 36 > K = ZeroMatrix(36, 36); // Stiffness matrix

    /* COMPUTATION OF GEOMETRIC STIFFNESS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Adding the geometric membrane stiffness */
    // Lower face
    CalculateAndAddMembraneKgeometric(K, rCartesianDerivatives, rIntegratedStress.SMembraneLower, GeometricLevel::LOWER);
    // Upper face
    CalculateAndAddMembraneKgeometric(K, rCartesianDerivatives, rIntegratedStress.SMembraneUpper, GeometricLevel::UPPER);

//    /* SHEAR CONTRIBUTION */
//    /* Adding the geometric shear stiffness */
//    // Lower face
//    CalculateAndAddShearKgeometric(K, rCartesianDerivatives, rIntegratedStress.SShearLower, GeometricLevel::LOWER);
//    // Upper face
//    CalculateAndAddShearKgeometric(K, rCartesianDerivatives, rIntegratedStress.SShearUpper, GeometricLevel::UPPER);

    /* NORMAL TRANSVERSE */
    /* Adding the geometric normal stiffness */
    CalculateAndAddNormalKgeometric(K, rCartesianDerivatives.TransversalCartesianDerivativesCenter, rIntegratedStress.SNormal);

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    IndexType index_i, index_j;
    for (IndexType i = 0; i < 36; ++i) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (index_i < 36) {
            for (IndexType j = 0; j < 36; ++j) {
                index_j = j < 18 ? j : id_vector[j - 18];
                if (index_j < 36) {
                    rLeftHandSideMatrix(index_i, index_j) += K(i, j);
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::ApplyEASLHS(
    MatrixType& rLeftHandSideMatrix,
    const EASComponents& rEAS
    )
{
    KRATOS_TRY;

    const BoundedMatrix<double, 36, 36 > lhs_aux = - prod(trans(rEAS.mHEAS), rEAS.mHEAS) / rEAS.mStiffAlpha;

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    // Note: IntegrationWeight already considered in the integration
    IndexType index_i, index_j;
    for (IndexType i = 0; i < 36; ++i) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (index_i < 36) {
            for (IndexType j = 0; j < 36; ++j) {
                index_j = j < 18 ? j : id_vector[j - 18];
                if (index_j < 36) {
                    rLeftHandSideMatrix(index_i, index_j) += lhs_aux(i, j);
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::ApplyEASRHS(
    BoundedMatrix<double, 36, 1 >& rRHSFull,
    const EASComponents& rEAS,
    double& AlphaEAS
    )
{
    KRATOS_TRY;

    /* Calculate the RHS */
    noalias(rRHSFull) -= trans(rEAS.mHEAS) * rEAS.mRHSAlpha / rEAS.mStiffAlpha;

    /* Update ALPHA_EAS */
    AlphaEAS -= rEAS.mRHSAlpha / rEAS.mStiffAlpha;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddExternalForces(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    Vector& rVolumeForce
    )
{
    KRATOS_TRY;

    const IndexType number_of_nodes = GetGeometry().PointsNumber();

    IndexType index;
    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        index = 3 * i;
        for ( IndexType j = 0; j < 3; ++j ) {
            rRightHandSideVector[index + j] += rVolumeForce[j]/static_cast<double>(number_of_nodes);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddInternalForces(
    VectorType& rRightHandSideVector,
    const StressIntegratedComponents& rIntegratedStress,
    const CommonComponents& rCommonComponents,
    const EASComponents& rEAS,
    double& AlphaEAS
    )
{
    KRATOS_TRY;

    BoundedMatrix<double, 36, 1 > rhs_full = ZeroMatrix(36, 1);

    IndexType aux_index = 0;
    for (IndexType i = 0; i < 18; ++i) {
        if (i == 9) aux_index += 9;

        /* Calculate residual forces */
        /* Apply membrane stress, adding the in-plane nodal force contribution */
        /* Nodes 1-3  and 7-9 */
        rhs_full(aux_index + i, 0)     += rIntegratedStress.SMembraneLower[0] * rCommonComponents.BMembraneLower(0, i); // xx
        rhs_full(aux_index + i, 0)     += rIntegratedStress.SMembraneLower[1] * rCommonComponents.BMembraneLower(1, i); // yy
        rhs_full(aux_index + i, 0)     += rIntegratedStress.SMembraneLower[2] * rCommonComponents.BMembraneLower(2, i); // xy

        /* Nodes 4-6  and 10-12 */
        rhs_full(aux_index + i + 9, 0) += rIntegratedStress.SMembraneUpper[0] * rCommonComponents.BMembraneUpper(0, i); // xx
        rhs_full(aux_index + i + 9, 0) += rIntegratedStress.SMembraneUpper[1] * rCommonComponents.BMembraneUpper(1, i); // yy
        rhs_full(aux_index + i + 9, 0) += rIntegratedStress.SMembraneUpper[2] * rCommonComponents.BMembraneUpper(2, i); // xy

        /* Apply transversal forces */
        /* Apply shear stress, adding the transverse nodal force contribution */
        rhs_full(i, 0) += rIntegratedStress.SShearLower[0] * rCommonComponents.BShearLower(0, i); // xz
        rhs_full(i, 0) += rIntegratedStress.SShearLower[1] * rCommonComponents.BShearLower(1, i); // yz
        rhs_full(i, 0) += rIntegratedStress.SShearUpper[0] * rCommonComponents.BShearUpper(0, i); // xz
        rhs_full(i, 0) += rIntegratedStress.SShearUpper[1] * rCommonComponents.BShearUpper(1, i); // yz

        /* Apply normal transverse stress */
        rhs_full(i, 0) += rIntegratedStress.SNormal * rCommonComponents.BNormal(0, i); // zz
    }

    /* Apply EAS stabilization */
    ApplyEASRHS(rhs_full, rEAS, AlphaEAS);

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    IndexType index_i;
    for (IndexType i = 0; i < 36; ++i) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (index_i < 36) {
            rRightHandSideVector[index_i] -= rhs_full(i, 0);
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetGeneralVariables(
    GeneralVariables& rVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType rPointNumber
    )
{
    KRATOS_ERROR_IF(rVariables.detF < 0) << "SPRISM ELEMENT: " << this->Id() << "INVERTED: |F| < 0  detF = " << rVariables.detF << std::endl;

    // Compute total F: FT
    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );

    rValues.SetDeterminantF(rVariables.detFT);
    rValues.SetDeformationGradientF(rVariables.FT);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);

    // Adding the standard prism shape functions
    rValues.SetShapeFunctionsValues(rVariables.N);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::InitializeSystemMatrices(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags
    )
{
    // Resizing as needed the LHS
    WeakPointerVectorNodesType& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const IndexType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);
    const IndexType mat_size = number_of_nodes * 3;

    if ( rCalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX) ) {// Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR) ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); // Resetting RHS
    }
}

/******************************* COMPUTE KINEMATICS ********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateKinematics(
    GeneralVariables& rVariables,
    const CommonComponents& rCommonComponents,
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType rPointNumber,
    const double AlphaEAS,
    const double ZetaGauss
    )
{
    KRATOS_TRY;

    const double L_1 = 0.5 * (1.0 - ZetaGauss);
    const double L_2 = 0.5 * (1.0 + ZetaGauss);

    const double factor_eas = std::exp(2.0 * AlphaEAS * ZetaGauss);  // EAS factor

    /* Assemble C */
    rVariables.C[0] = L_1 * rCommonComponents.CMembraneLower(0, 0) + L_2 * rCommonComponents.CMembraneUpper(0, 0); // xx
    rVariables.C[1] = L_1 * rCommonComponents.CMembraneLower(1, 0) + L_2 * rCommonComponents.CMembraneUpper(1, 0); // yy
    rVariables.C[2] = factor_eas * rCommonComponents.CNormal;                                            // zz
    rVariables.C[3] = L_1 * rCommonComponents.CMembraneLower(2, 0) + L_2 * rCommonComponents.CMembraneUpper(2, 0); // xy
    rVariables.C[4] = L_1 * rCommonComponents.CShearLower(1, 0)    + L_2 * rCommonComponents.CShearUpper(1, 0);    // yz
    rVariables.C[5] = L_1 * rCommonComponents.CShearLower(0, 0)    + L_2 * rCommonComponents.CShearUpper(0, 0);    // xz

    rVariables.detF = rVariables.C[0] * rVariables.C[1] * rVariables.C[2] + 2 * rVariables.C[3] * rVariables.C[4] * rVariables.C[5]
                    - rVariables.C[5] * rVariables.C[5] * rVariables.C[1] -     rVariables.C[4] * rVariables.C[4] * rVariables.C[0]
                    - rVariables.C[3] * rVariables.C[3] * rVariables.C[2];


    KRATOS_ERROR_IF(rVariables.detF < std::numeric_limits<double>::epsilon()) << "The determinant of C is zero or negative.  det(C): " << rVariables.detF << std::endl;

    rVariables.detF = std::sqrt(rVariables.detF);

    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        // PK2 stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

        // Jacobian Determinant for the isoparametric and numerical integration
        Matrix J0;
        GeometryUtils::JacobianOnInitialConfiguration(GetGeometry(), rIntegrationPoints[rPointNumber], J0);
        rVariables.detJ = MathUtils<double>::DetMat(J0);
    } else {
        // Cauchy stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

        //Determinant of the Deformation Gradient F0
        rVariables.detF0 = MathUtils<double>::Det3(mAuxContainer[rPointNumber]);
        rVariables.F0    = mAuxContainer[rPointNumber];
    }

    this->CbartoFbar(rVariables, rPointNumber);

    // Get the shape functions for the order of the integration method [N]
    const Matrix& N_container = rVariables.GetShapeFunctions();

    // Set Shape Functions Values for this integration point
    rVariables.N = row( N_container, rPointNumber);

    KRATOS_CATCH( "" );
}

/***************************** COMPUTE DELTA POSITION ******************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY;

    for ( IndexType i = 0; i < 6; ++i ) {
        const array_1d<double, 3 > & current_displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & previous_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, 1);

        for ( IndexType j = 0; j < 3; ++j )
            rDeltaPosition(i,j) = current_displacement[j] - previous_displacement[j];
    }

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CbartoFbar(
    GeneralVariables& rVariables,
    const int rPointNumber
    )
{
    KRATOS_TRY;

    /* We perform a polar decomposition of the CBar and F(regular) to obtain F_bar */

    /* Decompose C_bar */
    BoundedMatrix<double, 3, 3> eigen_vector_matrix,  eigen_values_matrix;

    // Assemble matrix C_bar
    const Matrix C_bar = MathUtils<double>::VectorToSymmetricTensor(rVariables.C);

    // Decompose matrix C_bar
    MathUtils<double>::EigenSystem<3>(C_bar, eigen_vector_matrix, eigen_values_matrix, 1e-24, 100);

    for (IndexType i = 0; i < 3; ++i)
        eigen_values_matrix(i, i) = std::sqrt(eigen_values_matrix(i, i));

    const Matrix U_bar = prod( eigen_values_matrix, eigen_vector_matrix );

    /* Decompose F */
    Matrix F = ZeroMatrix(3, 3);
    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[rPointNumber], mAuxContainer[rPointNumber] );
    } else {
        // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
        Matrix InvJ(3, 3);
        MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[rPointNumber], InvJ );
    }

    // Compute R
    Matrix R(3, 3);
    Matrix U(3, 3);
    ConstitutiveLawUtilities<6>::PolarDecomposition(F, R, U);

    /* Calculate F_bar */
    noalias(rVariables.F) = prod(R, U_bar);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateDeformationMatrix(
    Matrix& rB,
    const CommonComponents& rCommonComponents,
    const double ZetaGauss,
    const double AlphaEAS
    )
{
    KRATOS_TRY;

    rB.clear(); // Set all components to zero

    const double L1 = 0.5 * (1.0 - ZetaGauss);
    const double L2 = 0.5 * (1.0 + ZetaGauss);

    const double factor_eas = std::exp(2.0 * AlphaEAS * ZetaGauss); // EAS factor

    for (IndexType index = 0; index < 9; ++index) {
        /* Element nodes */ // Note: It's important to consider the Voigt notation order considered in Kratos
        // Lower face
        rB(0, index)      = L1 * rCommonComponents.BMembraneLower(0, index);  // xx
        rB(1, index)      = L1 * rCommonComponents.BMembraneLower(1, index);  // yy
        rB(2, index)      = factor_eas * rCommonComponents.BNormal(0, index);     // zz
        rB(3, index)      = L1 * rCommonComponents.BMembraneLower(2, index);  // xy
        rB(4, index)      = L1 * rCommonComponents.BShearLower(1, index) + L2 * rCommonComponents.BShearUpper(1, index); // yz
        rB(5, index)      = L1 * rCommonComponents.BShearLower(0, index) + L2 * rCommonComponents.BShearUpper(0, index); // xz
        // Upper face
        rB(0, index + 9)  = L2 * rCommonComponents.BMembraneUpper(0, index);  // xx
        rB(1, index + 9)  = L2 * rCommonComponents.BMembraneUpper(1, index);  // yy
        rB(2, index + 9)  = factor_eas * rCommonComponents.BNormal(0, index + 9); // zz
        rB(3, index + 9)  = L2 * rCommonComponents.BMembraneUpper(2, index);  // xy
        rB(4, index + 9)  = L1 * rCommonComponents.BShearLower(1, index + 9) + L2 * rCommonComponents.BShearUpper(1, index + 9); // yz
        rB(5, index + 9)  = L1 * rCommonComponents.BShearLower(0, index + 9) + L2 * rCommonComponents.BShearUpper(0, index + 9); // xz

        /* Neighbour nodes */
        // Lower face
        rB(0, index + 18) = L1 * rCommonComponents.BMembraneLower(0, index + 9); // xx
        rB(1, index + 18) = L1 * rCommonComponents.BMembraneLower(1, index + 9); // yy
        rB(3, index + 18) = L1 * rCommonComponents.BMembraneLower(2, index + 9); // xy
        // Upper face
        rB(0, index + 27) = L2 * rCommonComponents.BMembraneUpper(0, index + 9); // xx
        rB(1, index + 27) = L2 * rCommonComponents.BMembraneUpper(1, index + 9); // yy
        rB(3, index + 27) = L2 * rCommonComponents.BMembraneUpper(2, index + 9); // xy
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::InitializeGeneralVariables(GeneralVariables& rVariables)
{
    // StressMeasure_PK1             //stress related to reference configuration non-symmetric
    // StressMeasure_PK2             //stress related to reference configuration
    // StressMeasure_Kirchhoff       //stress related to current   configuration
    // StressMeasure_Cauchy          //stress related to current   configuration

    // StressMeasure
    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN))
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;
    else
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Doubles
    rVariables.detF  = 1.0;
    rVariables.detF0 = 1.0;
    rVariables.detFT = 1.0;
    rVariables.detJ  = 1.0;

    // Vectors
    rVariables.StrainVector = ZeroVector(6);
    rVariables.StressVector = ZeroVector(6);
    rVariables.C = ZeroVector(6);
    rVariables.N = ZeroVector(6);

    // Matrices
    rVariables.F  = IdentityMatrix(3);
    rVariables.F0 = IdentityMatrix(3);
    rVariables.FT = IdentityMatrix(3);
    rVariables.B  = ZeroMatrix(6, 36);

    rVariables.DN_DX = ZeroMatrix(6, 3);
    rVariables.ConstitutiveMatrix = ZeroMatrix(6, 6);

    // Reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( this->GetIntegrationMethod() ));

    // Jacobians
    rVariables.J.resize(1, false);
    rVariables.j.resize(1, false);
    rVariables.J[0] = ZeroMatrix(1, 1);
    rVariables.j[0] = ZeroMatrix(1, 1);

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, this->GetIntegrationMethod() );

    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN) == false ) {
        //Calculate Delta Position
        Matrix delta_position( 6 , 3);
        this->CalculateDeltaPosition(delta_position);
        rVariables.J = GetGeometry().Jacobian( rVariables.J, this->GetIntegrationMethod(), delta_position);
    }

    // Computing gradient
    const IndexType integration_point_number = GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );
    GeometryType::ShapeFunctionsGradientsType DN_DX(integration_point_number, ZeroMatrix(6, 3));
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    double detJ;
    Matrix inv_j;
    for (IndexType i_point = 0; i_point < integration_point_number; ++i_point) {
        MathUtils<double>::InvertMatrix( rVariables.j[i_point], inv_j, detJ );
        noalias(DN_DX[i_point]) = prod(DN_De[i_point], inv_j);
    }
    rVariables.SetShapeFunctionsGradients(DN_DX);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::FinalizeStepVariables(
    GeneralVariables & rVariables,
    const IndexType rPointNumber
    )
{
    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN) == false ) {
        // Update internal (historical) variables
        mAuxContainer[rPointNumber] = prod(rVariables.F, rVariables.F0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetHistoricalVariables(
    GeneralVariables& rVariables,
    const IndexType rPointNumber
    )
{
    /* Deformation Gradient F ( set to identity ) */
    const IndexType size =  rVariables.F.size1();

    rVariables.detF  = 1.0;
    rVariables.F     = IdentityMatrix(size);
}

/**************************** CALCULATE VOLUME CHANGE ******************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY;

    if ( mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
        rVolumeChange = 1.0;
    else
        rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    KRATOS_CATCH( "" );
}

/************************* CALCULATE VOLUME ACCELERATION ***************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateVolumeForce(
    Vector& rVolumeForce,
    GeneralVariables& rVariables,
    const double IntegrationWeight
    )
{
    KRATOS_TRY;

    array_1d<double,3> volume_acceleration = ZeroVector(3);
    if (GetProperties().Has( VOLUME_ACCELERATION ))
        volume_acceleration = GetProperties()[VOLUME_ACCELERATION];
    else if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        for (unsigned int i_node = 0; i_node < this->GetGeometry().size(); ++i_node)
            volume_acceleration += rVariables.N[i_node] * GetGeometry()[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    // Compute volume change
    double volume_change;
    this->CalculateVolumeChange( volume_change, rVariables );

    rVolumeForce += volume_acceleration * IntegrationWeight * volume_change * GetProperties()[DENSITY];

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseSolidElement);
    rSerializer.save("FinalizedStep",mFinalizedStep);
    rSerializer.save("HistoricalF0",mAuxContainer);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseSolidElement);
    rSerializer.load("FinalizedStep",mFinalizedStep);
    rSerializer.load("HistoricalF0",mAuxContainer);
}

} // Namespace Kratos.
