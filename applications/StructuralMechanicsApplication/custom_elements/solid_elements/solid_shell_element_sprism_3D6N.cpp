// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/global_pointer_variables.h"
#include "input_output/logger.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "solid_shell_element_sprism_3D6N.h"

namespace Kratos
{
/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, EAS_IMPLICIT_EXPLICIT,              4 ); // True means implicit
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, TOTAL_UPDATED_LAGRANGIAN,           5 ); // True means total lagrangian // TODO: change this using templates!!!
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, QUADRATIC_ELEMENT,                  6 ); // True means quadratic in-plane behaviour
KRATOS_CREATE_LOCAL_FLAG( SolidShellElementSprism3D6N, EXPLICIT_RHS_COMPUTATION,           7 ); // True means elastic behaviour for stabilization

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

/********************************** ASSIGNMENT OPERATOR ****************************/
/***********************************************************************************/

SolidShellElementSprism3D6N&  SolidShellElementSprism3D6N::operator=(SolidShellElementSprism3D6N const& rOther)
{
    BaseType::operator=(rOther);

    mIntegrationOrder = rOther.mIntegrationOrder;

    mAuxContainer.clear();
    mAuxContainer.resize(rOther.mAuxContainer.size());

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
        mAuxContainer[i] = rOther.mAuxContainer[i];
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
    return Kratos::make_intrusive<SolidShellElementSprism3D6N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//************************************************************************************
//************************************************************************************

Element::Pointer SolidShellElementSprism3D6N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SolidShellElementSprism3D6N>(NewId, pGeom, pProperties);
}

/*********************************** CLONE ******************************************/
/************************************************************************************/

Element::Pointer SolidShellElementSprism3D6N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    SolidShellElementSprism3D6N new_element( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    new_element.mIntegrationOrder = mIntegrationOrder;

    const unsigned int integration_point_number = mConstitutiveLawVector.size();

    if ( new_element.mConstitutiveLawVector.size() != integration_point_number) {
        new_element.mConstitutiveLawVector.resize(integration_point_number);
    }

    for(unsigned int i = 0; i < integration_point_number; ++i) {
        new_element.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    if (new_element.mAuxContainer.size() != mAuxContainer.size()) {
        new_element.mAuxContainer.resize(mAuxContainer.size());
    }

    for(unsigned int i = 0; i < mAuxContainer.size(); ++i) {
        new_element.mAuxContainer[i] = mAuxContainer[i];
    }

    return Kratos::make_intrusive<SolidShellElementSprism3D6N>(new_element);
}

//******************************* GETTING METHODS *********************************//
/***********************************************************************************/

void SolidShellElementSprism3D6N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Get the neighbours
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int number_of_nodes = r_geometry.size() + NumberOfActiveNeighbours(r_neighbour_nodes);
    const unsigned int dim = number_of_nodes * 3;

    if (rResult.size() != dim) {
        rResult.resize(dim, false);
    }

    // Nodes of the central element
    unsigned int index = 0;
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        rResult[index]     = r_geometry[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z).EquationId();
        index += 3;
    }

    // Adding the ids of the neighbouring nodes
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        if (HasNeighbour(i, r_neighbour_nodes[i])) {
            rResult[index]     = r_neighbour_nodes[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_neighbour_nodes[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_neighbour_nodes[i].GetDof(DISPLACEMENT_Z).EquationId();
            index += 3;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    // Get geometry
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();

    // Get the neighbours
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    rElementalDofList.resize(0);

    // Nodes of the central element
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        rElementalDofList.push_back(r_geometry[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_geometry[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(r_geometry[i].pGetDof(DISPLACEMENT_Z));
    }

    // Adding the dofs of the neighbouring nodes
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        if (HasNeighbour(i, r_neighbour_nodes[i])) {
            rElementalDofList.push_back(r_neighbour_nodes[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_neighbour_nodes[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_neighbour_nodes[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("");
}

/******************************** DISPLACEMENT **************************************/
/************************************************************************************/

void SolidShellElementSprism3D6N::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Get the neighbours
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int number_of_nodes = r_geometry.size() + NumberOfActiveNeighbours(r_neighbour_nodes);

    const unsigned int mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    unsigned int index = 0;

    // Nodes of the central element
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        const array_1d<double, 3>& r_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        for (unsigned int j = 0; j < 3; ++j) {
            rValues[index + j] = r_displacement[j];
        }
        index += 3;
    }

    // Neighbour nodes
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        if (HasNeighbour(i, r_neighbour_nodes[i])) {
            const array_1d<double, 3>& r_displacement = r_neighbour_nodes[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            for (unsigned int j = 0; j < 3; ++j) {
                rValues[index + j] = r_displacement[j];
            }
            index += 3;
        }
    }
}

/********************************** VELOCITY ****************************************/
/************************************************************************************/

void SolidShellElementSprism3D6N::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Get the neighbours
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int number_of_nodes = r_geometry.size() + NumberOfActiveNeighbours(r_neighbour_nodes);

    const unsigned int mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    unsigned int index = 0;

    // Nodes of the central element
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        const array_1d<double, 3>& r_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int j = 0; j < 3; ++j) {
            rValues[index + j] = r_velocity[j];
        }
        index += 3;
    }

    // Neighbour nodes
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        if (HasNeighbour(i, r_neighbour_nodes[i])) {
            const array_1d<double, 3>& r_velocity = r_neighbour_nodes[i].FastGetSolutionStepValue(VELOCITY, Step);
            for (unsigned int j = 0; j < 3; ++j) {
                rValues[index + j] = r_velocity[j];
            }
            index += 3;
        }
    }
}

/******************************** ACCELERATION **************************************/
/************************************************************************************/

void SolidShellElementSprism3D6N::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Get the neighbours
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int number_of_nodes = r_geometry.size() + NumberOfActiveNeighbours(r_neighbour_nodes);

    const unsigned int mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    unsigned int index = 0;

    // Nodes of the central element
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        const array_1d<double, 3>& r_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        for (unsigned int j = 0; j < 3; ++j) {
            rValues[index + j] = r_acceleration[j];
        }
        index += 3;
    }

    // Neighbour nodes
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        if (HasNeighbour(i, r_neighbour_nodes[i])) {
            const array_1d<double, 3>& r_acceleration = r_neighbour_nodes[i].FastGetSolutionStepValue(ACCELERATION, Step);
            for (unsigned int j = 0; j < 3; ++j) {
                rValues[index + j] = r_acceleration[j];
            }
            index += 3;
        }
    }
}

//****************************** COMPUTING METHODS ********************************//
/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR);

    MatrixType left_hand_side_matrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices(left_hand_side_matrix, rRightHandSideVector, local_system.CalculationFlags);

    //Set general_variables to Local system components
    local_system.SetLeftHandSideMatrix(left_hand_side_matrix);
    local_system.SetRightHandSideVector(rRightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
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
    const ProcessInfo& rCurrentProcessInfo
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

void SolidShellElementSprism3D6N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Get geometry
    const auto& r_geometry = GetGeometry();

    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    const unsigned int number_of_nodes = r_geometry.size() + NumberOfActiveNeighbours(r_neighbour_nodes);
    const unsigned int mat_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != mat_size) {
        rMassMatrix.resize(mat_size, mat_size, false);
    }

    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

    Matrix aux_matrix;
    const unsigned int aux_mat_size = r_geometry.size() * 3;
    BaseType::CalculateMassMatrix(aux_matrix, rCurrentProcessInfo);
    noalias(subrange(rMassMatrix, 0, aux_mat_size, 0, aux_mat_size)) = aux_matrix;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int mat_size = ( GetGeometry().size() + NumberOfActiveNeighbours(r_neighbour_nodes) ) * 3;

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        mat_size);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const MatrixType& rStiffnessMatrix,
    const MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(r_neighbour_nodes);

    // Resizing as needed the LHS
    const unsigned int mat_size = number_of_nodes * 3;

    if ( rDampingMatrix.size1() != mat_size ) {
        rDampingMatrix.resize( mat_size, mat_size, false );
    }

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Get Damping Coefficients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
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
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Get geometry
    const auto& r_geometry = GetGeometry();

    const unsigned int integration_point_number = msGeometryData.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number ) {
        rOutput.resize( integration_point_number );
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for (unsigned int ii = 0; ii < integration_point_number; ++ii ) {
            bool aux_bool;
            aux_bool = mConstitutiveLawVector[ii]->GetValue( rVariable, aux_bool);
            rOutput[ii] = aux_bool;
        }
    } else {
        /* Create and initialize element variables: */
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags& r_constitutive_law_options = cl_values.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if (mFinalizedStep) {
                this->GetHistoricalVariables(general_variables,point_number);
            }

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(general_variables, cl_values, point_number);

            bool aux_bool;
            aux_bool = mConstitutiveLawVector[point_number]->CalculateValue( cl_values, rVariable, aux_bool );
            rOutput[point_number] = aux_bool;
        }
    }

    if ( rOutput.size() != 6 ) {
        std::vector<bool> r_output_aux;
        r_output_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; ++iii) {
            rOutput[iii] = false;

            for (unsigned int i_gp = 0; i_gp < integration_point_number; i_gp++) {
                if (r_output_aux[i_gp]) {
                    rOutput[iii] = true;
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Get the integration points number
    const unsigned int integration_point_number = msGeometryData.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number ) {
        rOutput.resize( integration_point_number );
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    if ( rOutput.size() != 6 ) {
        std::vector<int> r_output_aux;
        r_output_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; ++iii) {
            rOutput[iii] = 0;

            for (unsigned int i_gp = 0; i_gp < integration_point_number; i_gp++) {
                rOutput[iii] += interpol(i_gp, iii) * r_output_aux[i_gp];
            }
        }
    }

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

    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Get the integration points number
    const unsigned int integration_point_number = msGeometryData.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number )
        rOutput.resize( integration_point_number );

    if ( rVariable == VON_MISES_STRESS ) {
        /* Create and initialize element variables: */
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags& r_constitutive_law_options = cl_values.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if (mFinalizedStep) {
                this->GetHistoricalVariables(general_variables,point_number);
            }

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(general_variables, cl_values, point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy (cl_values);

            const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(general_variables.StressVector); //reduced dimension stress tensor


            // In general coordinates:
            double sigma_equivalent =  (0.5)*((stress_tensor(0,0)-stress_tensor(1,1))*((stress_tensor(0,0)-stress_tensor(1,1)))+
                                              (stress_tensor(1,1)-stress_tensor(2,2))*((stress_tensor(1,1)-stress_tensor(2,2)))+
                                              (stress_tensor(2,2)-stress_tensor(0,0))*((stress_tensor(2,2)-stress_tensor(0,0)))+
                                            6*(stress_tensor(0,1)*stress_tensor(1,0)+stress_tensor(1,2)*stress_tensor(2,1)+stress_tensor(2,0)*stress_tensor(0,2)));

            if (sigma_equivalent < 0) {
                sigma_equivalent = 0;
            }

            sigma_equivalent = std::sqrt(sigma_equivalent);

            rOutput[point_number] =  sigma_equivalent;
        }
    } else if ( rVariable == NORM_ISOCHORIC_STRESS ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& r_constitutive_law_options = cl_values.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY, true);

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points,point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if (mFinalizedStep) {
                this->GetHistoricalVariables(general_variables,point_number);
            }

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(general_variables, cl_values, point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy (cl_values);

            const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(general_variables.StressVector); //reduced dimension stress tensor

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
        ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& r_constitutive_law_options = cl_values.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY, true);

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& r_alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if (mFinalizedStep) {
                this->GetHistoricalVariables(general_variables, point_number);
            }

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(general_variables, cl_values, point_number);

            double strain_energy = 0.0;

            // Compute stresses and constitutive parameters
            if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
                mConstitutiveLawVector[point_number]->CalculateMaterialResponseKirchhoff(cl_values);
            } else {
                mConstitutiveLawVector[point_number]->CalculateMaterialResponsePK2(cl_values);
            }

            mConstitutiveLawVector[point_number]->GetValue(STRAIN_ENERGY, strain_energy);

            rOutput[point_number] = general_variables.detJ * r_integration_points[point_number].Weight() * strain_energy;  // 1/2 * sigma * epsilon
        }
    } else if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
       CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    if ( rOutput.size() != 6 ) {
        std::vector<double> r_output_aux;
        r_output_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; ++iii) {
            rOutput[iii] = 0.0;

            for (unsigned int i_gp = 0; i_gp < integration_point_number; i_gp++) {
                rOutput[iii] += interpol(i_gp, iii) * r_output_aux[i_gp];
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Get the integration points
    const unsigned int integration_point_number = msGeometryData.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number ) {
        rOutput.resize( integration_point_number );
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    if ( rOutput.size() != 6 ) {
        std::vector<array_1d<double, 3>> r_output_aux;
        r_output_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; ++iii) {
            rOutput[iii] = ZeroVector(3);

            for (unsigned int i_gp = 0; i_gp < integration_point_number; i_gp++) {
                rOutput[iii] += interpol(i_gp, iii) * r_output_aux[i_gp];
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Get integration points
    const unsigned int integration_point_number = msGeometryData.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if (rOutput.size() != integration_point_number ) {
        rOutput.resize( integration_point_number );
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    if ( rOutput.size() != 6 ) {
        std::vector<array_1d<double, 6>> r_output_aux;
        r_output_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; ++iii) {
            rOutput[iii] = ZeroVector(6);

            for (unsigned int i_gp = 0; i_gp < integration_point_number; i_gp++) {
                rOutput[iii] += interpol(i_gp, iii) * r_output_aux[i_gp];
            }
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

    // Get geometry
    const auto& r_geometry = GetGeometry();

    const IndexType integration_point_number = msGeometryData.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number ) {
        rOutput.resize( integration_point_number );
    }

    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& r_constitutive_laws_options = cl_values.GetOptions();
        r_constitutive_laws_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_constitutive_laws_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& r_alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if (mFinalizedStep) {
                this->GetHistoricalVariables(general_variables,point_number);
            }

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(general_variables, cl_values, point_number);

            // Call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR)
                general_variables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
            else
                general_variables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

            mConstitutiveLawVector[point_number]->CalculateMaterialResponse(cl_values, general_variables.StressMeasure);

            if (rOutput[point_number].size() != general_variables.StressVector.size()) {
                rOutput[point_number].resize( general_variables.StressVector.size(), false);
            }
            rOutput[point_number] = general_variables.StressVector;
        }
    } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR || rVariable == HENCKY_STRAIN_VECTOR) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags& r_constitutive_law_options=cl_values.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Set the strain vector
        cl_values.SetStrainVector(general_variables.StrainVector);

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& r_alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if (mFinalizedStep) {
                this->GetHistoricalVariables(general_variables,point_number);
            }

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(general_variables, cl_values, point_number);

            // Compute Green-Lagrange Strain
            if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
                mConstitutiveLawVector[point_number]->CalculateMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_PK2);
            } else if (rVariable == ALMANSI_STRAIN_VECTOR) {
                mConstitutiveLawVector[point_number]->CalculateMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_Cauchy);
            } else if (rVariable == HENCKY_STRAIN_VECTOR) {
                mConstitutiveLawVector[point_number]->CalculateValue(cl_values, HENCKY_STRAIN_VECTOR, general_variables.StrainVector);
            }

            if (rOutput[point_number].size() != general_variables.StrainVector.size()) {
                rOutput[point_number].resize( general_variables.StrainVector.size(), false );
            }

            rOutput[point_number] = general_variables.StrainVector;
        }
    } else if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    if ( rOutput.size() != 6 ) {
        std::vector<Vector> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; ++iii) {
            rOutput[iii] = ZeroVector(rOutput[0].size());

            for (unsigned int i_gp = 0; i_gp < integration_point_number; i_gp++) {
                rOutput[iii] += interpol(i_gp, iii) * rOutput_aux[i_gp];
            }
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

    const IndexType integration_point_number = msGeometryData.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_point_number ) {
        rOutput.resize( integration_point_number );
    }

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
        std::vector<Vector> stress_vector;
        if( rVariable == CAUCHY_STRESS_TENSOR ) {
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
        } else {
            this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
        }

        // Loop integration points
        if ( rOutput.size() != stress_vector.size() ) {
            rOutput.resize( stress_vector.size() );
        }

        for (unsigned int point_number = 0; point_number < rOutput.size(); ++point_number ) {
            if (rOutput[point_number].size2() != 3) {
                rOutput[point_number].resize(3, 3, false);
            }
            rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR || rVariable == HENCKY_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ) {
            CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        } else if ( rVariable == ALMANSI_STRAIN_TENSOR ) {
            CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        } else if ( rVariable == HENCKY_STRAIN_TENSOR ) {
            CalculateOnIntegrationPoints( HENCKY_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        }

        // Loop integration points
        if ( rOutput.size() != StrainVector.size() ) {
            rOutput.resize( StrainVector.size() );
        }

        for (unsigned int point_number = 0; point_number < rOutput.size(); ++point_number ) {
            if (rOutput[point_number].size2() != 3) {
                rOutput[point_number].resize(3, 3, false);
            }

            rOutput[point_number] = MathUtils<double>::StrainVectorToTensor(StrainVector[point_number]);
        }
    } else if ( rVariable == CONSTITUTIVE_MATRIX ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& r_constitutive_laws_options = cl_values.GetOptions();
        r_constitutive_laws_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& r_alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(general_variables, cl_values, point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy(cl_values);

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
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

        double& r_alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

            if(rOutput[point_number].size2() != general_variables.F.size2() ) {
                rOutput[point_number].resize( general_variables.F.size1() , general_variables.F.size2() , false );
            }
            rOutput[point_number] = general_variables.F;
        }
    } else if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    if ( rOutput.size() != 6 ) {
        std::vector<Matrix> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; ++iii) {
            rOutput[iii] = ZeroMatrix(rOutput[0].size1(), rOutput[0].size2());

            for (unsigned int i_gp = 0; i_gp < integration_point_number; i_gp++) {
                rOutput[iii] += interpol(i_gp, iii) * rOutput_aux[i_gp];
            }
        }
    }

    KRATOS_CATCH( "" );
}

//**************************** ON INTEGRATION POINTS ******************************//
/******************************** SET DOUBLE VALUE *********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValuesOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/******************************** SET VECTOR VALUE *********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValuesOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    const std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValuesOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/******************************** SET MATRIX VALUE *********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValuesOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    const std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValuesOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/****************************** SET CONSTITUTIVE VALUE *****************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::SetValuesOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    const std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValuesOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//********************************* CHECK VALUES **********************************//
/***********************************************************************************/
/***********************************************************************************/

int SolidShellElementSprism3D6N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    /* Check the neighbours have been calculated */
    // Neighbour nodes
    KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_NODES)) << "The neighbour nodes are not calculated" << std::endl;
    if (this->Has(NEIGHBOUR_NODES)) {
        const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
        KRATOS_ERROR_IF(r_neighbour_nodes.size() == 0) << "The neighbour nodes calculated are empty" << std::endl;
    }

    const int check = BaseType::Check(rCurrentProcessInfo);

    /* Verify compatibility with the constitutive law */
    ConstitutiveLaw::Features law_features;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(law_features);

    // Check strain measure
    bool correct_strain_measure = false;
    for (unsigned int i = 0; i < law_features.mStrainMeasures.size(); ++i) {
        if (law_features.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient ||
            law_features.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal) {
            correct_strain_measure = true;
        }
    }

    KRATOS_ERROR_IF_NOT(correct_strain_measure) << "Constitutive law is not compatible with the element type SolidShellElementSprism3D6N" << std::endl;

    return check;

    KRATOS_CATCH( "" );
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    // Get geometry
    const auto& r_geometry = GetGeometry();

    /* Create and initialize element variables: */
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags& r_constitutive_law_options = cl_values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    /* Reading integration points */
    const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints(this->GetIntegrationMethod());

    double& r_alpha_eas = this->GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives this_cartesian_derivatives;
    this->CalculateCartesianDerivatives(this_cartesian_derivatives);

    /* Calculate common components (B, C) */
    CommonComponents common_components;
    common_components.clear();
    this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

    // Reading integration points
    const Properties& r_properties = GetProperties();
    const auto& r_N_values = msGeometryData.ShapeFunctionsValues(this->GetIntegrationMethod());
    for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
        const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

        // To take in account previous step writing
        if (mFinalizedStep) {
            this->GetHistoricalVariables(general_variables, point_number);
        }

        // Set general variables to constitutive law parameters
        this->SetGeneralVariables(general_variables, cl_values, point_number);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[point_number]->InitializeMaterialResponse(cl_values, GetStressMeasure());

        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->InitializeSolutionStep( r_properties, r_geometry, row( r_N_values, point_number ), rCurrentProcessInfo);
    }

    mFinalizedStep = false;
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Create and initialize element variables:
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

    // Get constitutive law flags:
    Flags& r_constitutive_law_options = cl_values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

    double& r_alpha_eas = this->GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives this_cartesian_derivatives;
    this->CalculateCartesianDerivatives(this_cartesian_derivatives);

    /* Calculate common components (B, C) */
    CommonComponents common_components;
    common_components.clear();
    this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

    // Reading integration points
    const Properties& r_properties = GetProperties();
    const auto& r_N_values = msGeometryData.ShapeFunctionsValues(this->GetIntegrationMethod());
    for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
        const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components,r_integration_points, point_number, r_alpha_eas, zeta_gauss);

        // Set general variables to constitutive law parameters
        this->SetGeneralVariables(general_variables, cl_values, point_number);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(cl_values, general_variables.StressMeasure);

        // Call the constitutive law to finalize the solution step
        mConstitutiveLawVector[point_number]->FinalizeSolutionStep( r_properties, r_geometry, row( r_N_values, point_number ), rCurrentProcessInfo);

        // Call the element internal variables update
        this->FinalizeStepVariables(general_variables, point_number);
    }

    mFinalizedStep = true;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::InitializeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::FinalizeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Call base class
    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    /* Create and initialize element variables: */
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags& r_constitutive_law_options = cl_values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if (mELementalFlags.IsNot(SolidShellElementSprism3D6N::EXPLICIT_RHS_COMPUTATION) ) { // Implicit calculation of the RHS
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    /* Reading integration points */
    const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints( this->GetIntegrationMethod() );

    /* Getting the alpha parameter of the EAS improvement */
    double& r_alpha_eas = this->GetValue(ALPHA_EAS);

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
    for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
        const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(general_variables.B, common_components, zeta_gauss, r_alpha_eas);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

        // Set general variables to constitutive law parameters
        this->SetGeneralVariables(general_variables, cl_values, point_number);

        // Compute stresses and constitutive parameters
        mConstitutiveLawVector[point_number]->CalculateMaterialResponse(cl_values, general_variables.StressMeasure);

        // Calculating weights for integration on the "reference configuration"
        const double integration_weight = r_integration_points[point_number].Weight() * general_variables.detJ;

        /* Integrate in Zeta */
        // EAS components
        IntegrateEASInZeta(general_variables, EAS, zeta_gauss, integration_weight);
    }

    /* Getting the increase of displacements */
    BoundedMatrix<double, 36, 1> delta_disp, current_disp, previous_disp;
    GetVectorCurrentPosition(current_disp);
    GetVectorPreviousPosition(previous_disp);

    // Calculates the increase of displacements
    noalias(delta_disp) = current_disp - previous_disp;

    /* Update alpha EAS */
    if (EAS.mStiffAlpha > std::numeric_limits<double>::epsilon()) { // Avoid division by zero
        r_alpha_eas -= prod(EAS.mHEAS, delta_disp)(0, 0) / EAS.mStiffAlpha;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        mFinalizedStep = true; // the creation is out of the time step, it must be true

        // Getting properties
        const Properties& r_properties = GetProperties();

        // Checking integration order
        if( r_properties.Has(INTEGRATION_ORDER) ) {
            mIntegrationOrder = r_properties.GetValue(INTEGRATION_ORDER);
            if (mIntegrationOrder < 0 || mIntegrationOrder > 5) {
                KRATOS_WARNING("SolidShellElementSprism3D6N") << "Integration order " << mIntegrationOrder << " is not available, using default integration order for the geometry" << std::endl;
                mIntegrationOrder = 0;
            }
        }
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints(this->GetIntegrationMethod());

        /* Constitutive Law initialisation */
        if ( mConstitutiveLawVector.size() != r_integration_points.size() ) {
            mConstitutiveLawVector.resize( r_integration_points.size() );
        }

        /* Implicit or explicit EAS update */
        if(r_properties.Has(CONSIDER_IMPLICIT_EAS_SPRISM_ELEMENT)) {
            mELementalFlags.Set(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT, r_properties.GetValue(CONSIDER_IMPLICIT_EAS_SPRISM_ELEMENT));
        } else {
            mELementalFlags.Set(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT, true);
        }

        /* Total or updated lagrangian */
        if(r_properties.Has(CONSIDER_TOTAL_LAGRANGIAN_SPRISM_ELEMENT)) {
            mELementalFlags.Set(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN, r_properties.GetValue(CONSIDER_TOTAL_LAGRANGIAN_SPRISM_ELEMENT));
        } else {
            mELementalFlags.Set(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN, true);
        }

        /* Quadratic or linear element */
        if(r_properties.Has(CONSIDER_QUADRATIC_SPRISM_ELEMENT)) {
            mELementalFlags.Set(SolidShellElementSprism3D6N::QUADRATIC_ELEMENT, r_properties.GetValue(CONSIDER_QUADRATIC_SPRISM_ELEMENT));
        } else {
            mELementalFlags.Set(SolidShellElementSprism3D6N::QUADRATIC_ELEMENT, true);
        }

        /* Explicit RHS computation */
        if(r_properties.Has(PURE_EXPLICIT_RHS_COMPUTATION)) {
            mELementalFlags.Set(SolidShellElementSprism3D6N::EXPLICIT_RHS_COMPUTATION, r_properties.GetValue(PURE_EXPLICIT_RHS_COMPUTATION));
        } else {
            mELementalFlags.Set(SolidShellElementSprism3D6N::EXPLICIT_RHS_COMPUTATION, false);
        }

        // Resizing the containers
        mAuxContainer.resize( r_integration_points.size() );

        if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) { // Jacobian inverses
            // Compute jacobian inverses and set the domain initial size
            const auto& r_geometry = this->GetGeometry();

            /* Calculating the inverse J0 */
            Matrix J0(3, 3);
            for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
                // Calculating and storing inverse of the jacobian and the parameters needed
                double aux_detJ;
                r_geometry.Jacobian(J0, r_integration_points[point_number]);
                MathUtils<double>::InvertMatrix( J0, mAuxContainer[point_number], aux_detJ );
            }
        } else { // Historic deformation gradient
            for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
                mAuxContainer[point_number] = IdentityMatrix(3);
            }
        }

        /* Initialize AlphaEAS */
        this->SetValue(ALPHA_EAS, 0.0);

        /* Material initialisation */
        InitializeMaterial();
    }

    KRATOS_CATCH("");
}

/***********************************PROTECTED***************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateElementalSystem(
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Get geometries
    const auto& r_geometry = this->GetGeometry();

    /* Create and initialize element variables: */
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags& r_constitutive_law_options = cl_values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( rLocalSystem.CalculationFlags.IsNot(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX) &&
         mELementalFlags.Is(SolidShellElementSprism3D6N::EXPLICIT_RHS_COMPUTATION) ) { // Explicit calculation of the RHS and calculation of the matrix is required
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    } else {
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    }

    /* Reading integration points */
    const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints(this->GetIntegrationMethod());

    /* Getting the alpha parameter of the EAS improvement */
    double& r_alpha_eas = this->GetValue(ALPHA_EAS);

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
    for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
        const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(general_variables.B, common_components, zeta_gauss, r_alpha_eas);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components, r_integration_points, point_number, r_alpha_eas, zeta_gauss);

        // Set general variables to constitutive law parameters
        this->SetGeneralVariables(general_variables, cl_values, point_number);

        // Compute stresses and constitutive parameters
        mConstitutiveLawVector[point_number]->CalculateMaterialResponse(cl_values, general_variables.StressMeasure);

        // Calculating weights for integration on the "reference configuration"
        const double integration_weight = r_integration_points[point_number].Weight() * general_variables.detJ;

        /* Integrate in Zeta */
        // Stresses
        IntegrateStressesInZeta(general_variables, rIntegratedStress, r_alpha_eas, zeta_gauss, integration_weight);
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
        this->CalculateAndAddRHS ( rLocalSystem, general_variables, volume_force, rIntegratedStress, common_components, EAS, r_alpha_eas );
    }

    if ( rLocalSystem.CalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX) ) { // Calculation of the matrix is required
        /* Contribution to the tangent stiffness matrix */
        this->CalculateAndAddLHS( rLocalSystem, general_variables, cl_values, rIntegratedStress, common_components, this_cartesian_derivatives, EAS, r_alpha_eas );
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

    // Get geometries and neighbour nodes
    const auto& r_geometry = this->GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const IndexType number_of_neighbours = NumberOfActiveNeighbours(r_neighbour_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i ) {
        const array_1d<double, 3>& r_current_position  = r_geometry[i].Coordinates();
        const array_1d<double, 3>& r_current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3>& r_previous_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, 1);
        const array_1d<double, 3> previous_position  = r_current_position - (r_current_displacement - r_previous_displacement);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous  Position  node[" << r_geometry[i].Id() << "]: " << previous_position << std::endl;
    }

    for (unsigned int i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3>& r_current_position  = r_neighbour_nodes[i].Coordinates();
        const array_1d<double, 3>& r_current_displacement  = r_neighbour_nodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3>& r_previous_displacement = r_neighbour_nodes[i].FastGetSolutionStepValue(DISPLACEMENT, 1);
        const array_1d<double, 3> previous_position  = r_current_position - (r_current_displacement - r_previous_displacement);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous  Position  neighbour node[" << r_neighbour_nodes[i].Id() << "]: "<<previous_position << std::endl;
    }

    for (unsigned int i = 0; i < number_of_nodes; ++i ) {
        const array_1d<double, 3>& r_current_position  = r_geometry[i].Coordinates();
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Position  node[" << r_geometry[i].Id()<<"]: " << r_current_position << std::endl;
    }

    for (unsigned int i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3>& r_current_position  = r_neighbour_nodes[i].Coordinates();
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Position neighbour node[" << r_neighbour_nodes[i].Id() <<"]: " << r_current_position << std::endl;
    }

    for (unsigned int i = 0; i < number_of_nodes; ++i ) {
        const array_1d<double, 3>& r_previous_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous Displacement node[" << r_geometry[i].Id() << "]: " << r_previous_displacement << std::endl;
    }

    for (unsigned int i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3>& r_previous_displacement = r_neighbour_nodes[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Previous Displacement neighbour node[" << r_neighbour_nodes[i].Id() << "]: " << r_previous_displacement << std::endl;
    }

    for (unsigned int i = 0; i < number_of_nodes; ++i ) {
        const array_1d<double, 3>& r_current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Displacement  node[" << r_geometry[i].Id() << "]: " << r_current_displacement << std::endl;
    }

    for (unsigned int i = 0; i < number_of_neighbours; ++i ) {
        const array_1d<double, 3>& r_current_displacement  = r_neighbour_nodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        KRATOS_INFO("SolidShellElementSprism3D6N") << " Current  Displacement  node[" << r_neighbour_nodes[i].Id() << "]: " << r_current_displacement << std::endl;
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
    const Node& NeighbourNode
    ) const
{
    if (NeighbourNode.Id() == GetGeometry()[Index].Id()) {
        return false;
    } else {
        if (mELementalFlags.Is(SolidShellElementSprism3D6N::QUADRATIC_ELEMENT)) {
            return true;
        } else {
            return false;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SolidShellElementSprism3D6N::NumberOfActiveNeighbours(const WeakPointerVectorNodesType& pNeighbourNodes) const
{
    std::size_t active_neighbours = 0;
    for (unsigned int i = 0; i < pNeighbourNodes.size(); ++i) {
        if (HasNeighbour(i, pNeighbourNodes[i])) {
           ++active_neighbours;
        }
    }
    return active_neighbours;
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetNodalCoordinates(
    BoundedMatrix<double, 12, 3>& rNodesCoordinates,
    const WeakPointerVectorNodesType& rNeighbourNodes,
    const Configuration ThisConfiguration
    ) const
{
    // Get geometry
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();

    // Get the number of neighbours
    rNodesCoordinates = ZeroMatrix(12, 3);
    const unsigned int number_of_neighbours = NumberOfActiveNeighbours(rNeighbourNodes);

     if (ThisConfiguration == Configuration::INITIAL) {
         /* Fill the aux matrix of coordinates */
         for (unsigned int i = 0; i < number_of_nodes; ++i) {
             const array_1d<double, 3>& r_initial_position = r_geometry[i].GetInitialPosition().Coordinates();
             for (unsigned int j = 0; j < 3; ++j) {
                 rNodesCoordinates(i, j) = r_initial_position[j];
             }
         }

         if (number_of_neighbours == number_of_nodes) { // All the possible neighbours
             for (unsigned int i = 0; i < number_of_nodes; ++i) {
                 const array_1d<double, 3>& r_initial_position = rNeighbourNodes[i].GetInitialPosition().Coordinates();
                 for (unsigned int j = 0; j < 3; ++j) {
                    rNodesCoordinates(i + number_of_nodes, j) = r_initial_position[j];
                 }
             }
         } else {
             for (unsigned int i = 0; i < number_of_nodes; ++i) {
                 if (HasNeighbour(i, rNeighbourNodes[i])) {
                     const array_1d<double, 3>& r_initial_position = rNeighbourNodes[i].GetInitialPosition().Coordinates();

                     for (unsigned int j = 0; j < 3; ++j) {
                        rNodesCoordinates(i + number_of_nodes, j) = r_initial_position[j];
                     }
                 } else {
                     for (unsigned int j = 0; j < 3; ++j) {
                        rNodesCoordinates(i + number_of_nodes, j) = 0.0;
                     }
                 }
             }
         }
     } else if (ThisConfiguration == Configuration::CURRENT) {
         /* Fill the aux matrix of coordinates */
         for (unsigned int i = 0; i < number_of_nodes; ++i) {
             const array_1d<double, 3>& r_current_position  = r_geometry[i].Coordinates();
             for (unsigned int j = 0; j < 3; ++j) {
                rNodesCoordinates(i, j) = r_current_position[j];
             }
         }

         if (number_of_neighbours == number_of_nodes) { // All the possible neighours
             for (unsigned int i = 0; i < number_of_nodes; ++i) {
                 const array_1d<double, 3>& r_current_position  = rNeighbourNodes[i].Coordinates();
                 for (unsigned int j = 0; j < 3; ++j) {
                    rNodesCoordinates(i + number_of_nodes, j) = r_current_position[j];
                 }
             }
         } else {
             for (unsigned int i = 0; i < number_of_nodes; ++i) {
                 if (HasNeighbour(i, rNeighbourNodes[i])) {
                     const array_1d<double, 3>& r_current_position  = rNeighbourNodes[i].Coordinates();
                     for (unsigned int j = 0; j < 3; ++j) {
                        rNodesCoordinates(i + number_of_nodes, j) = r_current_position[j];
                     }
                 } else {
                     for (unsigned int j = 0; j < 3; ++j) {
                        rNodesCoordinates(i + number_of_nodes, j) = 0.0;
                     }
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
    BoundedMatrix<double, 12, 3> nodes_coordinates; // Coordinates of the nodes
    const WeakPointerVectorNodesType& neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        this->GetNodalCoordinates(nodes_coordinates, neighbour_nodes, Configuration::INITIAL);
    } else {
        this->GetNodalCoordinates(nodes_coordinates, neighbour_nodes, Configuration::CURRENT);
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
    CalculateCartesianDerOnCenterTrans(rCartesianDerivatives, nodes_coordinates, this_orthogonal_base, GeometricLevel::CENTER); // Center
    CalculateCartesianDerOnCenterTrans(rCartesianDerivatives, nodes_coordinates, this_orthogonal_base, GeometricLevel::LOWER);  // Lower part
    CalculateCartesianDerOnCenterTrans(rCartesianDerivatives, nodes_coordinates, this_orthogonal_base, GeometricLevel::UPPER);  // Upper part

    //******************************** GAUSS POINTS *******************************

    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 0.5;
    local_coordinates[1] = 0.5;
    local_coordinates[2] = -1.0;

    /* Transversal derivatives */
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[0], nodes_coordinates, this_orthogonal_base, local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[3], nodes_coordinates, this_orthogonal_base, local_coordinates);
    local_coordinates[0] = 0.0;
    local_coordinates[2] = -1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[1], nodes_coordinates, this_orthogonal_base, local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[4], nodes_coordinates, this_orthogonal_base, local_coordinates);
    local_coordinates[0] = 0.5;
    local_coordinates[1] = 0.0;
    local_coordinates[2] = -1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[2], nodes_coordinates, this_orthogonal_base, local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(rCartesianDerivatives.TransversalCartesianDerivativesGauss[5], nodes_coordinates, this_orthogonal_base, local_coordinates);

    /* In-plane derivative */
    for (unsigned int i = 0; i < 3 ;++i) {
        if (HasNeighbour(i, neighbour_nodes[i])) { // Assuming that if the upper element has neighbours the lower has too
            CalculateCartesianDerOnGaussPlane(rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i    ], nodes_coordinates, this_orthogonal_base, i, GeometricLevel::LOWER);
            CalculateCartesianDerOnGaussPlane(rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i + 3], nodes_coordinates, this_orthogonal_base, i, GeometricLevel::UPPER);
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

    BoundedMatrix<double, 12, 3 > nodes_coordinates; // Coordinates of the nodes
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    this->GetNodalCoordinates(nodes_coordinates, r_neighbour_nodes, Configuration::CURRENT);

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
    for (unsigned int i = 0; i < 3; ++i) {
        CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i], nodes_coordinates, i, GeometricLevel::LOWER);
        CalculateAndAddBMembrane(rCommonComponents.BMembraneLower, rCommonComponents.CMembraneLower, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i], in_plane_gradient_F_gauss, i);
    }

    rCommonComponents.BMembraneLower *= 1.0/3.0;
    rCommonComponents.CMembraneLower *= 1.0/3.0;

    // Upper face
    for (unsigned int i = 0; i < 3; ++i) {
        CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[i + 3], nodes_coordinates, i, GeometricLevel::UPPER);
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
    CalculateTransverseGradientFinP(transverse_gradient_isoparametric, nodes_coordinates, GeometricLevel::LOWER);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(transverse_gradient.F0, rCartesianDerivatives.TransversalCartesianDerivativesGauss[0], nodes_coordinates);
    CalculateTransverseGradientF(transverse_gradient.F1, rCartesianDerivatives.TransversalCartesianDerivativesGauss[1], nodes_coordinates);
    CalculateTransverseGradientF(transverse_gradient.F2, rCartesianDerivatives.TransversalCartesianDerivativesGauss[2], nodes_coordinates);

    /* Shear contribution to the deformation matrix */
    CalculateAndAddBShear(rCommonComponents.BShearLower, rCommonComponents.CShearLower, rCartesianDerivatives, transverse_gradient, transverse_gradient_isoparametric, GeometricLevel::LOWER);

    // Upper face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(transverse_gradient_isoparametric, nodes_coordinates, GeometricLevel::UPPER);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(transverse_gradient.F0, rCartesianDerivatives.TransversalCartesianDerivativesGauss[3], nodes_coordinates);
    CalculateTransverseGradientF(transverse_gradient.F1, rCartesianDerivatives.TransversalCartesianDerivativesGauss[4], nodes_coordinates);
    CalculateTransverseGradientF(transverse_gradient.F2, rCartesianDerivatives.TransversalCartesianDerivativesGauss[5], nodes_coordinates);

    /* Shear contribution to the deformation matrix */
    CalculateAndAddBShear(rCommonComponents.BShearUpper, rCommonComponents.CShearUpper, rCartesianDerivatives, transverse_gradient, transverse_gradient_isoparametric, GeometricLevel::UPPER);

    /* NORMAL TRANSVERSE */
    /* Calculate f normal components */
    array_1d<double, 3 > F3;
    CalculateTransverseGradientF(F3, rCartesianDerivatives.TransversalCartesianDerivativesCenter, nodes_coordinates);

    /* Calculating the normal transverse strain-displacement matrix */
    CalculateAndAddBNormal(rCommonComponents.BNormal, rCommonComponents.CNormal, rCartesianDerivatives.TransversalCartesianDerivativesCenter, F3);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateLocalCoordinateSystem(
    OrthogonalBase& rOrthogonalBase,
    const OrthogonalBaseApproach ThisOrthogonalBaseApproach,
    const double ThisAngle
    )
{
    KRATOS_TRY;

    // Get the geometry of the element
    const auto& r_geometry = GetGeometry();

    /* Mid-surface vectors */
    double norm; // TODO: Use the geometry normal when available
    array_1d<double, 3 > vxe, vye;
    if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        vxe[0] = 0.5 * ((r_geometry[2].X0() + r_geometry[5].X0()) - (r_geometry[1].X0() + r_geometry[4].X0()));
        vxe[1] = 0.5 * ((r_geometry[2].Y0() + r_geometry[5].Y0()) - (r_geometry[1].Y0() + r_geometry[4].Y0()));
        vxe[2] = 0.5 * ((r_geometry[2].Z0() + r_geometry[5].Z0()) - (r_geometry[1].Z0() + r_geometry[4].Z0()));

        vye[0] = 0.5 * ((r_geometry[0].X0() + r_geometry[3].X0()) - (r_geometry[2].X0() + r_geometry[5].X0()));
        vye[1] = 0.5 * ((r_geometry[0].Y0() + r_geometry[3].Y0()) - (r_geometry[2].Y0() + r_geometry[5].Y0()));
        vye[2] = 0.5 * ((r_geometry[0].Z0() + r_geometry[3].Z0()) - (r_geometry[2].Z0() + r_geometry[5].Z0()));
    } else {
        vxe[0] = 0.5 * ((r_geometry[2].X() + r_geometry[5].X()) - (r_geometry[1].X() + r_geometry[4].X()));
        vxe[1] = 0.5 * ((r_geometry[2].Y() + r_geometry[5].Y()) - (r_geometry[1].Y() + r_geometry[4].Y()));
        vxe[2] = 0.5 * ((r_geometry[2].Z() + r_geometry[5].Z()) - (r_geometry[1].Z() + r_geometry[4].Z()));

        vye[0] = 0.5 * ((r_geometry[0].X() + r_geometry[3].X()) - (r_geometry[2].X() + r_geometry[5].X()));
        vye[1] = 0.5 * ((r_geometry[0].Y() + r_geometry[3].Y()) - (r_geometry[2].Y() + r_geometry[5].Y()));
        vye[2] = 0.5 * ((r_geometry[0].Z() + r_geometry[3].Z()) - (r_geometry[2].Z() + r_geometry[5].Z()));
    }

    MathUtils<double>::CrossProduct(rOrthogonalBase.Vzeta, vxe, vye);
    norm = norm_2(rOrthogonalBase.Vzeta);
    rOrthogonalBase.Vzeta /= norm;

    const double threshold = std::numeric_limits<double>::epsilon();
    double ortho_comp;

    /* Performing the calculation */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    if (ThisOrthogonalBaseApproach == OrthogonalBaseApproach::X) {
        ortho_comp = rOrthogonalBase.Vzeta[1] * rOrthogonalBase.Vzeta[1] + rOrthogonalBase.Vzeta[2] * rOrthogonalBase.Vzeta[2]; // Component in th Y-Z plane
        if (ortho_comp < threshold) { // If rOrthogonalBase.Vzeta is almost orthogonal to  Y-Z plane
            rOrthogonalBase.Veta[0] = - rOrthogonalBase.Vzeta[2]; // Choose rOrthogonalBase.Vxi orthogonal to global Y direction
            rOrthogonalBase.Veta[1] = 0.0;
            rOrthogonalBase.Veta[2] = rOrthogonalBase.Vzeta[0];

            norm = norm_2(rOrthogonalBase.Vxi);
            rOrthogonalBase.Vxi /= norm;
            MathUtils<double>::CrossProduct(rOrthogonalBase.Vxi, rOrthogonalBase.Veta, rOrthogonalBase.Vzeta);
        } else { // SELECT local y=rOrthogonalBase.Vxi in the global YZ plane
            rOrthogonalBase.Vxi[0] = 0.0;
            rOrthogonalBase.Vxi[1] = rOrthogonalBase.Vzeta[2];
            rOrthogonalBase.Vxi[2] = - rOrthogonalBase.Vzeta[1];

            norm = norm_2(rOrthogonalBase.Vxi);
            rOrthogonalBase.Vxi /= norm;

            rOrthogonalBase.Veta[0] = ortho_comp; // Choose rOrthogonalBase.Vxi orthogonal to global X direction
            rOrthogonalBase.Veta[1] = - rOrthogonalBase.Vzeta[0] * rOrthogonalBase.Vzeta[1];
            rOrthogonalBase.Veta[2] = - rOrthogonalBase.Vzeta[0] * rOrthogonalBase.Vzeta[2];

            norm = norm_2(rOrthogonalBase.Veta);
            rOrthogonalBase.Veta /= norm;
        }
    } else if (ThisOrthogonalBaseApproach == OrthogonalBaseApproach::Y) {
        ortho_comp = rOrthogonalBase.Vzeta[0] * rOrthogonalBase.Vzeta[0] + rOrthogonalBase.Vzeta[2] * rOrthogonalBase.Vzeta[2]; // Component in th Z-X plane
        if (ortho_comp < threshold) { // If vze is almost orthogonal to  Z-X plane
            rOrthogonalBase.Veta[0] =       0.0; // Choose rOrthogonalBase.Vxi orthogonal to global X direction
            rOrthogonalBase.Veta[1] =   rOrthogonalBase.Vzeta[2];
            rOrthogonalBase.Veta[2] = - rOrthogonalBase.Vzeta[1];

            norm = norm_2(rOrthogonalBase.Veta);
            rOrthogonalBase.Veta /= norm;
            MathUtils<double>::CrossProduct(rOrthogonalBase.Vxi, rOrthogonalBase.Veta, rOrthogonalBase.Vzeta);
        } else { // SELECT local z=rOrthogonalBase.Vxi in the global ZX plane
            rOrthogonalBase.Vxi[0] = - rOrthogonalBase.Vzeta[2]; // Choose rOrthogonalBase.Vxi orthogonal to global Y direction
            rOrthogonalBase.Vxi[1] = 0.0;
            rOrthogonalBase.Vxi[2] = - rOrthogonalBase.Vzeta[0];

            norm = norm_2(rOrthogonalBase.Vxi);
            rOrthogonalBase.Vxi /= norm;

            rOrthogonalBase.Veta[0] = - rOrthogonalBase.Vzeta[0] * rOrthogonalBase.Vzeta[1];
            rOrthogonalBase.Veta[1] = ortho_comp;
            rOrthogonalBase.Veta[2] = - rOrthogonalBase.Vzeta[2] * rOrthogonalBase.Vzeta[1];

            norm = norm_2(rOrthogonalBase.Veta);
            rOrthogonalBase.Veta /= norm;
        }
    } else if (ThisOrthogonalBaseApproach == OrthogonalBaseApproach::Z) {
        ortho_comp = rOrthogonalBase.Vzeta[0] * rOrthogonalBase.Vzeta[0] + rOrthogonalBase.Vzeta[1] * rOrthogonalBase.Vzeta[1]; // Component in th X-Y plane
        if (ortho_comp < threshold) { // If vze is almost orthogonal to  X-Y plane
            rOrthogonalBase.Veta[0] = 0.0; // Choose rOrthogonalBase.Vxi orthogonal to global X direction
            rOrthogonalBase.Veta[1] = rOrthogonalBase.Vzeta[2];
            rOrthogonalBase.Veta[2] = - rOrthogonalBase.Vzeta[1];

            norm = norm_2(rOrthogonalBase.Veta);
            rOrthogonalBase.Veta /= norm;
            MathUtils<double>::CrossProduct(rOrthogonalBase.Vxi, rOrthogonalBase.Veta, rOrthogonalBase.Vzeta);
        } else { // SELECT local x=rOrthogonalBase.Vxi in the global XY plane
            rOrthogonalBase.Vxi[0] = - rOrthogonalBase.Vzeta[1];
            rOrthogonalBase.Vxi[1] = rOrthogonalBase.Vzeta[0];
            rOrthogonalBase.Vxi[2] = 0.0;

            norm = norm_2(rOrthogonalBase.Vxi);
            rOrthogonalBase.Vxi /= norm;

            rOrthogonalBase.Veta[0] = - rOrthogonalBase.Vzeta[0] * rOrthogonalBase.Vzeta[2]; // Choose rOrthogonalBase.Vxi orthogonal to global Z direction
            rOrthogonalBase.Veta[1] = - rOrthogonalBase.Vzeta[1] * rOrthogonalBase.Vzeta[2];
            rOrthogonalBase.Veta[2] = ortho_comp;

            norm = norm_2(rOrthogonalBase.Veta);
            rOrthogonalBase.Veta /= norm;
        }
    } else {
        rOrthogonalBase.Vxi[0] = 1.0;
        rOrthogonalBase.Vxi[1] = 0.0;
        rOrthogonalBase.Vxi[2] = 0.0;

        rOrthogonalBase.Veta[0] = 0.0;
        rOrthogonalBase.Veta[1] = 1.0;
        rOrthogonalBase.Veta[2] = 0.0;
    }

    if (ThisAngle != 0.0) {
        // Compute angle between local system rOrthogonalBase.Vxi-rOrthogonalBase.Veta and L1
        const double cosa = std::cos(ThisAngle);
        const double sina = std::sin(ThisAngle);
        // Rotate local system rOrthogonalBase.Vxi-rOrthogonalBase.Veta to best fit L1-L2
        rOrthogonalBase.Vzeta = rOrthogonalBase.Vxi; // Reusing as auxiliary value
        rOrthogonalBase.Vxi =    cosa * rOrthogonalBase.Vxi    + sina * rOrthogonalBase.Veta;
        rOrthogonalBase.Veta = - sina * rOrthogonalBase.Vzeta  + cosa * rOrthogonalBase.Veta;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateIdVector(array_1d<IndexType, 18>& rIdVector)
{
    KRATOS_TRY;

    // Get the number of nodes
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();

    // Get the list of neighbour nodes
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    /* Compute ID vector */ // TODO: Optimize this
    unsigned int index = 18;
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        if (HasNeighbour(i, r_neighbour_nodes[i])) {
            for (unsigned int j = 0; j < 3; ++j) {
                rIdVector[i * 3 + j] = index;
                ++index;
            }
        } else {
            for (unsigned int j = 0; j < 3; ++j) {
                rIdVector[i * 3 + j] = 36;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::ComputeLocalDerivatives(
    BoundedMatrix<double, 6, 3 >& rLocalDerivativePatch,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    const double L_1 = 0.5 * (1.0 - rLocalCoordinates[2]);
    const double L_2 = 0.5 * (1.0 + rLocalCoordinates[2]);
//     const double zzeta = 1.0 - rLocalCoordinates[0] - rLocalCoordinates[1];

    /* Derivative in direction nu and xi */
    // Lower face
    rLocalDerivativePatch(0, 0) = - L_1;
    rLocalDerivativePatch(1, 0) =   L_1;
    rLocalDerivativePatch(2, 0) =   0.0;

    rLocalDerivativePatch(0, 1) = - L_1;
    rLocalDerivativePatch(1, 1) =   0.0;
    rLocalDerivativePatch(2, 1) =   L_1;

    // Upper face
    rLocalDerivativePatch(3, 0) = - L_2;
    rLocalDerivativePatch(4, 0) =   L_2;
    rLocalDerivativePatch(5, 0) =   0.0;

    rLocalDerivativePatch(3, 1) = - L_2;
    rLocalDerivativePatch(4, 1) =   0.0;
    rLocalDerivativePatch(5, 1) =   L_2;

//     /* Derivative in direction zeta */
//     rLocalDerivativePatch(0, 2) = - zzeta/2.0;
//     rLocalDerivativePatch(1, 2) = - rLocalCoordinates[0]/2.0;
//     rLocalDerivativePatch(2, 2) = - rLocalCoordinates[1]/2.0;
//     rLocalDerivativePatch(3, 2) =   zzeta/2.0;
//     rLocalDerivativePatch(4, 2) =   rLocalCoordinates[0]/2.0;
//     rLocalDerivativePatch(5, 2) =   rLocalCoordinates[1]/2.0;

    /* Derivative in direction zeta */
    rLocalDerivativePatch(0, 2) = - 1.0 + rLocalCoordinates[1] + rLocalCoordinates[0];
    rLocalDerivativePatch(1, 2) = - rLocalCoordinates[0];
    rLocalDerivativePatch(2, 2) = - rLocalCoordinates[1];
    rLocalDerivativePatch(3, 2) =   1.0 - rLocalCoordinates[1] - rLocalCoordinates[0];
    rLocalDerivativePatch(4, 2) =   rLocalCoordinates[0];
    rLocalDerivativePatch(5, 2) =   rLocalCoordinates[1];
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::ComputeLocalDerivativesQuadratic(
    BoundedMatrix<double, 4, 2 >& rLocalDerivativePatch,
    const unsigned int NodeGauss
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
    GeometryType::JacobiansType& rJ,
    std::vector< Matrix >& rJinv,
    Vector& rDetJ,
    const unsigned int PointNumber,
    const double ZetaGauss
    )
{
    /* Fill the aux matrix of coordinates */
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    BoundedMatrix<double, 3, 6 > nodes_coordinates;
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3>& r_current_position  = r_geometry[i].Coordinates();
        for (unsigned int j = 0; j < 3; ++j) {
            nodes_coordinates(j, i) = r_current_position[j];
        }
    }

    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 1.0/3.0;
    local_coordinates[1] = 1.0/3.0;
    local_coordinates[2] = ZetaGauss;

    /* Local derivatives patch */
    BoundedMatrix<double, 6, 3 > rLocalDerivativePatch;
    ComputeLocalDerivatives(rLocalDerivativePatch, local_coordinates);

    /* Compute Jacobian */
    noalias(rJ[PointNumber]) = prod(nodes_coordinates, rLocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    MathUtils<double>::InvertMatrix(rJ[PointNumber], rJinv[PointNumber], rDetJ[PointNumber]);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateJacobian(
    double& rDetJ,
    BoundedMatrix<double, 3, 3>& rJ,
    BoundedMatrix<double, 6, 3>& rLocalDerivativePatch,
    const BoundedMatrix<double, 12, 3>& rNodesCoordinates,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    // Get the number of nodes
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();

    /* Auxiliary coordinates of the nodes */
    BoundedMatrix<double, 3, 6> nodes_coord_aux;

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            nodes_coord_aux(j, i) = rNodesCoordinates(i, j);
        }
    }

    /* Local derivatives patch */
    ComputeLocalDerivatives(rLocalDerivativePatch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(rJ) = prod(nodes_coord_aux, rLocalDerivativePatch);

    /* Compute determinant */
    rDetJ = MathUtils<double>::Det3(rJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateJacobianAndInv(
    BoundedMatrix<double, 3, 3>& rJ,
    BoundedMatrix<double, 3, 3>& rJinv,
    BoundedMatrix<double, 6, 3>& rLocalDerivativePatch,
    const BoundedMatrix<double, 3, 6>& rNodesCoordinates,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Local derivatives patch */
    ComputeLocalDerivatives(rLocalDerivativePatch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(rJ) = prod(rNodesCoordinates, rLocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    MathUtils<double>::InvertMatrix(rJ, rJinv, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateJacobianAndInv(
    BoundedMatrix<double, 3, 3>& rJ,
    BoundedMatrix<double, 3, 3>& rJinv,
    const BoundedMatrix<double, 3, 6>& rNodesCoordinates,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Local derivatives patch */
    BoundedMatrix<double, 6, 3> local_derivative_patch;
    ComputeLocalDerivatives(local_derivative_patch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(rJ) = prod(rNodesCoordinates, local_derivative_patch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    MathUtils<double>::InvertMatrix(rJ, rJinv, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerivativesOnCenterPlane(
    BoundedMatrix<double, 2, 4 >& rCartesianDerivativesCenter,
    const OrthogonalBase& rOrthogonalBase,
    const GeometricLevel Part
    )
{
    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Get the index of the nodes
    // 0- Lower face
    // 1- Upper face
    // 2- Lower face
    // 3- Upper face
    // 4- Lower face
    // 5- Upper face
    // 6- Lower face
    const unsigned int index = Part == GeometricLevel::UPPER ? 3 : 0;

    double norm0, norm;
    array_1d<double, 3 > vxe, vye;
    if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        vxe[0] = r_geometry[2 + index].X0() - r_geometry[1 + index].X0();
        vxe[1] = r_geometry[2 + index].Y0() - r_geometry[1 + index].Y0();
        vxe[2] = r_geometry[2 + index].Z0() - r_geometry[1 + index].Z0();

        vye[0] = r_geometry[0 + index].X0() - r_geometry[2 + index].X0();
        vye[1] = r_geometry[0 + index].Y0() - r_geometry[2 + index].Y0();
        vye[2] = r_geometry[0 + index].Z0() - r_geometry[2 + index].Z0();
    } else {
        vxe[0] = r_geometry[2 + index].X() - r_geometry[1 + index].X();
        vxe[1] = r_geometry[2 + index].Y() - r_geometry[1 + index].Y();
        vxe[2] = r_geometry[2 + index].Z() - r_geometry[1 + index].Z();

        vye[0] = r_geometry[0 + index].X() - r_geometry[2 + index].X();
        vye[1] = r_geometry[0 + index].Y() - r_geometry[2 + index].Y();
        vye[2] = r_geometry[0 + index].Z() - r_geometry[2 + index].Z();
    }

    array_1d<double, 3 > t1g, t2g, t3g;
    MathUtils<double>::CrossProduct(t3g, vxe, vye);
    norm0 = norm_2(t3g);
    t3g /= norm0;

    MathUtils<double>::CrossProduct(t2g, t3g, rOrthogonalBase.Vxi);
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

    rCartesianDerivativesCenter = ZeroMatrix(2, 4);
    for (unsigned int i = 0; i < 3; ++i) {
       rCartesianDerivativesCenter(0, i) = - b[i];
       rCartesianDerivativesCenter(1, i) =   a[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerOnGaussPlane(
    BoundedMatrix<double, 2, 4 >& rInPlaneCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3>& rNodesCoordinates,
    const OrthogonalBase& rOrthogonalBase,
    const unsigned int NodeGauss,
    const GeometricLevel Part
    )
{
    const unsigned int index = Part == GeometricLevel::UPPER ? 3 : 0;

    /* Local derivatives patch */
    BoundedMatrix<double, 4, 2 > local_derivative_patch;
    ComputeLocalDerivativesQuadratic(local_derivative_patch,NodeGauss);

    /* Auxiliary coordinates of the nodes */
    BoundedMatrix<double, 3, 4 > nodes_coord_aux;

    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            nodes_coord_aux(j, i) = rNodesCoordinates(i + index, j);
        }
    }

    for (unsigned int j = 0; j < 3; ++j) {
        nodes_coord_aux(j, 3) = rNodesCoordinates(NodeGauss + 6 + index, j);
    }

    /* Compute local derivatives */
    const BoundedMatrix<double, 3, 2> Xd = prod(nodes_coord_aux, local_derivative_patch);

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
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t1g, t2g, t3g, rOrthogonalBase.Vxi, Xdxi, Xdeta);

    /* Compute Jacobian */
    BoundedMatrix<double, 2, 2 > jac;
    jac(0, 0) = inner_prod(Xdxi,  t1g);
    jac(0, 1) = inner_prod(Xdxi,  t2g);
    jac(1, 0) = inner_prod(Xdeta, t1g);
    jac(1, 1) = inner_prod(Xdeta, t2g);

    /* Compute the inverse of the Jacobian */
    double aux_det;
    BoundedMatrix<double, 2, 2 > JinvPlane;
    MathUtils<double>::InvertMatrix(jac, JinvPlane, aux_det);

    /* Compute the Cartesian derivatives */
    noalias(rInPlaneCartesianDerivativesGauss) = prod(JinvPlane, trans(local_derivative_patch));
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerOnGaussTrans(
    BoundedMatrix<double, 6, 1>& rTransversalCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3>& rNodesCoordinates,
    const OrthogonalBase& rOrthogonalBase,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Compute local derivatives */
    double det;
    BoundedMatrix<double, 3, 3 > Xd;
    BoundedMatrix<double, 6, 3 > local_derivatives_patch;
    CalculateJacobian(det, Xd, local_derivatives_patch, rNodesCoordinates, rLocalCoordinates);

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
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t, rOrthogonalBase.Vxi, Xdxi, Xdeta);

    /* Compute Jacobian */
    const BoundedMatrix<double, 3, 3 > jac = prod(t, Xd);

    /* Compute inverse of the Jaccobian (just third column) */
    BoundedMatrix<double, 3 ,1> JinvTrans;
    JinvTrans(0, 0) =   (jac(0, 1) * jac(1, 2) - jac(0, 2) * jac(1, 1)) / det;
    JinvTrans(1, 0) = - (jac(0, 0) * jac(1, 2) - jac(0, 2) * jac(1, 0)) / det;
    JinvTrans(2, 0) =   (jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1)) / det;

    /* Compute Cartesian derivatives */
    noalias(rTransversalCartesianDerivativesGauss) = prod(local_derivatives_patch, JinvTrans);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateCartesianDerOnCenterTrans(
    CartesianDerivatives& rCartesianDerivatives,
    const BoundedMatrix<double, 12, 3>& rNodesCoordinates,
    const OrthogonalBase& rOrthogonalBase,
    const GeometricLevel Part
    )
{
    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 1.0/3.0;
    local_coordinates[1] = 1.0/3.0;

    if (Part == GeometricLevel::CENTER) {
        local_coordinates[2] =   0.0;
    } else if (Part == GeometricLevel::LOWER) {
        local_coordinates[2] = - 1.0;
    } else if (Part == GeometricLevel::UPPER) {
        local_coordinates[2] =   1.0;
    }

    /* Auxiliary coordinates of the nodes */
    BoundedMatrix<double, 3, 6 > nodes_coord_aux;
    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            nodes_coord_aux(j, i) = rNodesCoordinates(i, j);
        }
    }

    /* Auxiliary components to calculate the Jacobian and his inverse */
    BoundedMatrix<double, 3, 3 > J, Jinv;

    if (Part == GeometricLevel::CENTER) {
        /* Calculate the Jacobian and his inverse */
        BoundedMatrix<double, 6, 3 > local_derivatives_patch;
        CalculateJacobianAndInv(J, Jinv, local_derivatives_patch, nodes_coord_aux, local_coordinates);

        // Compute cartesian (y3) derivatives of the shape functions necessary to compute f_3
        /* Compute Cartesian derivatives */
        const BoundedMatrix<double, 6, 3 > transverse_cartesian_derivatives_gauss_aux = prod(local_derivatives_patch, Jinv);

        for (unsigned int i = 0; i < 6 ; ++i) {
            rCartesianDerivatives.TransversalCartesianDerivativesCenter(i, 0) = inner_prod(rOrthogonalBase.Vzeta, row(transverse_cartesian_derivatives_gauss_aux, i));
        }
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
             rCartesianDerivatives.JInvPlaneLower(0, 0) = inner_prod(Xdxi,  rOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneLower(0, 1) = inner_prod(Xdeta, rOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneLower(1, 0) = inner_prod(Xdxi,  rOrthogonalBase.Veta);
             rCartesianDerivatives.JInvPlaneLower(1, 1) = inner_prod(Xdeta, rOrthogonalBase.Veta);
         } else if (Part == GeometricLevel::UPPER) {
             rCartesianDerivatives.JInvPlaneUpper(0, 0) = inner_prod(Xdxi,  rOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneUpper(0, 1) = inner_prod(Xdeta, rOrthogonalBase.Vxi);
             rCartesianDerivatives.JInvPlaneUpper(1, 0) = inner_prod(Xdxi,  rOrthogonalBase.Veta);
             rCartesianDerivatives.JInvPlaneUpper(1, 1) = inner_prod(Xdeta, rOrthogonalBase.Veta);
         }
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateInPlaneGradientFGauss(
    BoundedMatrix<double, 3, 2>& rInPlaneGradientFGauss,
    const BoundedMatrix<double, 2, 4>& rInPlaneCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3>& rNodesCoordinates,
    const unsigned int NodeGauss,
    const GeometricLevel Part
    )
{
    /* Auxiliary operators */
    const IndexType index = Part == GeometricLevel::UPPER ? 3 : 0;
    BoundedMatrix<double, 3, 3> nodes_coord_aux;
    BoundedMatrix<double, 3, 2> in_plane_cartesian_derivatives_gauss_aux;

    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            nodes_coord_aux(j, i) = rNodesCoordinates(i + index, j);
        }

        in_plane_cartesian_derivatives_gauss_aux(i, 0) = rInPlaneCartesianDerivativesGauss(0, i);
        in_plane_cartesian_derivatives_gauss_aux(i, 1) = rInPlaneCartesianDerivativesGauss(1, i);
    }

    noalias(rInPlaneGradientFGauss) = prod(nodes_coord_aux, in_plane_cartesian_derivatives_gauss_aux);

    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    if (HasNeighbour(NodeGauss, r_neighbour_nodes[NodeGauss])) {
        for (unsigned int j = 0; j < 3 ; ++j) {
            rInPlaneGradientFGauss(j, 0) += rNodesCoordinates(NodeGauss + 6 + index, j) * rInPlaneCartesianDerivativesGauss(0, 3);
            rInPlaneGradientFGauss(j, 1) += rNodesCoordinates(NodeGauss + 6 + index, j) * rInPlaneCartesianDerivativesGauss(1, 3);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateTransverseGradientF(
    array_1d<double, 3>& rTransverseGradientF,
    const BoundedMatrix<double, 6, 1>& rTransversalCartesianDerivativesGauss,
    const BoundedMatrix<double, 12, 3>& rNodesCoordinates
    )
{
    for (unsigned int j = 0; j < 3; ++j) {
        rTransverseGradientF[j] = 0.0;
    }

    // Get number of nodes
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            rTransverseGradientF[j] += rTransversalCartesianDerivativesGauss(i, 0) * rNodesCoordinates(i, j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateTransverseGradientFinP(
    TransverseGradientIsoParametric& rTransverseGradientIsoParametric,
    const BoundedMatrix<double, 12, 3>& rNodesCoordinates,
    const GeometricLevel Part
    )
{
    const unsigned int index = Part == GeometricLevel::UPPER ? 3 : 0;

    for (unsigned int i = 0; i < 3; ++i) {
        rTransverseGradientIsoParametric.Ft[i]   = rNodesCoordinates(2 + index, i) - rNodesCoordinates(1 + index, i);
        rTransverseGradientIsoParametric.Fxi[i]  = rNodesCoordinates(0 + index, i) - rNodesCoordinates(2 + index, i);
        rTransverseGradientIsoParametric.Feta[i] = rNodesCoordinates(1 + index, i) - rNodesCoordinates(0 + index, i);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddBMembrane(
    BoundedMatrix<double, 3, 18>& rBMembrane,
    BoundedMatrix<double, 3, 1>& rCMembrane,
    const BoundedMatrix<double, 2, 4>& rInPlaneCartesianDerivativesGauss,
    const BoundedMatrix<double, 3, 2>& rInPlaneGradientFGauss,
    const unsigned int NodeGauss
    )
{
    for (unsigned int i = 0; i < 4; ++i) {
        unsigned int base = i * 3;
        if (i == 3) {
            base += NodeGauss * 3;
        }
        for (unsigned int j = 0; j < 3; ++j) {
            rBMembrane(0, base + j) += rInPlaneCartesianDerivativesGauss(0, i) * rInPlaneGradientFGauss(j, 0);
            rBMembrane(1, base + j) += rInPlaneCartesianDerivativesGauss(1, i) * rInPlaneGradientFGauss(j, 1);
            rBMembrane(2, base + j) += rInPlaneCartesianDerivativesGauss(1, i) * rInPlaneGradientFGauss(j, 0)
                                    +  rInPlaneCartesianDerivativesGauss(0, i) * rInPlaneGradientFGauss(j, 1);
        }
    }

    /* Calculate de componets of Cauchy tensor */
    // In plane auxiliary components
    array_1d<double, 3 > aux_deformation_gradient_F1, aux_deformation_gradient_F2;

    for (unsigned int i = 0; i < 3; ++i) {
        aux_deformation_gradient_F1[i] = rInPlaneGradientFGauss(i, 0);
        aux_deformation_gradient_F2[i] = rInPlaneGradientFGauss(i, 1);
    }

    rCMembrane(0, 0) += inner_prod(aux_deformation_gradient_F1, aux_deformation_gradient_F1);
    rCMembrane(1, 0) += inner_prod(aux_deformation_gradient_F2, aux_deformation_gradient_F2);
    rCMembrane(2, 0) += inner_prod(aux_deformation_gradient_F1, aux_deformation_gradient_F2);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddMembraneKgeometric(
    BoundedMatrix<double, 36, 36>& rKgeometricmembrane,
    const CartesianDerivatives& rCartesianDerivatives,
    const array_1d<double, 3>& rSMembrane,
    const GeometricLevel Part
    )
{
    const unsigned int index = static_cast<unsigned int>(Part);
    const unsigned int auxiliary_index = Part == GeometricLevel::UPPER ? 3 : 0;

    BoundedMatrix<double, 6, 6> H = ZeroMatrix(6, 6);

    unsigned int ii, jj;
    for (unsigned int i = 0; i < 4; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            // Gauss 1
            ii = i;
            jj = j;
            H(ii, jj) += rSMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](0, j)
                       + rSMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](1, j)
                       + rSMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](1, j)
                                        + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 0](0, j));

            // Gauss 2
            ii = (i ==  3) ? 4 : i;
            jj = (j ==  3) ? 4 : j;

            H(ii, jj) += rSMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](0, j)
                       + rSMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](1, j)
                       + rSMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](1, j)
                                        + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 1](0, j));

            // Gauss 3
            ii = (i ==  3) ? 5 : i;
            jj = (j ==  3) ? 5 : j;

            H(ii, jj) += rSMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](0, j)
                       + rSMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](1, j)
                       + rSMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](1, j)
                                        + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliary_index + 2](0, j));
        }
    }

    H *= 1.0/3.0;

    // Assembling in Kgeometricmembrane
    unsigned int rowindex, colindex;
    for (unsigned int i = 0; i < 6; ++i) {
        rowindex = (i < 3) ? i * 3 + index : i * 3 + index + 9;
        for (unsigned int j = i; j < 6; ++j) {
            colindex = (j < 3) ? j * 3 + index : j * 3 + index + 9;
            for(unsigned int ii = 0; ii < 3; ++ii) {
                rKgeometricmembrane(rowindex + ii,colindex + ii) += H (i, j);
                if (rowindex != colindex) { // Skip diagonal
                    rKgeometricmembrane(colindex + ii, rowindex + ii) += H (i, j); // Symmetric part
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddBShear(
    BoundedMatrix<double, 2, 18>& rBShear,
    BoundedMatrix<double, 2, 1>& rCShear,
    const CartesianDerivatives& rCartesianDerivatives,
    const TransverseGradient& rTransverseGradient,
    const TransverseGradientIsoParametric& rTransverseGradientIsoParametric,
    const GeometricLevel Part
    )
{
    const unsigned int index = static_cast<unsigned int>(Part);
    const unsigned int auxiliary_index = Part == GeometricLevel::UPPER ? 3 : 0;

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

    BoundedMatrix<double, 3, 18> aux_b_shear = ZeroMatrix(3, 18);

    /* First contribution*/
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        unsigned int base = i * 3;
        for (unsigned int j = 0; j < 3; ++j) {
            aux_b_shear(0, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliary_index + 0](i, 0) * rTransverseGradientIsoParametric.Ft[j];
            aux_b_shear(1, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliary_index + 1](i, 0) * rTransverseGradientIsoParametric.Fxi[j];
            aux_b_shear(2, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliary_index + 2](i, 0) * rTransverseGradientIsoParametric.Feta[j];
        }
    }

    /* Second contibution */
    for (unsigned int i = 0; i < 3; ++i) {
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
    noalias(rBShear) = prod(aux_prod, aux_b_shear);

    // Calculating the components of C
    BoundedMatrix<double, 3, 1 > aux_c_shear;
    aux_c_shear(0, 0) = inner_prod(rTransverseGradientIsoParametric.Ft  , rTransverseGradient.F0);
    aux_c_shear(1, 0) = inner_prod(rTransverseGradientIsoParametric.Fxi , rTransverseGradient.F1);
    aux_c_shear(2, 0) = inner_prod(rTransverseGradientIsoParametric.Feta, rTransverseGradient.F2);

    noalias(rCShear) = prod(aux_prod, aux_c_shear);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddShearKgeometric(
    BoundedMatrix<double, 36, 36>& rKgeometricshear,
    const CartesianDerivatives& rCartesianDerivatives,
    const array_1d<double, 2>& rSShear,
    const GeometricLevel Part
    )
{
//     const IndexType index = static_cast<IndexType>(Part);
    const IndexType auxiliary_index = Part == GeometricLevel::UPPER ? 3 : 0;

    const BoundedMatrix<double, 2, 2 >& JInvPlane = Part == GeometricLevel::UPPER ? rCartesianDerivatives.JInvPlaneUpper : rCartesianDerivatives.JInvPlaneLower;

    const double Q1 = 1.0/3.0 * (rSShear[0] * JInvPlane(0, 0) + rSShear[1] * JInvPlane(0, 1));
    const double Q2 = 1.0/3.0 * (rSShear[0] * JInvPlane(1, 0) + rSShear[1] * JInvPlane(1, 1));

//    array_1d<double, 3 > q;
//    q[0] = -Q1 + Q2;
//    q[1] = -(Q1 + 2.0 * Q2);
//    q[2] = (2.0 * Q1 + Q2);

//    int delta;
//    if (index == 9)
//        delta = 3;
//    else
//        delta = 0;

//    for (unsigned int i = 0; i < 3; ++i) { // For each DOF
//        /* First assembling */
//        Kgeometricshear(i + index + 3, i + index + 3) -= q[0] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliary_index](1 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) += q[0] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliary_index](2 + delta, 0);

//        /* Second assembling */
//        Kgeometricshear(i + index, i + index)         += q[1] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliary_index](0 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) -= q[1] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliary_index](2 + delta, 0);
//        /* Third assembling */
//        Kgeometricshear(i + index, i + index)         -= q[2] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliary_index](0 + delta, 0);
//        Kgeometricshear(i + index + 3, i + index + 3) += q[2] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliary_index](1 + delta, 0);
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
    for (unsigned int k = 0; k < 3; k++) {
        IndexType l = 0; // Initializes DOF associated to N_3
        for (unsigned int i = 0; i < 6; ++i) { //  For each node
            if (k == 0) {
                value = (-Q1 + Q2) *  rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliary_index](i, 0);
            } else if (k == 1) {
                value = -(Q1 + 2.0 * Q2) * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliary_index](i, 0);
            } else if (k == 2) {
                value = (2.0 * Q1 + Q2) * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliary_index](i, 0);
            }

            for (unsigned int j = 0; j < 3; ++j) { // For each DOF (diagonal only)
                rKgeometricshear(n1[k] + j, l + j) += value;
                rKgeometricshear(l + j, n1[k] + j) += value;
            }

            for (unsigned int j = 0; j < 3; ++j) { // For each DOF (diagonal only)
                rKgeometricshear(n2[k] + j, l + j) -= value;
                rKgeometricshear(l + j, n2[k] + j) -= value;
            }

            l += 3; // Increment DOF position I
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddBNormal(
    BoundedMatrix<double, 1, 18>& rBNormal,
    double& rCNormal,
    const BoundedMatrix<double, 6, 1>& rTransversalCartesianDerivativesCenter,
    const array_1d<double, 3>& rTransversalDeformationGradientF
    )
{
    unsigned int base;
    for (unsigned int i = 0; i < 6; ++i) {
        base = i * 3;
        for (unsigned int j = 0; j < 3; ++j) {
            rBNormal(0, base + j) = rTransversalCartesianDerivativesCenter(i, 0) * rTransversalDeformationGradientF[j];
        }
    }

    rCNormal = inner_prod(rTransversalDeformationGradientF, rTransversalDeformationGradientF);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateAndAddNormalKgeometric(
    BoundedMatrix<double, 36, 36>& rKgeometricnormal,
    const BoundedMatrix<double, 6, 1>& rTransversalCartesianDerivativesCenter,
    const double SNormal
    )
{
    BoundedMatrix<double, 6, 6 > H = ZeroMatrix(6, 6);
    for (unsigned int i = 0; i < 6; ++i) {
        const double aux = SNormal * rTransversalCartesianDerivativesCenter(i, 0);
        for (unsigned int j = 0; j < 6; ++j) {
            H(i, j) =  aux * rTransversalCartesianDerivativesCenter(j, 0);
        }
    }

    noalias(H) = SNormal * prod(rTransversalCartesianDerivativesCenter, trans(rTransversalCartesianDerivativesCenter));

    unsigned int rowindex, colindex;
    for (unsigned int i = 0; i < 6; ++i) {
        rowindex = i * 3;
        for (unsigned int j = 0; j < 6; ++j) {
            colindex = j * 3;
            for (unsigned int ii = 0; ii < 3; ++ii) {
                rKgeometricnormal(rowindex + ii, colindex + ii) += H(i, j);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetVectorCurrentPosition(BoundedMatrix<double, 36, 1>& rVectorCurrentPosition)
{
    KRATOS_TRY;

    // Get neighbour nodes
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    // Get geometry
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();

    /* Element nodes */
    for (unsigned int index = 0; index < number_of_nodes; ++index) {
        const array_1d<double,3>& r_current_position = r_geometry[index].Coordinates();
        for (unsigned int j = 0; j < 3; ++j) {
            rVectorCurrentPosition(index * 3 + j, 0) = r_current_position[j];
        }
    }

    /* Neighbour nodes */
    const SizeType number_of_neighbours = NumberOfActiveNeighbours(r_neighbour_nodes);

    if (number_of_neighbours == 6) { // All the possible neighbours
        for (unsigned int index = 0; index < number_of_nodes; ++index) {
            const array_1d<double,3>& r_current_position = r_neighbour_nodes[index].Coordinates();
            for (unsigned int j = 0; j < 3; ++j) {
                rVectorCurrentPosition(18 + index * 3 + j, 0) = r_current_position[j];
            }
        }
    } else {
        for (unsigned int index = 0; index < number_of_nodes; ++index) {
            if (HasNeighbour(index, r_neighbour_nodes[index])) {
                const array_1d<double,3>& r_current_position = r_neighbour_nodes[index].Coordinates();
                for (unsigned int j = 0; j < 3; ++j) {
                    rVectorCurrentPosition(18 + index * 3 + j, 0) = r_current_position[j];
                }
            } else {
                for (unsigned int j = 0; j < 3; ++j) {
                    rVectorCurrentPosition(18 + index * 3 + j, 0) = 0.0;
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetVectorPreviousPosition(BoundedMatrix<double, 36, 1>& rVectorPreviousPosition)
{
    KRATOS_TRY;

    // Get neighbour nodes
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    // Get geometry
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();

    /* Element nodes */
    for (unsigned int index = 0; index < number_of_nodes; ++index) {
        const array_1d<double,3>& r_previous_position = r_geometry[index].GetInitialPosition().Coordinates()
                                                      + r_geometry[index].FastGetSolutionStepValue(DISPLACEMENT, 1);
        for (unsigned int j = 0; j < 3; ++j) {
            rVectorPreviousPosition(index * 3 + j, 0) = r_previous_position[j];
        }
    }

    /* Neighbour nodes */
    const SizeType number_of_neighbours = NumberOfActiveNeighbours(r_neighbour_nodes);

    if (number_of_neighbours == 6) { // All the possible neighbours
        for (unsigned int index = 0; index < number_of_nodes; ++index) {
            const array_1d<double, 3>& r_previous_position = r_neighbour_nodes[index].GetInitialPosition().Coordinates()
                                                           + r_neighbour_nodes[index].FastGetSolutionStepValue(DISPLACEMENT, 1);

            for (unsigned int j = 0; j < 3; ++j) {
                rVectorPreviousPosition(18 + index * 3 + j, 0) = r_previous_position[j];
            }
        }
    } else {
        for (unsigned int index = 0; index < number_of_nodes; ++index) {
            if (HasNeighbour(index, r_neighbour_nodes[index])) {
                const array_1d<double, 3>& r_previous_position = r_neighbour_nodes[index].GetInitialPosition().Coordinates()
                                                              + r_neighbour_nodes[index].FastGetSolutionStepValue(DISPLACEMENT, 1);
                for (unsigned int j = 0; j < 3; ++j) {
                    rVectorPreviousPosition(18 + index * 3 + j, 0) = r_previous_position[j];
                }
            } else {
                for (unsigned int j = 0; j < 3; ++j) {
                    rVectorPreviousPosition(18 + index * 3 + j, 0) = 0.0;
                }
            }
        }
    }

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

    // Auxiliary values
    BoundedMatrix<double, 1, 36 > B3;
    BoundedMatrix<double, 1,  6 > D3;

    if (mELementalFlags.Is(SolidShellElementSprism3D6N::EXPLICIT_RHS_COMPUTATION) ) { // Explicit calculation of the RHS
        // Kirchoff PK2 tangent tensor
        const Properties& r_properties = GetProperties();
        const double young_modulus = r_properties[YOUNG_MODULUS];
        const double poisson_ratio = r_properties[POISSON_RATIO];
        const double c1 = young_modulus / (( 1.0 + poisson_ratio ) * ( 1.0 - 2.0 * poisson_ratio ) );
        const double c2 = c1 * ( 1 - poisson_ratio );
        const double c3 = c1 * poisson_ratio;

        // Calculate EAS stiffness
        rEAS.mStiffAlpha += IntegrationWeight * ZetaGauss * ZetaGauss * rVariables.C[2] * (c2 * rVariables.C[2] + 2.0 * rVariables.StressVector[2]);

        // Assign D3
        D3(0, 0) = c3;
        D3(0, 1) = c3;
        D3(0, 2) = c2;
        D3(0, 3) = 0.0;
        D3(0, 4) = 0.0;
        D3(0, 5) = 0.0;

    } else { // Implicit solution uses the tangent tensor
        // Calculate EAS stiffness
        rEAS.mStiffAlpha += IntegrationWeight * ZetaGauss * ZetaGauss * rVariables.C[2] * (rVariables.ConstitutiveMatrix(2, 2) * rVariables.C[2] + 2.0 * rVariables.StressVector[2]);

        // Assign D3
        for (unsigned int i = 0; i < 6; ++i) {
            D3(0, i) = rVariables.ConstitutiveMatrix(2, i);
        }
    }
    for (unsigned int i = 0; i < 36; ++i) {
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
    ConstitutiveLaw::Parameters& rConstitutiveValues,
    const StressIntegratedComponents& rIntegratedStress,
    const CommonComponents& rCommonComponents,
    const CartesianDerivatives& rCartesianDerivatives,
    const EASComponents& rEAS,
    double& rAlphaEAS
    )
{
    /* Contributions of the stiffness matrix calculated on the reference configuration */
    if (rLocalSystem.CalculationFlags.Is( SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX_WITH_COMPONENTS )) {
        std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
        const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

        for (unsigned int i = 0; i < rLeftHandSideVariables.size(); ++i) {
            bool calculated = false;
            /* Calculate the Material Stiffness Matrix */
            if( rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX ) {
                /* Reading integration points */
                const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints(this->GetIntegrationMethod());

                for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                    const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

                    /* Assemble B */
                    this->CalculateDeformationMatrix(rVariables.B, rCommonComponents, zeta_gauss, rAlphaEAS);

                    // Compute element kinematics C, F ...
                    this->CalculateKinematics(rVariables, rCommonComponents, r_integration_points, point_number, rAlphaEAS, zeta_gauss);

                    // Set general variables to constitutive law parameters
                    this->SetGeneralVariables(rVariables, rConstitutiveValues, point_number);

                    // Compute stresses and constitutive parameters
                    mConstitutiveLawVector[point_number]->CalculateMaterialResponse(rConstitutiveValues, rVariables.StressMeasure);

                    // Calculating weights for integration on the "reference configuration"
                    const double integration_weight = r_integration_points[point_number].Weight() * rVariables.detJ;

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
            if (mELementalFlags.Is(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT)) {
                /* Apply EAS stabilization */
                ApplyEASLHS(rLeftHandSideMatrices[i], rEAS);
            }

            KRATOS_ERROR_IF_NOT(calculated) << " ELEMENT can not supply the required local system variable: " << rLeftHandSideVariables[i] << std::endl;
        }
    } else {
        MatrixType& LeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        /* Reading integration points */
        const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints(this->GetIntegrationMethod());

        /* Calculate the Material Stiffness Matrix */
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            const double zeta_gauss = 2.0 * r_integration_points[point_number].Z() - 1.0;

            /* Assemble B */
            this->CalculateDeformationMatrix(rVariables.B, rCommonComponents, zeta_gauss, rAlphaEAS);

            // Compute element kinematics C, F ...
            this->CalculateKinematics(rVariables, rCommonComponents, r_integration_points, point_number, rAlphaEAS, zeta_gauss);

            // Set general variables to constitutive law parameters
            this->SetGeneralVariables(rVariables, rConstitutiveValues, point_number);

            // Compute stresses and constitutive parameters
            mConstitutiveLawVector[point_number]->CalculateMaterialResponse(rConstitutiveValues, rVariables.StressMeasure);

            // Calculating weights for integration on the "reference configuration"
            const double integration_weight = r_integration_points[point_number].Weight() * rVariables.detJ;

            /* Operation performed: add Km to the LefsHandSideMatrix */
            this->CalculateAndAddKuum( LeftHandSideMatrix, rVariables, integration_weight);
        }

        /* Calculate the Geometric Stiffness Matrix */
        /* Operation performed: add Kg to the LefsHandSideMatrix */
        this->CalculateAndAddKuug( LeftHandSideMatrix, rIntegratedStress, rCartesianDerivatives );

        /* Implicit or explicit EAS update*/
        if (mELementalFlags.Is(SolidShellElementSprism3D6N::EAS_IMPLICIT_EXPLICIT)) {
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
    double& rAlphaEAS
    )
{
    /* Contribution of the internal and external forces */
    if (rLocalSystem.CalculationFlags.Is( SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS )) {
        std::vector<VectorType>& RightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
        const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
        for (unsigned int i = 0; i < rRightHandSideVariables.size(); ++i) {
            bool calculated = false;
            if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR ) {
                /* Operation performed: RightHandSideVector += ExtForce */
                this->CalculateAndAddExternalForces( RightHandSideVectors[i], rVariables, rVolumeForce );
                calculated = true;
            }

            if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR ) {
                /* Operation performed: RightHandSideVector -= IntForce */
                this->CalculateAndAddInternalForces( RightHandSideVectors[i], rIntegratedStress, rCommonComponents, rEAS, rAlphaEAS );
                calculated = true;
            }

            KRATOS_ERROR_IF_NOT(calculated) << " ELEMENT can not supply the required local system variable: " << rRightHandSideVariables[i] << std::endl;
        }
    } else {
        VectorType& RightHandSideVector = rLocalSystem.GetRightHandSideVector();

        /* Operation performed: RightHandSideVector += ExtForce */
        this->CalculateAndAddExternalForces( RightHandSideVector, rVariables, rVolumeForce );

        /* Operation performed: RightHandSideVector -= IntForce */
        this->CalculateAndAddInternalForces( RightHandSideVector, rIntegratedStress, rCommonComponents, rEAS, rAlphaEAS );
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

    unsigned int index_i, index_j;
    for (unsigned int i = 0; i < 36; ++i) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (index_i < 36) {
            for (unsigned int j = 0; j < 36; ++j) {
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

    /* Auxiliary stiffness matrix */
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

    unsigned int index_i, index_j;
    for (unsigned int i = 0; i < 36; ++i) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (index_i < 36) {
            for (unsigned int j = 0; j < 36; ++j) {
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
    unsigned int index_i, index_j;
    for (unsigned int i = 0; i < 36; ++i) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (index_i < 36) {
            for (unsigned int j = 0; j < 36; ++j) {
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
    double& rAlphaEAS
    )
{
    KRATOS_TRY;

    /* Calculate the RHS */
    noalias(rRHSFull) -= trans(rEAS.mHEAS) * rEAS.mRHSAlpha / rEAS.mStiffAlpha;

    /* Update ALPHA_EAS */
    rAlphaEAS -= rEAS.mRHSAlpha / rEAS.mStiffAlpha;

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

    // Number of nodes
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int index;
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        index = 3 * i;
        for (unsigned int j = 0; j < 3; ++j) {
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
    double& rAlphaEAS
    )
{
    KRATOS_TRY;

    BoundedMatrix<double, 36, 1 > rhs_full = ZeroMatrix(36, 1);

    unsigned int aux_index = 0;
    for (unsigned int i = 0; i < 18; ++i) {
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
    ApplyEASRHS(rhs_full, rEAS, rAlphaEAS);

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    unsigned int index_i;
    for (unsigned int i = 0; i < 36; ++i) {
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
    ConstitutiveLaw::Parameters& rConstitutiveValues,
    const unsigned int PointNumber
    )
{
    KRATOS_ERROR_IF(rVariables.detF < 0) << "SPRISM ELEMENT: " << this->Id() << "INVERTED: |F| < 0  detF = " << rVariables.detF << std::endl;

    // Compute total F: FT
    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );

    rConstitutiveValues.SetDeterminantF(rVariables.detFT);
    rConstitutiveValues.SetDeformationGradientF(rVariables.FT);
    rConstitutiveValues.SetStrainVector(rVariables.StrainVector);
    rConstitutiveValues.SetStressVector(rVariables.StressVector);
    rConstitutiveValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);

    // Adding the standard prism shape functions
    rConstitutiveValues.SetShapeFunctionsValues(rVariables.N);
    rConstitutiveValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
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
    const WeakPointerVectorNodesType& r_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(r_neighbour_nodes);
    const unsigned int mat_size = number_of_nodes * 3;

    if (rCalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_LHS_MATRIX)) {// Calculation of the matrix is required
        if (rLeftHandSideMatrix.size1() != mat_size) {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); // Resetting LHS
    }

    // Resizing as needed the RHS
    if (rCalculationFlags.Is(SolidShellElementSprism3D6N::COMPUTE_RHS_VECTOR)) { // Calculation of the matrix is required
        if (rRightHandSideVector.size() != mat_size) {
            rRightHandSideVector.resize(mat_size, false);
        }

        rRightHandSideVector = ZeroVector(mat_size); // Resetting RHS
    }
}

/******************************* COMPUTE KINEMATICS ********************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateKinematics(
    GeneralVariables& rVariables,
    const CommonComponents& rCommonComponents,
    const IntegrationPointsArrayType& rIntegrationPoints,
    const unsigned int PointNumber,
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

    if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        // PK2 stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

        // Jacobian Determinant for the isoparametric and numerical integration
        Matrix J0;
        auto& r_geometry = GetGeometry();
        GeometryUtils::JacobianOnInitialConfiguration(r_geometry, rIntegrationPoints[PointNumber], J0);
        rVariables.detJ = MathUtils<double>::Det(J0);
    } else {
        // Cauchy stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

        //Determinant of the Deformation Gradient F0
        rVariables.detF0 = MathUtils<double>::Det3(mAuxContainer[PointNumber]);
        noalias(rVariables.F0) = mAuxContainer[PointNumber];
    }

    this->CbartoFbar(rVariables, PointNumber);

    // Get the shape functions for the order of the integration method [N]
    const Matrix& r_N_container = rVariables.GetShapeFunctions();

    // Set Shape Functions Values for this integration point
    noalias(rVariables.N) = row( r_N_container, PointNumber);

    KRATOS_CATCH( "" );
}

/***************************** COMPUTE DELTA POSITION ******************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CalculateDeltaPosition(Matrix& rDeltaPosition)
{
    KRATOS_TRY;

    // Get geometry
    const auto& r_geometry = GetGeometry();

    // Iterate over the nodes to calculate the delta position
    const unsigned int number_of_nodes = r_geometry.size();
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3>& r_current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3>& r_previous_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, 1);

        for (unsigned int j = 0; j < 3; ++j) {
            rDeltaPosition(i,j) = r_current_displacement[j] - r_previous_displacement[j];
        }
    }

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::CbartoFbar(
    GeneralVariables& rVariables,
    const unsigned int PointNumber
    )
{
    KRATOS_TRY;

    /* We perform a polar decomposition of the CBar and F(regular) to obtain F_bar */

    // Assemble matrix C_bar
    const Matrix C_bar = MathUtils<double>::VectorToSymmetricTensor(rVariables.C);

    // Decompose matrix C_bar, get U_bar
    Matrix U_bar;
    MathUtils<double>::MatrixSquareRoot(C_bar, U_bar, 1e-24, 100);

    /* Decompose F */
    Matrix F = ZeroMatrix(3, 3);
    if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[PointNumber], mAuxContainer[PointNumber] );
    } else {
        // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
        Matrix InvJ(3, 3);
        MathUtils<double>::InvertMatrix( rVariables.J[PointNumber], InvJ, rVariables.detJ);

        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[PointNumber], InvJ );
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

    for (unsigned int index = 0; index < 9; ++index) {
        /* Element nodes */ // Note: It's important to consider the Voigt notation order considered in Kratos
        // Lower face
        rB(0, index)      = L1 * rCommonComponents.BMembraneLower(0, index);  // xx
        rB(1, index)      = L1 * rCommonComponents.BMembraneLower(1, index);  // yy
        rB(2, index)      = factor_eas * rCommonComponents.BNormal(0, index); // zz
        rB(3, index)      = L1 * rCommonComponents.BMembraneLower(2, index);  // xy
        rB(4, index)      = L1 * rCommonComponents.BShearLower(1, index) + L2 * rCommonComponents.BShearUpper(1, index); // yz
        rB(5, index)      = L1 * rCommonComponents.BShearLower(0, index) + L2 * rCommonComponents.BShearUpper(0, index); // xz
        // Upper face
        rB(0, index + 9)  = L2 * rCommonComponents.BMembraneUpper(0, index);      // xx
        rB(1, index + 9)  = L2 * rCommonComponents.BMembraneUpper(1, index);      // yy
        rB(2, index + 9)  = factor_eas * rCommonComponents.BNormal(0, index + 9); // zz
        rB(3, index + 9)  = L2 * rCommonComponents.BMembraneUpper(2, index);      // xy
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
    if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;
    } else {
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
    }

    // Get geometry
    const GeometryType& r_geometry = GetGeometry();

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
    rVariables.SetShapeFunctions(msGeometryData.ShapeFunctionsValues( this->GetIntegrationMethod() ));

    // Jacobians
    const IntegrationPointsArrayType& r_integration_points = msGeometryData.IntegrationPoints(this->GetIntegrationMethod());
    const IndexType integration_point_number = r_integration_points.size();
    rVariables.j.resize(integration_point_number, false);

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        r_geometry.Jacobian(rVariables.j[point_number], r_integration_points[point_number]);
    }

    if (mELementalFlags.IsNot(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        //Calculate Delta Position
        rVariables.J.resize(integration_point_number, false);
        Matrix delta_position(6, 3);
        this->CalculateDeltaPosition(delta_position);
        for (unsigned int point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            r_geometry.Jacobian(rVariables.J[point_number], r_integration_points[point_number], delta_position);
        }
    } else {
        rVariables.J.resize(1, false);
    }

    // Computing gradient
    GeometryType::ShapeFunctionsGradientsType DN_DX(integration_point_number, ZeroMatrix(6, 3));
    const GeometryType::ShapeFunctionsGradientsType& DN_De = msGeometryData.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    double detJ;
    Matrix inv_j;
    for (unsigned int i_point = 0; i_point < integration_point_number; ++i_point) {
        MathUtils<double>::InvertMatrix( rVariables.j[i_point], inv_j, detJ );
        noalias(DN_DX[i_point]) = prod(DN_De[i_point], inv_j);
    }

    // Setting shape functions gradients
    rVariables.SetShapeFunctionsGradients(DN_DX);
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::FinalizeStepVariables(
    GeneralVariables & rVariables,
    const unsigned int PointNumber
    )
{
    if (mELementalFlags.IsNot(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        // Update internal (historical) variables
        mAuxContainer[PointNumber] = prod(rVariables.F, rVariables.F0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SolidShellElementSprism3D6N::GetHistoricalVariables(
    GeneralVariables& rVariables,
    const unsigned int PointNumber
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

    if (mELementalFlags.Is(SolidShellElementSprism3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        rVolumeChange = 1.0;
    } else {
        rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);
    }

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

    // Ge the properties and geometry
    const auto& r_properties = GetProperties();
    const auto& r_geometry = GetGeometry();

    array_1d<double,3> volume_acceleration = ZeroVector(3);
    if (r_properties.Has( VOLUME_ACCELERATION )) {
        volume_acceleration = r_properties[VOLUME_ACCELERATION];
    } else if( r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        for (unsigned int i_node = 0; i_node < r_geometry.size(); ++i_node) {
            volume_acceleration += rVariables.N[i_node] * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
        }
    }

    // Compute volume change
    double volume_change;
    this->CalculateVolumeChange( volume_change, rVariables );

    rVolumeForce += volume_acceleration * IntegrationWeight * volume_change * r_properties[DENSITY];

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

const SolidShellElementSprism3D6N::IntegrationPointsContainerType SolidShellElementSprism3D6N::AllIntegrationPoints()
{
    IntegrationPointsContainerType r_integration_points =
    {
        {
            Quadrature<PrismGaussLegendreIntegrationPointsInAxis1, 3, IntegrationPoint<3>>::GenerateIntegrationPoints(),
            Quadrature<PrismGaussLegendreIntegrationPointsInAxis2, 3, IntegrationPoint<3>>::GenerateIntegrationPoints(),
            Quadrature<PrismGaussLegendreIntegrationPointsInAxis3, 3, IntegrationPoint<3>>::GenerateIntegrationPoints(),
            Quadrature<PrismGaussLegendreIntegrationPointsInAxis4, 3, IntegrationPoint<3>>::GenerateIntegrationPoints(),
            Quadrature<PrismGaussLegendreIntegrationPointsInAxis5, 3, IntegrationPoint<3>>::GenerateIntegrationPoints()
        }
    };
    return r_integration_points;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix SolidShellElementSprism3D6N::CalculateShapeFunctionsIntegrationPointsValues(const int ThisMethod)
{
    IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
    IntegrationPointsArrayType r_integration_points = all_integration_points[ThisMethod];

    // Number of integration points
    const unsigned int integration_points_number = r_integration_points.size();

    // Number of nodes in current geometry
    const unsigned int points_number = 6;

    // Setting up return matrix
    Matrix shape_function_values(integration_points_number, points_number);

    // Loop over all integration points
    for (unsigned int pnt = 0; pnt < integration_points_number; pnt++) {
        shape_function_values(pnt, 0) = (1.0
                                         - r_integration_points[pnt].X()
                                         - r_integration_points[pnt].Y()
                                         - r_integration_points[pnt].Z()
                                         + ( r_integration_points[pnt].X() * r_integration_points[pnt].Z())
                                         + ( r_integration_points[pnt].Y() * r_integration_points[pnt].Z()));
        shape_function_values(pnt, 1) = r_integration_points[pnt].X()
                                         - ( r_integration_points[pnt].X() * r_integration_points[pnt].Z());
        shape_function_values(pnt, 2) = r_integration_points[pnt].Y()
                                         - ( r_integration_points[pnt].Y() * r_integration_points[pnt].Z());
        shape_function_values(pnt, 3) = r_integration_points[pnt].Z()
                                         - ( r_integration_points[pnt].X() * r_integration_points[pnt].Z())
                                         - ( r_integration_points[pnt].Y() * r_integration_points[pnt].Z());
        shape_function_values(pnt, 4) = ( r_integration_points[pnt].X() * r_integration_points[pnt].Z());
        shape_function_values(pnt, 5) = ( r_integration_points[pnt].Y() * r_integration_points[pnt].Z());
    }

    return shape_function_values;
}

/***********************************************************************************/
/***********************************************************************************/

const SolidShellElementSprism3D6N::ShapeFunctionsValuesContainerType SolidShellElementSprism3D6N::AllShapeFunctionsValues()
{
    ShapeFunctionsValuesContainerType shape_functions_values =
    {
        {
            CalculateShapeFunctionsIntegrationPointsValues(0),
            CalculateShapeFunctionsIntegrationPointsValues(1),
            CalculateShapeFunctionsIntegrationPointsValues(2),
            CalculateShapeFunctionsIntegrationPointsValues(3),
            CalculateShapeFunctionsIntegrationPointsValues(4)
        }
    };
    return shape_functions_values;
}

/***********************************************************************************/
/***********************************************************************************/

SolidShellElementSprism3D6N::ShapeFunctionsGradientsType SolidShellElementSprism3D6N::CalculateShapeFunctionsIntegrationPointsLocalGradients(const int ThisMethod )
{
    IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
    IntegrationPointsArrayType r_integration_points = all_integration_points[ThisMethod];

    // Number of integration points
    const unsigned int integration_points_number = r_integration_points.size();
    ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );

    // Initialising container
    // Loop over all integration points
    Matrix result(6, 3);
    for (unsigned int pnt = 0; pnt < integration_points_number; pnt++ ) {
        result(0, 0) = -1.0 + r_integration_points[pnt].Z();
        result(0, 1) = -1.0 + r_integration_points[pnt].Z();
        result(0, 2) = -1.0 + r_integration_points[pnt].X() + r_integration_points[pnt].Y();
        result(1, 0) =  1.0 - r_integration_points[pnt].Z();
        result(1, 1) =  0.0;
        result(1, 2) =  -r_integration_points[pnt].X();
        result(2, 0) =  0.0;
        result(2, 1) =  1.0 - r_integration_points[pnt].Z();
        result(2, 2) =  -r_integration_points[pnt].Y();
        result(3, 0) =  -r_integration_points[pnt].Z();
        result(3, 1) =  -r_integration_points[pnt].Z();
        result(3, 2) =  1.0 - r_integration_points[pnt].X() - r_integration_points[pnt].Y();
        result(4, 0) =  r_integration_points[pnt].Z();
        result(4, 1) =  0.0;
        result(4, 2) =  r_integration_points[pnt].X();
        result(5, 0) =  0.0;
        result(5, 1) =  r_integration_points[pnt].Z();
        result(5, 2) =  r_integration_points[pnt].Y();
        d_shape_f_values[pnt] = result;
    }

    return d_shape_f_values;
}

/***********************************************************************************/
/***********************************************************************************/

const SolidShellElementSprism3D6N::ShapeFunctionsLocalGradientsContainerType SolidShellElementSprism3D6N::AllShapeFunctionsLocalGradients()
{
    ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
    {
        {
            CalculateShapeFunctionsIntegrationPointsLocalGradients(0),
            CalculateShapeFunctionsIntegrationPointsLocalGradients(1),
            CalculateShapeFunctionsIntegrationPointsLocalGradients(2),
            CalculateShapeFunctionsIntegrationPointsLocalGradients(3),
            CalculateShapeFunctionsIntegrationPointsLocalGradients(4)
        }
    };
    return shape_functions_local_gradients;
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

/***********************************************************************************/
/***********************************************************************************/

const PrismGaussLegendreIntegrationPointsInAxis1::IntegrationPointsArrayType& PrismGaussLegendreIntegrationPointsInAxis1::IntegrationPoints()
{
    static constexpr double one_third = 1.0 / 3.0;
    static const IntegrationPointsArrayType s_integration_points{{
        IntegrationPointType( one_third , one_third, 0.2113248654051871177454 , 0.25 ),
        IntegrationPointType( one_third , one_third, 0.7886751345948128822546 , 0.25 )
    }};
    return s_integration_points;
}

/***********************************************************************************/
/***********************************************************************************/

const PrismGaussLegendreIntegrationPointsInAxis2::IntegrationPointsArrayType& PrismGaussLegendreIntegrationPointsInAxis2::IntegrationPoints()
{
    static constexpr double one_third = 1.0 / 3.0;
    static const IntegrationPointsArrayType s_integration_points{{
        IntegrationPointType( one_third , one_third, 0.1127016653792583114821 , 0.1388888888888888889 ),
        IntegrationPointType( one_third , one_third, 0.5 , 0.2222222222222222222 ),
        IntegrationPointType( one_third , one_third, 0.8872983346207416885180 , 0.1388888888888888889 )
    }};
    return s_integration_points;
}

/***********************************************************************************/
/***********************************************************************************/

const PrismGaussLegendreIntegrationPointsInAxis3::IntegrationPointsArrayType& PrismGaussLegendreIntegrationPointsInAxis3::IntegrationPoints()
{
    static constexpr double one_third = 1.0 / 3.0;
    static const IntegrationPointsArrayType s_integration_points{{
        IntegrationPointType( one_third , one_third, 0.0469100770306680036012 , 0.0592317212640472718 ),
        IntegrationPointType( one_third , one_third, 0.2307653449471584544819 , 0.1196571676248416170 ),
        IntegrationPointType( one_third , one_third, 0.5 , 0.1422222222222222222 ),
        IntegrationPointType( one_third , one_third, 0.7692346550528415455182 , 0.1196571676248416170 ),
        IntegrationPointType( one_third , one_third, 0.9530899229693319963988 , 0.0592317212640472718 )
    }};
    return s_integration_points;
}

/***********************************************************************************/
/***********************************************************************************/

const PrismGaussLegendreIntegrationPointsInAxis4::IntegrationPointsArrayType& PrismGaussLegendreIntegrationPointsInAxis4::IntegrationPoints()
{
    static constexpr double one_third = 1.0 / 3.0;
    static const IntegrationPointsArrayType s_integration_points{{
        IntegrationPointType( one_third , one_third, 0.0254460438286207377369 , 0.0261224489795918367347 ),
        IntegrationPointType( one_third , one_third, 0.1292344072003027800681 , 0.069926347872319166975  ),
        IntegrationPointType( one_third , one_third, 0.2970774243113014165467 , 0.09545751262627973624   ),
        IntegrationPointType( one_third , one_third, 0.5 , 0.1044897959183673469388 ),
        IntegrationPointType( one_third , one_third, 0.7029225756886985834533 , 0.09545751262627973624   ),
        IntegrationPointType( one_third , one_third, 0.8707655927996972199320 , 0.069926347872319166975  ),
        IntegrationPointType( one_third , one_third, 0.9745539561713792622631 , 0.0261224489795918367347 )
    }};
    return s_integration_points;
}

/***********************************************************************************/
/***********************************************************************************/

const PrismGaussLegendreIntegrationPointsInAxis5::IntegrationPointsArrayType& PrismGaussLegendreIntegrationPointsInAxis5::IntegrationPoints()
{
    static constexpr double one_third = 1.0 / 3.0;
    static const IntegrationPointsArrayType s_integration_points{{
        IntegrationPointType( one_third , one_third, 0.0108856709269715035981 , 0.0139171417790434166207 ),
        IntegrationPointType( one_third , one_third, 0.0564687001159523504624 , 0.031395092366226156159  ),
        IntegrationPointType( one_third , one_third, 0.1349239972129753379533 , 0.0465725527319335628565 ),
        IntegrationPointType( one_third , one_third, 0.2404519353965940920372 , 0.058298441147997619980  ),
        IntegrationPointType( one_third , one_third, 0.3652284220238275138343 , 0.065701136127561665545  ),
        IntegrationPointType( one_third , one_third, 0.5 , 0.0682312716944751576786 ),
        IntegrationPointType( one_third , one_third, 0.6347715779761724861657 , 0.065701136127561665545  ),
        IntegrationPointType( one_third , one_third, 0.7595480646034059079628 , 0.058298441147997619980  ),
        IntegrationPointType( one_third , one_third, 0.8650760027870246620467 , 0.0465725527319335628565 ),
        IntegrationPointType( one_third , one_third, 0.9435312998840476495376 , 0.031395092366226156159  ),
        IntegrationPointType( one_third , one_third, 0.9891143290730284964020 , 0.0139171417790434166207 )
    }};
    return s_integration_points;
}

} // Namespace Kratos.
