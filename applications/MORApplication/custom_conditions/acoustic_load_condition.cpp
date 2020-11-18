//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt

// System includes


// External includes


// Project includes
#include "custom_conditions/acoustic_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"
#include "includes/checks.h"

#include "mor_application_variables.h"


namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

AcousticLoadCondition::AcousticLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

AcousticLoadCondition::AcousticLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer AcousticLoadCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticLoadCondition>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer AcousticLoadCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer AcousticLoadCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<AcousticLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

AcousticLoadCondition::~AcousticLoadCondition()
{
}

/***********************************************************************************/
/***********************************************************************************/
/**
 * @brief Sets on rResult the ID's of the element degrees of freedom
 * The dofs are ordered for each node as displacement - pressure - (rotation)
 * @param rResult The vector containing the equation id
 * @param rCurrentProcessInfo The current process info instance
 */
void AcousticLoadCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != number_of_nodes) {
        rResult.resize(number_of_nodes, false);
    }

    for( SizeType i=0; i<number_of_nodes; ++i ) {
        rResult[i] = GetGeometry()[i].GetDof(PRESSURE, pos + i).EquationId();
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
void AcousticLoadCondition::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    ElementalDofList.resize(0);
    ElementalDofList.reserve(number_of_nodes);

    for( SizeType i=0; i<number_of_nodes; ++i ) {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof(PRESSURE) );
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticLoadCondition::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticLoadCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticLoadCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

int AcousticLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_node)

        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node)
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

double AcousticLoadCondition::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticLoadCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    if ( rRightHandSideVector.size( ) != number_of_nodes ) {
        rRightHandSideVector.resize( number_of_nodes, false );
    }
    noalias( rRightHandSideVector ) = ZeroVector( number_of_nodes );

     // Reading integration points and local gradients
    const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
    const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

    double load_on_condition = 0.0;
    if( this->Has( ACOUSTIC_LOAD ) ) {
        load_on_condition += this->GetValue( ACOUSTIC_LOAD );
    }
    if( this->GetProperties().Has(ACOUSTIC_LOAD) ) {
        load_on_condition += this->GetProperties()[ACOUSTIC_LOAD];
    }
    const double frequency2 = std::pow( rCurrentProcessInfo[FREQUENCY], 2 );

    Vector load(number_of_nodes, load_on_condition);
    for( IndexType i=0; i<number_of_nodes; ++i ) {
        if( r_geometry[i].SolutionStepsDataHas( ACOUSTIC_LOAD ) ) {
            load[i] += r_geometry[i].FastGetSolutionStepValue( ACOUSTIC_LOAD );
        }
        load[i] *= frequency2;
    }

    for( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

        const double detJ = r_geometry.DeterminantOfJacobian( integration_points[point_number] );
        const double integration_weight =
            GetIntegrationWeight(integration_points, point_number, detJ);

        const Vector& rN = row(Ncontainer,point_number);

        for( IndexType i=0; i<number_of_nodes; ++i ) {
            rRightHandSideVector[i] += integration_weight * rN[i] * load[i];
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
void AcousticLoadCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = GetGeometry().size();

    if( rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes ) {
        rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);

    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticLoadCondition::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = GetGeometry().size();

    if( rMassMatrix.size1() != number_of_nodes || rMassMatrix.size2() != number_of_nodes ) {
        rMassMatrix.resize(number_of_nodes, number_of_nodes, false);
    }

    noalias(rMassMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticLoadCondition::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = GetGeometry().size();

    if( rDampingMatrix.size1() != number_of_nodes || rDampingMatrix.size2() != number_of_nodes ) {
        rDampingMatrix.resize(number_of_nodes, number_of_nodes, false);
    }

    noalias(rDampingMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);
}

} // Namespace Kratos


