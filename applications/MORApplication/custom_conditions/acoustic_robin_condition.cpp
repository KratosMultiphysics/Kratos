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
#include "custom_conditions/acoustic_robin_condition.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"
#include "includes/checks.h"

#include "mor_application_variables.h"


namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

AcousticRobinCondition::AcousticRobinCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

AcousticRobinCondition::AcousticRobinCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer AcousticRobinCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticRobinCondition>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer AcousticRobinCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticRobinCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer AcousticRobinCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<AcousticRobinCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

AcousticRobinCondition::~AcousticRobinCondition()
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
void AcousticRobinCondition::EquationIdVector(
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
void AcousticRobinCondition::GetDofList(
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

void AcousticRobinCondition::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticRobinCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticRobinCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

int AcousticRobinCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE)
    KRATOS_CHECK_VARIABLE_KEY(ADMITTANCE)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADMITTANCE,r_node)

        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node)
    }

    // Check if the required variables are available
    KRATOS_ERROR_IF( !this->GetProperties().Has(ADMITTANCE) ) << "The required variable ADMITTANCE is not available in condition "
        << this->Id() << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

double AcousticRobinCondition::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticRobinCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = GetGeometry().size();

    if ( rRightHandSideVector.size( ) != number_of_nodes ) {
        rRightHandSideVector.resize( number_of_nodes, false );
    }
    noalias( rRightHandSideVector ) = ZeroVector( number_of_nodes );
}

/***********************************************************************************/
/***********************************************************************************/
void AcousticRobinCondition::CalculateLocalSystem(
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

    if ( rRightHandSideVector.size( ) != number_of_nodes ) {
        rRightHandSideVector.resize( number_of_nodes, false );
    }
    noalias( rRightHandSideVector ) = ZeroVector( number_of_nodes );
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticRobinCondition::CalculateMassMatrix(
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

void AcousticRobinCondition::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    if( rDampingMatrix.size1() != number_of_nodes || rDampingMatrix.size2() != number_of_nodes ) {
        rDampingMatrix.resize(number_of_nodes, number_of_nodes, false);
    }
    noalias(rDampingMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Reading integration points and local gradients
    const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
    const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

    const double admittance = this->GetProperties()[ADMITTANCE];

    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

        const double detJ = r_geometry.DeterminantOfJacobian( integration_points[point_number] );

        const double integration_weight =
            GetIntegrationWeight(integration_points, point_number, detJ);

        const Vector& rN = row(Ncontainer,point_number);

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {

            for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                const double NiNj_weight = rN[i] * rN[j] * integration_weight;
                rDampingMatrix( i, j ) += NiNj_weight * admittance;
            }
        }

    }
}

} // Namespace Kratos


