// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Mahmoud Zidan
//

// System includes

// External includes

// Project includes
#include "custom_conditions/displacement_control_condition.h"
#include "includes/checks.h"

namespace Kratos
{

// Constructor void
DisplacementControlCondition::DisplacementControlCondition(IndexType NewId)
    : Condition(NewId) {}

// Constructor using an array of nodes
DisplacementControlCondition::DisplacementControlCondition(IndexType NewId, const NodesArrayType& rThisNodes)
    : Condition(NewId, rThisNodes) {}

DisplacementControlCondition::DisplacementControlCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry) {}

// Constructor using an array of nodes with properties
DisplacementControlCondition::DisplacementControlCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties) {}

///Copy constructor
DisplacementControlCondition::DisplacementControlCondition(DisplacementControlCondition const& rOther)
    : Condition(rOther) {}

// Destructor
DisplacementControlCondition::~DisplacementControlCondition() {}

/***********************************************************************************/
/***********************************************************************************/

/// Assignment operator.
DisplacementControlCondition& DisplacementControlCondition::operator=(DisplacementControlCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);
    Flags::operator=(rOther);
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DisplacementControlCondition::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<DisplacementControlCondition>(NewId, GetGeometry().Create(rThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DisplacementControlCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<DisplacementControlCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DisplacementControlCondition::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<DisplacementControlCondition>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_cond->SetData(GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType block_size = GetBlockSize();
    if (rResult.size() != number_of_nodes * block_size) {
        rResult.resize(number_of_nodes * block_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        int index = i * 2;
        if (GetValue(POINT_LOAD_X) != 0.0) // prescribed displacement in x direction
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        else if (GetValue(POINT_LOAD_Y) != 0.0) // prescribed displacement in y direction
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        else if (GetValue(POINT_LOAD_Z) != 0.0) // prescribed displacement in z direction
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(LOAD_FACTOR).EquationId();
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::GetDofList(
    DofsVectorType& rConditionlDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType block_size = GetBlockSize();

    if (rConditionlDofList.size() != number_of_nodes * block_size) {
        rConditionlDofList.resize(number_of_nodes * block_size);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        SizeType index = i * block_size;
        if (GetValue(POINT_LOAD_X) != 0.0) // prescribed displacement in x direction
            rConditionlDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        else if (GetValue(POINT_LOAD_Y) != 0.0) // prescribed displacement in y direction
            rConditionlDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        else if (GetValue(POINT_LOAD_Z) != 0.0) // prescribed displacement in z direction
            rConditionlDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
        rConditionlDofList[index + 1] = GetGeometry()[i].pGetDof(LOAD_FACTOR);
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType mat_size = number_of_nodes * GetBlockSize();

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        SizeType index = i * GetBlockSize();
        if (GetValue(POINT_LOAD_X) != 0.0) // prescribed displacement in x direction
            rValues[index] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
        else if (GetValue(POINT_LOAD_Y) != 0.0) // prescribed displacement in y direction
            rValues[index] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
        else if (GetValue(POINT_LOAD_Z) != 0.0) // prescribed displacement in z direction
            rValues[index] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z, Step);
        rValues[index] = GetGeometry()[i].FastGetSolutionStepValue(LOAD_FACTOR, Step);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = false;
    const bool calculate_residual_vector_flag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag );
}

/***********************************************************************************/
/***********************************************************************************/
void DisplacementControlCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    //calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_residual_vector_flag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag );
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    if(rMassMatrix.size1() != 0) {
        rMassMatrix.resize(0, 0, false);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    if(rDampingMatrix.size1() != 0) {
        rDampingMatrix.resize(0, 0, false);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    // currently only works for one node
    SizeType number_of_nodes = 1;
    SizeType mat_size = number_of_nodes * GetBlockSize();

    if ( CalculateStiffnessMatrixFlag == true )
    {
        if ( rLeftHandSideMatrix.size1() != mat_size)
        {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
        if (GetValue(POINT_LOAD_X) != 0.0) // prescribed displacement in x direction
            rLeftHandSideMatrix(0, 1) -= GetValue(POINT_LOAD_X);
        else if (GetValue(POINT_LOAD_Y) != 0.0) // prescribed displacement in y direction
            rLeftHandSideMatrix(0, 1) -= GetValue(POINT_LOAD_Y);
        else if (GetValue(POINT_LOAD_Z) != 0.0) // prescribed displacement in z direction
            rLeftHandSideMatrix(0, 1) -= GetValue(POINT_LOAD_Z);
        rLeftHandSideMatrix(1, 0) += 1.0;
    }

    if ( CalculateResidualVectorFlag == true )
    {
        if ( rRightHandSideVector.size() != mat_size )
        {
            rRightHandSideVector.resize(mat_size, false);
        }
        noalias(rRightHandSideVector) = ZeroVector(mat_size);

        rRightHandSideVector(0) += GetGeometry()[0].FastGetSolutionStepValue(LOAD_FACTOR);
        rRightHandSideVector(1) += GetValue(PRESCRIBED_DISPLACEMENT);
        if (GetValue(POINT_LOAD_X) != 0.0) // prescribed displacement in x direction
        {
            rRightHandSideVector(0) *= GetValue(POINT_LOAD_X);
            rRightHandSideVector(1) -= GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
        }
        else if (GetValue(POINT_LOAD_Y) != 0.0) // prescribed displacement in y direction
        {
            rRightHandSideVector(0) *= GetValue(POINT_LOAD_Y);
            rRightHandSideVector(1) -= GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
        }
        else if (GetValue(POINT_LOAD_Z) != 0.0) // prescribed displacement in z direction
        {
            rRightHandSideVector(0) *= GetValue(POINT_LOAD_Z);
            rRightHandSideVector(1) -= GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

int DisplacementControlCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(LOAD_FACTOR)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : GetGeometry().Points()) {
        KRATOS_CHECK_DOF_IN_NODE(LOAD_FACTOR, r_node)
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementControlCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos
