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
#include "custom_conditions/pressure_output_condition.h"
#include "includes/checks.h"

namespace Kratos
{
// template< class TVariableType, TVariableType& TVariable >
PressureOutputCondition::PressureOutputCondition(PressureOutputCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

// template<class TVariable>
PressureOutputCondition& PressureOutputCondition::operator=(PressureOutputCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PressureOutputCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<PressureOutputCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PressureOutputCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<PressureOutputCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PressureOutputCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<PressureOutputCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void PressureOutputCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    if (rResult.size() != number_of_nodes) {
        rResult.resize(number_of_nodes, false);
    }

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        rResult[i] = GetGeometry()[i].GetDof(PRESSURE,pos).EquationId();
    }
    KRATOS_CATCH("")
}

// /***********************************************************************************/
// /***********************************************************************************/
void PressureOutputCondition::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();

    ElementalDofList.resize(0);
    ElementalDofList.reserve(number_of_nodes);

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof(PRESSURE));
    }
    KRATOS_CATCH("")
}

// /***********************************************************************************/
// /***********************************************************************************/

void PressureOutputCondition::CalculateRightHandSide(
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
void PressureOutputCondition::CalculateLocalSystem(
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

void PressureOutputCondition::CalculateMassMatrix(
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

void PressureOutputCondition::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    if(rDampingMatrix.size1() != 0) {
        rDampingMatrix.resize(0, 0, false);
    }
}

// /***********************************************************************************/
// /***********************************************************************************/

void PressureOutputCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();

    // Resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != number_of_nodes )
        {
            rLeftHandSideMatrix.resize( number_of_nodes, number_of_nodes, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( number_of_nodes, number_of_nodes ); //resetting LHS
    }

    //resize and compute the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the vector is required
    {
        if ( rRightHandSideVector.size( ) != number_of_nodes )
        {
            rRightHandSideVector.resize( number_of_nodes, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( number_of_nodes ); //resetting RHS

        if( rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 301 ) {
            for ( SizeType ii = 0; ii < number_of_nodes; ++ii) {
                double output_factor = 0.0;
                if( this->Has(SCALAR_OUTPUT) ) {
                    output_factor = this->GetValue( SCALAR_OUTPUT );
                } else if( GetGeometry()[ii].SolutionStepsDataHas( SCALAR_OUTPUT ) ) {
                    output_factor = GetGeometry()[ii].FastGetSolutionStepValue( SCALAR_OUTPUT );
                }

                rRightHandSideVector[ii] += output_factor;
            }
        }
    }


    KRATOS_CATCH( "" )
}

// /***********************************************************************************/
// /***********************************************************************************/

int PressureOutputCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(SCALAR_OUTPUT)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        // KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(SCALAR_OUTPUT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node)
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void PressureOutputCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void PressureOutputCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos
