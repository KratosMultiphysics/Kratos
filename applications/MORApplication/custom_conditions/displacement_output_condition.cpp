// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/displacement_output_condition.h"
#include "includes/checks.h"

namespace Kratos
{
// template< class TVariableType, TVariableType& TVariable >
DisplacementOutputCondition::DisplacementOutputCondition(DisplacementOutputCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

// template<class TVariable>
DisplacementOutputCondition& DisplacementOutputCondition::operator=(DisplacementOutputCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DisplacementOutputCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<DisplacementOutputCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DisplacementOutputCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<DisplacementOutputCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DisplacementOutputCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<DisplacementOutputCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementOutputCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    // const SizeType block_size = this->GetBlockSize();
    if (rResult.size() != dim * number_of_nodes) {
        rResult.resize(number_of_nodes * dim, false);
    }

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    if(dim == 2) {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * dim;
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            // if (this->HasRotDof())
            //     rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z,pos + 2).EquationId();
        }
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * dim;
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
        }
    }
    KRATOS_CATCH("")
}

// /***********************************************************************************/
// /***********************************************************************************/
void DisplacementOutputCondition::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim =  GetGeometry().WorkingSpaceDimension();
    // const SizeType block_size = this->GetBlockSize();
    ElementalDofList.resize(0);
    ElementalDofList.reserve(number_of_nodes * dim);

    if(dim == 2) {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementOutputCondition::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dim;

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 > & r_output = GetGeometry()[i].FastGetSolutionStepValue(COMPONENT_OUTPUT, Step);
        SizeType index = i * dim;
        for(SizeType k = 0; k < dim; ++k) {
            rValues[index + k] = r_output[k];
        }
    }
}

// /***********************************************************************************/
// /***********************************************************************************/

// // void DisplacementOutputCondition::GetFirstDerivativesVector(
// //     Vector& rValues,
// //     int Step
// //     )
// // {
// //     const SizeType number_of_nodes = GetGeometry().size();
// //     const SizeType dim = GetGeometry().WorkingSpaceDimension();
// //     const SizeType mat_size = number_of_nodes * dim;

// //     if (rValues.size() != mat_size) {
// //         rValues.resize(mat_size, false);
// //     }

// //     for (SizeType i = 0; i < number_of_nodes; ++i) {
// //         const array_1d<double, 3 > & Velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
// //         const SizeType index = i * dim;
// //         for(SizeType k = 0; k<dim; ++k) {
// //             rValues[index + k] = Velocity[k];
// //         }
// //     }
// // }

// /***********************************************************************************/
// /***********************************************************************************/

// // void DisplacementOutputCondition::GetSecondDerivativesVector(
// //     Vector& rValues,
// //     int Step
// //     )
// // {
// //     const SizeType number_of_nodes = GetGeometry().size();
// //     const SizeType dim = GetGeometry().WorkingSpaceDimension();
// //     const SizeType mat_size = number_of_nodes * dim;

// //     if (rValues.size() != mat_size) {
// //         rValues.resize(mat_size, false);
// //     }

// //     for (SizeType i = 0; i < number_of_nodes; ++i) {
// //         const array_1d<double, 3 >& r_acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
// //         const SizeType index = i * dim;
// //         for(SizeType k = 0; k < dim; ++k) {
// //             rValues[index + k] = r_acceleration[k];
// //         }
// //     }
// // }

// /***********************************************************************************/
// /***********************************************************************************/

void DisplacementOutputCondition::CalculateRightHandSide(
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
void DisplacementOutputCondition::CalculateLocalSystem(
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

void DisplacementOutputCondition::CalculateMassMatrix(
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

void DisplacementOutputCondition::CalculateDampingMatrix(
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

void DisplacementOutputCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int mat_size = number_of_nodes * dim;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mat_size )
        {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the vector is required
    {
        if ( rRightHandSideVector.size( ) != mat_size )
        {
            rRightHandSideVector.resize( mat_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }

    if( rCurrentProcessInfo[BUILD_LEVEL] == 301 )
    {
        array_1d<double, 3 > output_factor = ZeroVector(3);
        if( this->Has( COMPONENT_OUTPUT ) )
        {
            noalias(output_factor) = this->GetValue( COMPONENT_OUTPUT );
        }

        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            const unsigned int base = ii*dim;

            // if( GetGeometry()[ii].SolutionStepsDataHas( POINT_LOAD ) )
            // {
            //     noalias(output_factor) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_LOAD );
            // }

            for(unsigned int k = 0; k < dim; ++k)
            {
                rRightHandSideVector[base + k] += output_factor[k];
            }
        }
    }

    KRATOS_CATCH( "" )
}

// /***********************************************************************************/
// /***********************************************************************************/

// template< class TVariableType, TVariableType& TVariable >
int DisplacementOutputCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(COMPONENT_OUTPUT)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(COMPONENT_OUTPUT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    }

    return 0;
}

// /***********************************************************************************/
// /***********************************************************************************/

// // double DisplacementOutputCondition::GetIntegrationWeight(
// //     const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
// //     const SizeType PointNumber,
// //     const double detJ
// //     ) const
// // {
// //     return IntegrationPoints[PointNumber].Weight() * detJ;
// // }

// /***********************************************************************************/
// /***********************************************************************************/

// // void DisplacementOutputCondition::AddExplicitContribution(
// //     const VectorType& rRHS,
// //     const Variable<VectorType>& rRHSVariable,
// //     Variable<array_1d<double,3> >& rDestinationVariable,
// //     const ProcessInfo& rCurrentProcessInfo
// //     )
// // {
// //     KRATOS_TRY;

// //     const SizeType number_of_nodes = GetGeometry().PointsNumber();
// //     const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

// //     if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
// //         for(SizeType i=0; i< number_of_nodes; ++i) {
// //             SizeType index = dimension * i;

// //             array_1d<double, 3 >& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
// //             for(SizeType j=0; j<dimension; ++j) {
// //                 #pragma omp atomic
// //                 r_force_residual[j] += rRHS[index + j];
// //             }
// //         }
// //     }

// //     KRATOS_CATCH( "" )
// // }

/***********************************************************************************/
/***********************************************************************************/

void DisplacementOutputCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void DisplacementOutputCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos
