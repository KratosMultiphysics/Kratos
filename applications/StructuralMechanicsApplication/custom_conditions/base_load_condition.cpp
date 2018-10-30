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
#include "includes/checks.h"
#include "custom_conditions/base_load_condition.h"

namespace Kratos
{
BaseLoadCondition::BaseLoadCondition(BaseLoadCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

BaseLoadCondition& BaseLoadCondition::operator=(BaseLoadCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer BaseLoadCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<BaseLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer BaseLoadCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<BaseLoadCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer BaseLoadCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_shared<BaseLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::Initialize()
{
    // TODO: Add somethig if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add somethig if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add somethig if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add somethig if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add somethig if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType block_size = this->GetBlockSize();
    if (rResult.size() != dim * number_of_nodes) {
        rResult.resize(number_of_nodes * block_size, false);
    }

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    if(dim == 2) {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * block_size;
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            if (this->HasRotDof())
                rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z,pos + 2).EquationId();
        }
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * block_size;
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
void BaseLoadCondition::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim =  GetGeometry().WorkingSpaceDimension();
    const SizeType block_size = this->GetBlockSize();
    ElementalDofList.resize(0);
    ElementalDofList.reserve(number_of_nodes * block_size);

    if(dim == 2) {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            if (this->HasRotDof())
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Z));
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

void BaseLoadCondition::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dim;

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 > & r_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        SizeType index = i * dim;
        for(SizeType k = 0; k < dim; ++k) {
            rValues[index + k] = r_displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dim;

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 > & Velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * dim;
        for(SizeType k = 0; k<dim; ++k) {
            rValues[index + k] = Velocity[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dim;

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& r_acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * dim;
        for(SizeType k = 0; k < dim; ++k) {
            rValues[index + k] = r_acceleration[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::CalculateRightHandSide(
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
void BaseLoadCondition::CalculateLocalSystem(
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

void BaseLoadCondition::CalculateMassMatrix(
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

void BaseLoadCondition::CalculateDampingMatrix(
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

void BaseLoadCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_ERROR << "You are calling the CalculateAll from the base class for loads" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

int BaseLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    const SizeType number_of_nodes = this->GetGeometry().size();
    for ( SizeType i = 0; i < number_of_nodes; ++i ) {
        NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

double BaseLoadCondition::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    )
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::AddExplicitContribution(
    const VectorType& rRHS,
    const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double,3> >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
        for(SizeType i=0; i< number_of_nodes; ++i) {
            SizeType index = dimension * i;

            array_1d<double, 3 >& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for(SizeType j=0; j<dimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHS[index + j];
            }
        }
    }

    KRATOS_CATCH( "" )
}
/***********************************************************************************/
/***********************************************************************************/

void BaseLoadCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

void BaseLoadCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos
