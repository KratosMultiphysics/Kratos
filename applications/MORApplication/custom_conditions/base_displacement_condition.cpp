// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: MOR_application/license.txt
//
//  Main authors:    Ramsubramanian
//                   Ricky Aristio
//

// System includes

// External includes

// Project includes
#include "custom_conditions/base_displacement_condition.h"
#include "includes/checks.h"

namespace Kratos
{
BaseDisplacementCondition::BaseDisplacementCondition(BaseDisplacementCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

BaseDisplacementCondition& BaseDisplacementCondition::operator=(BaseDisplacementCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer BaseDisplacementCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<BaseDisplacementCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer BaseDisplacementCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<BaseDisplacementCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer BaseDisplacementCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<BaseDisplacementCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseDisplacementCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    
    if (rResult.size() != number_of_nodes) {
        rResult.resize(number_of_nodes, false);
    }

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(ACOUSTIC_PRESSURE);

    
        for (SizeType i = 0; i < number_of_nodes; ++i) {
           
            rResult[i] = GetGeometry()[i].GetDof(ACOUSTIC_PRESSURE,pos).EquationId();
           
        }
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
void BaseDisplacementCondition::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    ElementalDofList.resize(0);
    ElementalDofList.reserve(number_of_nodes);
    for (SizeType i = 0; i < number_of_nodes; ++i) {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof(ACOUSTIC_PRESSURE));  
    }
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseDisplacementCondition::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType mat_size = number_of_nodes;

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const double r_acousticpressure = GetGeometry()[i].FastGetSolutionStepValue(ACOUSTIC_PRESSURE, Step);
        rValues[i] = r_acousticpressure;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseDisplacementCondition::GetFirstDerivativesVector(
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

void BaseDisplacementCondition::GetSecondDerivativesVector(
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

void BaseDisplacementCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = false;
    const bool calculate_mass_matrix_flag = false;
    const bool calculate_residual_vector_flag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag, calculate_mass_matrix_flag );
}

/***********************************************************************************/
/***********************************************************************************/
void BaseDisplacementCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    //calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_mass_matrix_flag = true;
    const bool calculate_residual_vector_flag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag, calculate_mass_matrix_flag );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseDisplacementCondition::CalculateMassMatrix(
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

void BaseDisplacementCondition::CalculateDampingMatrix(
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

void BaseDisplacementCondition :: CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag,
        const bool CalculateMassMatrixFlag
        )
{
    KRATOS_ERROR << "You are calling the CalculateAll from the base class for loads" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

int BaseDisplacementCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(ACOUSTIC_PRESSURE)
    KRATOS_CHECK_VARIABLE_KEY(ACOUSTIC_DISPLACEMENT)
    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACOUSTIC_PRESSURE,r_node)
        KRATOS_CHECK_DOF_IN_NODE(ACOUSTIC_PRESSURE, r_node)
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

double BaseDisplacementCondition::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseDisplacementCondition::AddExplicitContribution(
    const VectorType& rRHS,
    const Variable<VectorType>& rRHSVariable,
    Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == ACOUSTIC_PRESSURE_RESIDUAL ) {
        for(SizeType i=0; i< number_of_nodes; ++i) {  

            double& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(ACOUSTIC_PRESSURE_RESIDUAL);
         
                r_force_residual += rRHS[i];          
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseDisplacementCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseDisplacementCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos
