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
#include "custom_conditions/acoustic_structure_coupling_condition.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"
#include "includes/checks.h"


namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

template<std::size_t TDim>
AcousticStructureCouplingCondition<TDim>::AcousticStructureCouplingCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
AcousticStructureCouplingCondition<TDim>::AcousticStructureCouplingCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer AcousticStructureCouplingCondition<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticStructureCouplingCondition<TDim>>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer AcousticStructureCouplingCondition<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticStructureCouplingCondition<TDim>>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer AcousticStructureCouplingCondition<TDim>::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<AcousticStructureCouplingCondition<TDim>>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

template<std::size_t TDim>
AcousticStructureCouplingCondition<TDim>::~AcousticStructureCouplingCondition()
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
template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType block_size = this->GetBlockSize();
    if (rResult.size() != block_size * number_of_nodes) {
        rResult.resize(number_of_nodes * block_size, false);
    }

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
    const SizeType pos_p = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if(dim == 2) {
        std::cout << "dim2!\n";
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = this->HasRotDof() ? i * (block_size - 2) : i * (block_size - 1);

            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[i + dim * number_of_nodes] = GetGeometry()[i].GetDof(PRESSURE,pos_p).EquationId();

            if (this->HasRotDof())
                rResult[i + dim * number_of_nodes + number_of_nodes] = GetGeometry()[i].GetDof(ROTATION_Z,pos + 2).EquationId();
        }
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * (block_size - 1);

            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
            rResult[i + dim * number_of_nodes] = GetGeometry()[i].GetDof(PRESSURE,pos_p).EquationId();
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::GetDofList(
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
        }
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(PRESSURE));
        }
        if (this->HasRotDof()){
            for (SizeType i = 0; i < number_of_nodes; ++i) {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Z));
            }
        }
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(PRESSURE));
        }
        if (this->HasRotDof()){
            for (SizeType i = 0; i < number_of_nodes; ++i) {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Y));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Z));
            }
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

template<std::size_t TDim>
int AcousticStructureCouplingCondition<TDim>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        if( TDim == 3)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node)
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
double AcousticStructureCouplingCondition<TDim>::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    MatrixType temp = Matrix();

    const bool calculate_stiffness_matrix_flag = false;
    const bool calculate_mass_matrix_flag = false;
    const bool calculate_vector_flag = true;
    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo,
        calculate_stiffness_matrix_flag, calculate_mass_matrix_flag, calculate_vector_flag);
}

/***********************************************************************************/
/***********************************************************************************/
template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_mass_matrix_flag = false;
    const bool calculate_vector_flag = true;
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
        calculate_stiffness_matrix_flag, calculate_mass_matrix_flag, calculate_vector_flag);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const bool calculate_stiffness_matrix_flag = false;
    const bool calculate_mass_matrix_flag = true;
    const bool calculate_vector_flag = false;
    VectorType temp = Vector();
    CalculateAll(rMassMatrix, temp, rCurrentProcessInfo,
        calculate_stiffness_matrix_flag, calculate_mass_matrix_flag, calculate_vector_flag);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    if(rDampingMatrix.size1() != 0) {
        rDampingMatrix.resize(0, 0, false);
    }
}


template<std::size_t TDim>
void AcousticStructureCouplingCondition<TDim>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateMassMatrixFlag,
    const bool CalculateVectorFlag
    )
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag || CalculateMassMatrixFlag ) {
        if ( rLeftHandSideMatrix.size1() != mat_size ) {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateVectorFlag ) {
        if ( rRightHandSideVector.size( ) != mat_size ) {
            rRightHandSideVector.resize( mat_size, false );
        }
        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }

    if ( CalculateStiffnessMatrixFlag || CalculateMassMatrixFlag ) {

        // Reading integration points and local gradients
        const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
        const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

        Matrix J0(dimension, dimension);
        Matrix J(dimension, 1);
        array_1d<double, 3> tangent_xi, tangent_eta;
        array_1d<double, 3> normal;

        int sign;
        CalculateMassMatrixFlag ? sign = -1 : sign = 1;

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            r_geometry.Jacobian(J, point_number, integration_method);
            GetLocalAxis1(tangent_xi, J);
            GetLocalAxis2(tangent_eta, J);
            MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);

            const double detJ = r_geometry.DeterminantOfJacobian( integration_points[point_number] );

            const double integration_weight =
                GetIntegrationWeight(integration_points, point_number, detJ);

            const Vector& rN = row(Ncontainer,point_number);

            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                const SizeType index_i = i * dimension;

                for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                    const SizeType index_j = j + (dimension * number_of_nodes);
                    const double NiNj_weight = rN[i] * rN[j] * integration_weight;
                    for ( IndexType k = 0; k < dimension; ++k ) {
                        rLeftHandSideMatrix( index_i + k, index_j ) += NiNj_weight * normal(k) * sign;
                    }
                }
            }

        }
        if( CalculateMassMatrixFlag ) {
            rLeftHandSideMatrix = trans(rLeftHandSideMatrix);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<2>::GetLocalAxis1(
    array_1d<double, 3>& rLocalAxis,
    const Matrix& rJacobian
    ) const
{
    rLocalAxis[0] = rJacobian(0, 0);
    rLocalAxis[1] = rJacobian(1, 0);
    rLocalAxis[2] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<3>::GetLocalAxis1(
    array_1d<double, 3>& rLocalAxis,
    const Matrix& rJacobian
    ) const
{
    rLocalAxis[0] = rJacobian(0, 0);
    rLocalAxis[1] = rJacobian(1, 0);
    rLocalAxis[2] = rJacobian(2, 0);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<2>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis, const Matrix& rJacobian) const
{
    rLocalAxis[0] = 0.0;
    rLocalAxis[1] = 0.0;
    rLocalAxis[2] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<3>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis, const Matrix& rJacobian) const
{
    rLocalAxis[0] = rJacobian(0, 1);
    rLocalAxis[1] = rJacobian(1, 1);
    rLocalAxis[2] = rJacobian(2, 1);
}

/***********************************************************************************/
/***********************************************************************************/

template class AcousticStructureCouplingCondition<2>;
template class AcousticStructureCouplingCondition<3>;

} // Namespace Kratos


