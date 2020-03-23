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
#include "custom_conditions/acoustic_structure_line_condition.h"
#include "utilities/math_utils.h"
#include "utilities/beam_math_utilities.hpp"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"
#include "includes/checks.h"


namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

template<std::size_t TDim>
AcousticStructureLineCondition<TDim>::AcousticStructureLineCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
AcousticStructureLineCondition<TDim>::AcousticStructureLineCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer AcousticStructureLineCondition<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticStructureLineCondition<TDim>>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer AcousticStructureLineCondition<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticStructureLineCondition<TDim>>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer AcousticStructureLineCondition<TDim>::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<AcousticStructureLineCondition<TDim>>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

template<std::size_t TDim>
AcousticStructureLineCondition<TDim>::~AcousticStructureLineCondition()
{
}

/***********************************************************************************/
/***********************************************************************************/

// template<std::size_t TDim>
// void AcousticStructureLineCondition<TDim>::GetValueOnIntegrationPoints(
//     const Variable<array_1d<double, 3>>& rVariable,
//     std::vector< array_1d<double, 3>>& rOutput,
//     const ProcessInfo& rCurrentProcessInfo
//     )
// {
//     KRATOS_TRY;

//     this->CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo );

//     KRATOS_CATCH( "" );
// }

/***********************************************************************************/
/***********************************************************************************/

// template<std::size_t TDim>
// void AcousticStructureLineCondition<TDim>::CalculateOnIntegrationPoints(
//     const Variable<array_1d<double, 3>>& rVariable,
//     std::vector< array_1d<double, 3>>& rOutput,
//     const ProcessInfo& rCurrentProcessInfo
//     )
// {
//     KRATOS_TRY;

//     const auto& r_geometry = this->GetGeometry();
//     const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
//     const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);

//     if ( rOutput.size() != r_integration_points.size() )
//         rOutput.resize( r_integration_points.size() );

//     if (rVariable == NORMAL) {
//         // Declaring tangent and Jacobian
//         array_1d<double, 3> tangent_xi, tangent_eta;
//         Matrix J(TDim, 1);

//         // Getting LOCAL_AXIS_2
//         GetLocalAxis2(tangent_eta);

//         // Iterate over the Gauss points
//         for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
//             r_geometry.Jacobian(J, point_number, integration_method);

//             // Definition of the tangent
//             GetLocalAxis1(tangent_xi, J);

//             // Computing normal
//             MathUtils<double>::UnitCrossProduct(rOutput[point_number], tangent_xi, tangent_eta);
//         }
//     } else {
//         for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
//             rOutput[point_number] = ZeroVector(3);
//         }
//     }

//     KRATOS_CATCH( "" );
// }

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::EquationIdVector(
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
            // const SizeType index = i * (block_size - 1);
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[i + dim * number_of_nodes] = GetGeometry()[i].GetDof(PRESSURE,pos_p).EquationId();

            if (this->HasRotDof())
                rResult[i + dim * number_of_nodes + number_of_nodes] = GetGeometry()[i].GetDof(ROTATION_Z,pos + 2).EquationId();
        }
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * block_size;
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(PRESSURE,pos_p).EquationId();
        }
    }
    std::cout << "condition equation id vector\n";
    KRATOS_WATCH(rResult)
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    std::cout << "GetDofListGetDofListGetDofListGetDofListGetDofListGetDofList!\n";

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
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(PRESSURE));
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    std::cout << "GETVALUESVECTOR-----------------------\n";
    // const SizeType number_of_nodes = GetGeometry().size();
    // const SizeType dim = GetGeometry().WorkingSpaceDimension();
    // const SizeType mat_size = number_of_nodes * dim;

    // if (rValues.size() != mat_size) {
    //     rValues.resize(mat_size, false);
    // }

    // for (SizeType i = 0; i < number_of_nodes; ++i) {
    //     const array_1d<double, 3 > & r_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
    //     SizeType index = i * dim;
    //     for(SizeType k = 0; k < dim; ++k) {
    //         rValues[index + k] = r_displacement[k];
    //     }
    // }
}

/***********************************************************************************/
/***********************************************************************************/

// template<std::size_t TDim>
// void AcousticStructureLineCondition<TDim>::GetFirstDerivativesVector(
//     Vector& rValues,
//     int Step
//     )
// {
//     const SizeType number_of_nodes = GetGeometry().size();
//     const SizeType dim = GetGeometry().WorkingSpaceDimension();
//     const SizeType mat_size = number_of_nodes * dim;

//     if (rValues.size() != mat_size) {
//         rValues.resize(mat_size, false);
//     }

//     for (SizeType i = 0; i < number_of_nodes; ++i) {
//         const array_1d<double, 3 > & Velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
//         const SizeType index = i * dim;
//         for(SizeType k = 0; k<dim; ++k) {
//             rValues[index + k] = Velocity[k];
//         }
//     }
// }

/***********************************************************************************/
/***********************************************************************************/

// template<std::size_t TDim>
// void AcousticStructureLineCondition<TDim>::GetSecondDerivativesVector(
//     Vector& rValues,
//     int Step
//     )
// {
//     const SizeType number_of_nodes = GetGeometry().size();
//     const SizeType dim = GetGeometry().WorkingSpaceDimension();
//     const SizeType mat_size = number_of_nodes * dim;

//     if (rValues.size() != mat_size) {
//         rValues.resize(mat_size, false);
//     }

//     for (SizeType i = 0; i < number_of_nodes; ++i) {
//         const array_1d<double, 3 >& r_acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
//         const SizeType index = i * dim;
//         for(SizeType k = 0; k < dim; ++k) {
//             rValues[index + k] = r_acceleration[k];
//         }
//     }
// }

template<std::size_t TDim>
int AcousticStructureLineCondition<TDim>::Check( const ProcessInfo& rCurrentProcessInfo )
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
double AcousticStructureLineCondition<TDim>::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

// template<std::size_t TDim>
// void AcousticStructureLineCondition<TDim>::AddExplicitContribution(
//     const VectorType& rRHS,
//     const Variable<VectorType>& rRHSVariable,
//     Variable<array_1d<double,3> >& rDestinationVariable,
//     const ProcessInfo& rCurrentProcessInfo
//     )
// {
//     KRATOS_TRY;

//     const SizeType number_of_nodes = GetGeometry().PointsNumber();
//     const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

//     if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
//         for(SizeType i=0; i< number_of_nodes; ++i) {
//             SizeType index = dimension * i;

//             array_1d<double, 3 >& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
//             for(SizeType j=0; j<dimension; ++j) {
//                 #pragma omp atomic
//                 r_force_residual[j] += rRHS[index + j];
//             }
//         }
//     }

//     KRATOS_CATCH( "" )
// }

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // std::cout << "RESIZE AND SET TO ZERO!!\n";
    // Calculation flags
    // const bool calculate_stiffness_matrix_flag = false;
    // const bool calculate_residual_vector_flag = true;
    MatrixType temp = Matrix();

    // CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag );

    const bool calculate_stiffness_matrix_flag = false;
    const bool calculate_mass_matrix_flag = false;
    const bool calculate_vector_flag = true;
    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo,
        calculate_stiffness_matrix_flag, calculate_mass_matrix_flag, calculate_vector_flag);
}

/***********************************************************************************/
/***********************************************************************************/
template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::CalculateLocalSystem(
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
    //calculation flags
    // const bool calculate_stiffness_matrix_flag = true;
    // const bool calculate_residual_vector_flag = true;

    /*
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    // if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
    if ( rLeftHandSideMatrix.size1() != mat_size ) {
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );
    }
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    // }
    KRATOS_WATCH(mat_size)

    // // Resizing as needed the RHS
    // if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
    //     if ( rRightHandSideVector.size( ) != mat_size ) {
    //         rRightHandSideVector.resize( mat_size, false );
    //     }
    //     noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    // }


    // Reading integration points and local gradients
    const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
    const GeometryType::ShapeFunctionsGradientsType& rDN_De = r_geometry.ShapeFunctionsLocalGradients(integration_method);
    const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

    KRATOS_WATCH(integration_points)
    KRATOS_WATCH(rDN_De)
    KRATOS_WATCH(Ncontainer)

    Matrix J0(dimension, dimension);
    Matrix J(dimension, 1);
    array_1d<double, 3> tangent_xi, tangent_eta;
    array_1d<double, 3> normal;


    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        r_geometry.Jacobian(J, point_number, integration_method);
        GetLocalAxis1(tangent_xi, J);
        GetLocalAxis2(tangent_eta);
        MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);
        // KRATOS_WATCH(normal)
        // GeometryUtils::JacobianOnInitialConfiguration(
        //     r_geometry, integration_points[point_number], J0);
        // KRATOS_WATCH(J0)
        // KRATOS_WATCH(J1)
        const double detJ = r_geometry.DeterminantOfJacobian( integration_points[point_number] );
        // const double detJ = MathUtils<double>::DetMat(J);
        // KRATOS_WATCH(J)
        // KRATOS_WATCH(detJ)
        const double integration_weight =
            GetIntegrationWeight(integration_points, point_number, detJ);
        const Vector& rN = row(Ncontainer,point_number);

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const SizeType index_i = i * dimension;

            for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                const SizeType index_j = j + (dimension * number_of_nodes);
                const double NiNj_weight = rN[i] * rN[j] * integration_weight;

                for ( IndexType k = 0; k < dimension; ++k )
                    rLeftHandSideMatrix( index_i + k, index_j + k ) += NiNj_weight * normal(k);
            }
        }

        KRATOS_WATCH(rLeftHandSideMatrix)
    }
    */
    // // Sizing work matrices
    // Vector pressure_on_nodes = ZeroVector( number_of_nodes );

    // // Pressure applied to the element itself
    // double pressure_on_condition = 0.0;
    // if (TDim == 2) {
    //     if( this->Has( PRESSURE ) ) {
    //         pressure_on_condition += this->GetValue( PRESSURE );
    //     }
    //     if( this->Has( NEGATIVE_FACE_PRESSURE ) ) {
    //         pressure_on_condition += this->GetValue( NEGATIVE_FACE_PRESSURE );
    //     }
    //     if( this->Has( POSITIVE_FACE_PRESSURE ) ) {
    //         pressure_on_condition -= this->GetValue( POSITIVE_FACE_PRESSURE );
    //     }

    //     for ( IndexType i = 0; i < pressure_on_nodes.size(); i++ ) {
    //         pressure_on_nodes[i] = pressure_on_condition;
    //         if( r_geometry[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) ) {
    //             pressure_on_nodes[i] += r_geometry[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
    //         }
    //         if( r_geometry[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) ) {
    //             pressure_on_nodes[i] -= r_geometry[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
    //         }
    //     }
    // }

    // // Vector with a loading applied to the elemnt
    // array_1d<double, 3 > line_load = ZeroVector(3);
    // if( this->Has( LINE_LOAD ) ) {
    //     noalias(line_load) = this->GetValue( LINE_LOAD );
    // }

    // // Declaring tangent and Jacobian
    // array_1d<double, 3> tangent_xi, tangent_eta;
    // Matrix J(TDim, 1);

    // // Iterate over the Gauss points
    // for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
    //     const double det_j = r_geometry.DeterminantOfJacobian( integration_points[point_number] );
    //     const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);

    //     // Calculating the pressure on the gauss point
    //     double gauss_pressure = 0.0;
    //     for ( IndexType ii = 0; ii < number_of_nodes; ii++ ) {
    //         gauss_pressure += rNcontainer( point_number, ii ) * pressure_on_nodes[ii];
    //     }

    //     // Definition of the tangent
    //     if ( gauss_pressure != 0.0 ) {
    //         // Definition of the tangent
    //         r_geometry.Jacobian(J, point_number, integration_method);
    //         GetLocalAxis1(tangent_xi, J);
    //     }

    //     // Adding contributions to the LHS matrix
    //     if ( CalculateStiffnessMatrixFlag ) {
    //         if ( gauss_pressure != 0.0 ) {
    //             CalculateAndSubKp( rLeftHandSideMatrix, tangent_xi, rDN_De[point_number], row( rNcontainer, point_number ), gauss_pressure, integration_weight );
    //         }
    //     }
    //     // Adding contributions to the residual vector
    //     if ( CalculateResidualVectorFlag ) {
    //         if ( gauss_pressure != 0.0 ) {
    //             array_1d<double, 3> normal;

    //             // Getting LOCAL_AXIS_2
    //             GetLocalAxis2(tangent_eta);

    //             // Computing normal
    //             MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);

    //             CalculateAndAddPressureForce( rRightHandSideVector, row( rNcontainer, point_number ), normal, gauss_pressure, integration_weight );
    //         }
    //     }

    //     array_1d<double,3> gauss_load = line_load;
    //     for (IndexType ii = 0; ii < number_of_nodes; ++ii) {
    //         if( r_geometry[ii].SolutionStepsDataHas( LINE_LOAD ) ) {
    //             noalias(gauss_load) += ( rNcontainer( point_number, ii )) * r_geometry[ii].FastGetSolutionStepValue( LINE_LOAD );
    //         } else if( r_geometry[ii].Has( LINE_LOAD ) ) {
    //             noalias(gauss_load) += ( rNcontainer( point_number, ii )) * r_geometry[ii].GetValue( LINE_LOAD );
    //         }
    //     }

    //     for (IndexType ii = 0; ii < number_of_nodes; ++ii) {
    //         const IndexType base = ii * block_size;

    //         for (IndexType k = 0; k < TDim; ++k) {
    //             rRightHandSideVector[base + k] += integration_weight * rNcontainer( point_number, ii ) * gauss_load[k];
    //         }
    //     }
    // }


    // CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::CalculateMassMatrix(
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
    // if(rMassMatrix.size1() != 0) {
    //     rMassMatrix.resize(0, 0, false);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    if(rDampingMatrix.size1() != 0) {
        rDampingMatrix.resize(0, 0, false);
    }
}


template<std::size_t TDim>
void AcousticStructureLineCondition<TDim>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateMassMatrixFlag,
    // const bool CalculateStiffnessMatrixFlag,
    const bool CalculateVectorFlag
    )
{
    KRATOS_TRY;
    std::cout << "damn alright!\n";
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
    // KRATOS_WATCH(mat_size)

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
        // const GeometryType::ShapeFunctionsGradientsType& rDN_De = r_geometry.ShapeFunctionsLocalGradients(integration_method);
        const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

        // KRATOS_WATCH(integration_points)
        // KRATOS_WATCH(rDN_De)
        KRATOS_WATCH(dimension)

        Matrix J0(dimension, dimension);
        Matrix J(dimension, 1);
        array_1d<double, 3> tangent_xi, tangent_eta;
        array_1d<double, 3> normal;

        int sign;
        CalculateMassMatrixFlag ? sign = -1 : sign = 1;

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            r_geometry.Jacobian(J, point_number, integration_method);
            GetLocalAxis1(tangent_xi, J);
            GetLocalAxis2(tangent_eta);
            MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);

            // KRATOS_WATCH(normal)
            // GeometryUtils::JacobianOnInitialConfiguration(
            //     r_geometry, integration_points[point_number], J0);
            // KRATOS_WATCH(J0)
            // KRATOS_WATCH(J1)
            const double detJ = r_geometry.DeterminantOfJacobian( integration_points[point_number] );
            // const double detJ = MathUtils<double>::DetMat(J);
            // KRATOS_WATCH(J)
            // KRATOS_WATCH(detJ)
            const double integration_weight =
                GetIntegrationWeight(integration_points, point_number, detJ);
            const Vector& rN = row(Ncontainer,point_number);

            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                const SizeType index_i = i * dimension;

                for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                    // const SizeType index_j = j * dimension + dimension; //j +
                    const SizeType index_j = j + (dimension * number_of_nodes); //j +
                    const double NiNj_weight = rN[i] * rN[j] * integration_weight;
                    for ( IndexType k = 0; k < dimension; ++k ) {
                        // KRATOS_WATCH(index_i+k)
                        // KRATOS_WATCH(index_j+k)
                        rLeftHandSideMatrix( index_i + k, index_j ) += NiNj_weight * normal(k) * sign;
                    }
                }
            }

        }
        std::cout << "yo?\n";
        if( CalculateMassMatrixFlag ) {
            // rLeftHandSideMatrix *= -1.0;
            rLeftHandSideMatrix = trans(rLeftHandSideMatrix);
        }

        // KRATOS_WATCH(CalculateStiffnessMatrixFlag)
        // KRATOS_WATCH(CalculateMassMatrixFlag)
        KRATOS_WATCH(CalculateMassMatrixFlag)
        KRATOS_WATCH(this->Id())
        KRATOS_WATCH(rLeftHandSideMatrix)
    }
    std::cout << "condition end\n";

            // MATLAB STYLE
            //     for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            //     const SizeType index_i = i * dimension;

            //     for ( IndexType j = 0; j < number_of_nodes; ++j ) {
            //         // const SizeType index_j = j * dimension; //j +
            //         const SizeType index_j = j + (dimension * number_of_nodes); //j +
            //         // const SizeType index_j = (dimension * number_of_nodes);
            //         const double NiNj_weight = rN[i] * rN[j] * integration_weight;
            //         for ( IndexType k = 0; k < dimension; ++k ) {
            //             // KRATOS_WATCH(index_i+k)
            //             // KRATOS_WATCH(index_j+k)
            //             rLeftHandSideMatrix( index_i + k, index_j ) += NiNj_weight * normal(k) * sign;
            //         }
            //     }
            // }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

// template<std::size_t TDim>
// void AcousticStructureLineCondition<TDim>::CalculateAndSubKp(
//     Matrix& rK,
//     const array_1d<double, 3>& rTangentXi,
//     const Matrix& rDN_De,
//     const Vector& rN,
//     const double Pressure,
//     const double IntegrationWeight
//     ) const
// {
//     KRATOS_TRY

//     // Getting geometry
//     const auto& r_geometry = GetGeometry();

//     // Local stiffness
//     Matrix Kij( TDim, TDim );

//     // Getting cross tangent matrix
//     BoundedMatrix<double, TDim, TDim> cross_tangent_xi;
//     GetCrossTangentMatrix(cross_tangent_xi, rTangentXi);

//     // Getting geometry
//     const SizeType number_of_nodes = r_geometry.size();

//     for ( IndexType i = 0; i < number_of_nodes; i++ ) {
//         const IndexType row_index = i * TDim;

//         for ( IndexType j = 0; j < number_of_nodes; j++ ) {
//             const IndexType col_index = j * TDim;

//             const double coeff = Pressure * rN[i] * rDN_De( j, 0 ) * IntegrationWeight;
//             noalias(Kij) = coeff * cross_tangent_xi;

//             MathUtils<double>::AddMatrix( rK, Kij, row_index, col_index );
//         }
//     }

//     KRATOS_CATCH( "" )
// }

/***********************************************************************************/
/***********************************************************************************/

// template<std::size_t TDim>
// void AcousticStructureLineCondition<TDim>::CalculateAndAddPressureForce(
//     Vector& rRightHandSideVector,
//     const Vector& rN,
//     const array_1d<double, 3>& rNormal,
//     double Pressure,
//     double IntegrationWeight
//     ) const
// {
//     const SizeType number_of_nodes = this->GetGeometry().size();
//     const SizeType block_size = this->GetBlockSize();

//     for ( IndexType i = 0; i < number_of_nodes; ++i ) {
//         const IndexType index = block_size * i;

//         const double coeff = Pressure * rN[i] * IntegrationWeight;

//         for ( IndexType j = 0; j < TDim; ++j ) {
//             rRightHandSideVector[index + j] -= coeff * rNormal[j];
//         }
//     }
// }

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureLineCondition<2>::GetLocalAxis1(
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
void AcousticStructureLineCondition<3>::GetLocalAxis1(
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
void AcousticStructureLineCondition<2>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis) const
{
    rLocalAxis[0] = 0.0;
    rLocalAxis[1] = 0.0;
    rLocalAxis[2] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureLineCondition<3>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis) const
{
    KRATOS_ERROR_IF(!Has(LOCAL_AXIS_2)) << "The variable LOCAL_AXIS_2 is needed to compute the normal" << std::endl;
    noalias(rLocalAxis) = this->GetValue(LOCAL_AXIS_2);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureLineCondition<2>::GetCrossTangentMatrix(
    BoundedMatrix<double, 2, 2>& rCrossTangentMatrix,
    const array_1d<double, 3>& rTangentXi
    ) const
{
    const auto& r_properties = GetProperties();
    const double h0 = r_properties.Has(THICKNESS) ? r_properties.GetValue(THICKNESS): 1.0;
    rCrossTangentMatrix( 0, 0 ) =  0.0;
    rCrossTangentMatrix( 0, 1 ) =  h0;
    rCrossTangentMatrix( 1, 0 ) = -h0;
    rCrossTangentMatrix( 1, 1 ) =  0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureLineCondition<3>::GetCrossTangentMatrix(
    BoundedMatrix<double, 3, 3>& rCrossTangentMatrix,
    const array_1d<double, 3>& rTangentXi
    ) const
{
    BeamMathUtils<double>::VectorToSkewSymmetricTensor(rTangentXi, rCrossTangentMatrix);
}

/***********************************************************************************/
/***********************************************************************************/

template class AcousticStructureLineCondition<2>;
template class AcousticStructureLineCondition<3>;

} // Namespace Kratos


