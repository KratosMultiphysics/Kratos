//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
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

template<std::size_t TDim, bool TIsMapping>
AcousticStructureCouplingCondition<TDim, TIsMapping>::AcousticStructureCouplingCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
AcousticStructureCouplingCondition<TDim, TIsMapping>::AcousticStructureCouplingCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
Condition::Pointer AcousticStructureCouplingCondition<TDim, TIsMapping>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticStructureCouplingCondition<TDim, TIsMapping>>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
Condition::Pointer AcousticStructureCouplingCondition<TDim, TIsMapping>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AcousticStructureCouplingCondition<TDim, TIsMapping>>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
Condition::Pointer AcousticStructureCouplingCondition<TDim, TIsMapping>::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<AcousticStructureCouplingCondition<TDim, TIsMapping>>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
AcousticStructureCouplingCondition<TDim, TIsMapping>::~AcousticStructureCouplingCondition()
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
template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    //if has info about mapped eq ids provide this
    //..
    //..


    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType block_size = this->GetBlockSize();

    if( !TIsMapping ) {
        if (rResult.size() != block_size * number_of_nodes) {
            rResult.resize(number_of_nodes * block_size, false);
        }

        const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
        const SizeType pos_p = this->GetGeometry()[0].GetDofPosition(PRESSURE);

        if(dim == 2) {
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
                //this was not correct:
                // const SizeType index = i * (block_size - 1);

                //this should be fine for disp(xyz) - pressure - rot(xyz)
                const SizeType index = i * dim;

                rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
                rResult[i + dim * number_of_nodes] = GetGeometry()[i].GetDof(PRESSURE,pos_p).EquationId();
            }
        }
    } else {
        // std::cout << "i am a mapping element yay\n";
        // KRATOS_WATCH(*this)
        auto& mapped_nodes = this->GetValue(MAPPING_NODES);
        const SizeType n_mapped_nodes = mapped_nodes.size();

        KRATOS_ERROR_IF( n_mapped_nodes < 1 ) << "No mapping nodes defined!" << std::endl;



        if( this->GetGeometry()[0].HasDofFor(DISPLACEMENT_X) ) {
            //mapping from coarse structure to fine fluid


        } else {
            //mapping from coarse fluid to fine structure
            // block size - 1 = disp/rot dofs * number of mapped nodes + 1 dof per node here
            // if (rResult.size() != n_mapped_nodes * (this->GetBlockSize() - 1) + number_of_nodes) {
            //     rResult.resize(n_mapped_nodes * (this->GetBlockSize() - 1) + number_of_nodes, false);
            // }
            rResult.resize(0);
            rResult.reserve(this->GetSystemSize());
            // rResult.reserve(n_mapped_nodes * (this->GetBlockSize() - 1) + number_of_nodes);
            if(dim == 2) {
                for (SizeType i = 0; i < n_mapped_nodes; ++i) {
                    // const SizeType index = this->HasRotDof() ? i * (block_size - 2) : i * (block_size - 1);

                    rResult.push_back(mapped_nodes[i]->GetDof(DISPLACEMENT_X).EquationId());
                    rResult.push_back(mapped_nodes[i]->GetDof(DISPLACEMENT_Y).EquationId());
                    // rResult[i + dim * number_of_nodes] = mapped_nodes[i].GetDof(PRESSURE,pos_p).EquationId();

                }
                for( SizeType i=0; i<number_of_nodes; ++i ) {
                    rResult.push_back(GetGeometry()[i].GetDof(PRESSURE).EquationId());
                }
                if ( this-> HasRotDof() ) {
                    // std::cout << "MAPPED NODES HAVE ROTATION DOF!!!\n";
                    // KRATOS_WATCH(rResult)
                    for( SizeType i=0; i<n_mapped_nodes; ++i ) {
                        rResult.push_back(mapped_nodes[i]->GetDof(ROTATION_Z).EquationId());
                    }
                }
            } else {
                for (SizeType i = 0; i < n_mapped_nodes; ++i) {
                    // const SizeType index = this->HasRotDof() ? i * (block_size - 2) : i * (block_size - 1);

                    rResult.push_back(mapped_nodes[i]->GetDof(DISPLACEMENT_X).EquationId());
                    rResult.push_back(mapped_nodes[i]->GetDof(DISPLACEMENT_Y).EquationId());
                    rResult.push_back(mapped_nodes[i]->GetDof(DISPLACEMENT_Z).EquationId());
                    // rResult[i + dim * number_of_nodes] = mapped_nodes[i].GetDof(PRESSURE,pos_p).EquationId();

                }
                for( SizeType i=0; i<number_of_nodes; ++i ) {
                    rResult.push_back(GetGeometry()[i].GetDof(PRESSURE).EquationId());
                }
            }
        }
        // KRATOS_WATCH(rResult)

        // KRATOS_WATCH(this->GetValue(MAPPING_EQ_ID))
        // auto& r_mapping_eq_id = this->GetValue(MAPPING_EQ_ID);
        // const SizeType result_size = r_mapping_eq_id.size();
        // if (rResult.size() != result_size) {
        //     rResult.resize(result_size, false);
        // }

        // for( SizeType i=0; i<result_size; ++i ) {
        //     rResult[i] = static_cast<unsigned long>(r_mapping_eq_id[i]);
        // }
    }
    // KRATOS_WATCH(rResult)

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    // std::cout << "get dof list\n";

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim =  GetGeometry().WorkingSpaceDimension();
    const SizeType block_size = this->GetBlockSize();
    ElementalDofList.resize(0);

    if( !TIsMapping ) {
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
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof(PRESSURE));
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_ERROR << "Condition not prepared for time step analysis" << std::endl;
}

template<std::size_t TDim, bool TIsMapping>
int AcousticStructureCouplingCondition<TDim, TIsMapping>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    if( !TIsMapping ) {
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    }
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE)

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_node)
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node)

        if( !TIsMapping ) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
            if( TDim == 3)
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
        }
    }

    // Check required data for mapping
    if( TIsMapping ) {
        KRATOS_CHECK_VARIABLE_KEY(MAPPING_NODES)
        KRATOS_CHECK_VARIABLE_KEY(MAPPING_FACTOR)
        // KRATOS_WATCH(MAPPING_FACTOR)
        //TODO throw error for non-matching grids
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
double AcousticStructureCouplingCondition<TDim, TIsMapping>::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::CalculateRightHandSide(
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
template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::CalculateLocalSystem(
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

template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::CalculateMassMatrix(
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

template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType block_size = this->GetBlockSize();
    SizeType mat_size;

    if( !TIsMapping )
        mat_size = number_of_nodes * block_size;
    else
        mat_size = this->GetSystemSize();

    if( rDampingMatrix.size1() != mat_size || rDampingMatrix.size2() != mat_size ) {
        rDampingMatrix.resize(mat_size, mat_size, false);
    }

    noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);
}


template<std::size_t TDim, bool TIsMapping>
void AcousticStructureCouplingCondition<TDim, TIsMapping>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateMassMatrixFlag,
    const bool CalculateVectorFlag
    )
{
    KRATOS_TRY;
    // std::cout << "condition calculate all\n";
    // KRATOS_WATCH(*this)
    // KRATOS_WATCH(CalculateStiffnessMatrixFlag)
    // KRATOS_WATCH(CalculateVectorFlag)

    //if info about mapped stuff, call new function to build evtl rectangular lhs (or add this after the normal computation, should be easier)

    //else
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;
    const SizeType mapping_mat_size = this->GetSystemSize();
    // } else {
    //     // auto& mapped_nodes = this->GetValue(MAPPING_NODES);
    //     // const SizeType n_mapped_nodes = mapped_nodes.size();
    //     // mat_size = n_mapped_nodes * (this->GetBlockSize() - 1) + number_of_nodes;
    //     // if (mapped_nodes[0]->HasDofFor(ROTATION_Z)) {
    //     //     mat_size +=
    //     mat_size = this->GetSystemSize();
    //     KRATOS_WATCH(this->HasRotDof())
    //     KRATOS_WATCH(mat_size)
    // }
    // std::cout << "ok0\n";

    // if ( CalculateStiffnessMatrixFlag || CalculateMassMatrixFlag ) {
    //     if ( rLeftHandSideMatrix.size1() != mat_size ) {
    //         rLeftHandSideMatrix.resize( mat_size, mat_size, false );
    //     }
    //     noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    // }
    // KRATOS_WATCH(mat_size)
    // Resizing as needed the RHS or LHS
    if ( CalculateVectorFlag ) {
        if( !TIsMapping ) {
            if ( rRightHandSideVector.size( ) != mat_size ) {
                rRightHandSideVector.resize( mat_size, false );
            }
            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        } else {
            if ( rRightHandSideVector.size( ) != mapping_mat_size ) {
                rRightHandSideVector.resize( mapping_mat_size, false );
            }
            noalias( rRightHandSideVector ) = ZeroVector( mapping_mat_size ); //resetting RHS
        }
    }
    if ( CalculateStiffnessMatrixFlag || CalculateMassMatrixFlag ) {
        if ( rLeftHandSideMatrix.size1() != mat_size ) {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }
    // KRATOS_WATCH(rLeftHandSideMatrix)
    // std::cout << "ok1\n";

    if ( CalculateStiffnessMatrixFlag || CalculateMassMatrixFlag ) {
        // std::cout << "begin all\n";

        // Reading integration points and local gradients
        const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
        const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

        Matrix J(dimension, r_geometry.LocalSpaceDimension());
        array_1d<double, 3> tangent_xi, tangent_eta;
        array_1d<double, 3> normal;

        int sign = 1;
        if( CalculateMassMatrixFlag ) {
            sign = -1;
        } else if( TIsMapping && CalculateStiffnessMatrixFlag ) {
            sign = -1;
        }
        // CalculateMassMatrixFlag ? sign = -1 : sign = 1;

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
        // std::cout << "after all\n";

        // rLeftHandSideMatrix *= -1.0;

        // if( TIsMapping ) {
        //     //the normals of origin and destination must have the same direction
        //     // if( CalculateStiffnessMatrixFlag ) {
        //     //     rLeftHandSideMatrix *= -1.0;
        //     // }
        //     // rLeftHandSideMatrix *= -1.0;

        //     // //treat non-matching grids
        //     // if( mat_size != mapping_mat_size ) {

        //     //     //treat coarse structure to fine fluid!
        //     //     //TODO

        //     //     //coarse fluid to fine structure
        //     //     std::cout << "before:\n";
        //     //     KRATOS_WATCH(rLeftHandSideMatrix)
        //     //     MatrixType tmp_matrix = ZeroMatrix(mapping_mat_size,mapping_mat_size);
        //     //     KRATOS_WATCH(this->GetValue(MAPPING_FACTOR))




        //     //     rLeftHandSideMatrix = tmp_matrix;
        //     //     std::cout << "after:\n";
        //     //     KRATOS_WATCH(rLeftHandSideMatrix)
        //     // }
        // }

        if( CalculateMassMatrixFlag ) {
            rLeftHandSideMatrix = trans(rLeftHandSideMatrix);
            // if( !TIsMapping ) {
            //     rLeftHandSideMatrix *= -1.0;
            // }
            // KRATOS_WATCH(rLeftHandSideMatrix)
        } else {
            // if( TIsMapping ) {
            //     rLeftHandSideMatrix *= -1.0;
            // }
        }
    }
    // std::cout << "ok2\n";

    // if( TIsMapping ) {
    //     KRATOS_WATCH(rLeftHandSideMatrix)
    //     const SizeType system_size = this->GetSystemSize();
    //     MatrixType tmp_matrix = ZeroMatrix(system_size,system_size);
    //     rLeftHandSideMatrix = tmp_matrix;
    //     // rLeftHandSideMatrix.resize(1,1,false);
    //     // rLeftHandSideMatrix[0,0] = 1;
    //     KRATOS_WATCH(rLeftHandSideMatrix)
    //     // rRightHandSideVector.resize( system_size, false );
    //     // noalias( rRightHandSideVector ) = ZeroVector( system_size ); //resetting RHS
    //     // KRATOS_WATCH(rRightHandSideVector)
    // }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<2, false>::GetLocalAxis1(
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
void AcousticStructureCouplingCondition<3, false>::GetLocalAxis1(
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
void AcousticStructureCouplingCondition<2, false>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis, const Matrix& rJacobian) const
{
    rLocalAxis[0] = 0.0;
    rLocalAxis[1] = 0.0;
    rLocalAxis[2] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<3, false>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis, const Matrix& rJacobian) const
{
    rLocalAxis[0] = rJacobian(0, 1);
    rLocalAxis[1] = rJacobian(1, 1);
    rLocalAxis[2] = rJacobian(2, 1);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<2, true>::GetLocalAxis1(
    array_1d<double, 3>& rLocalAxis,
    const Matrix& rJacobian
    ) const
{
    rLocalAxis[0] = -rJacobian(0, 0);
    rLocalAxis[1] = -rJacobian(1, 0);
    rLocalAxis[2] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<3, true>::GetLocalAxis1(
    array_1d<double, 3>& rLocalAxis,
    const Matrix& rJacobian
    ) const
{
    rLocalAxis[0] = -rJacobian(0, 0);
    rLocalAxis[1] = -rJacobian(1, 0);
    rLocalAxis[2] = -rJacobian(2, 0);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<2, true>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis, const Matrix& rJacobian) const
{
    rLocalAxis[0] = 0.0;
    rLocalAxis[1] = 0.0;
    rLocalAxis[2] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AcousticStructureCouplingCondition<3, true>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis, const Matrix& rJacobian) const
{
    rLocalAxis[0] = -rJacobian(0, 1);
    rLocalAxis[1] = -rJacobian(1, 1);
    rLocalAxis[2] = -rJacobian(2, 1);
}

/***********************************************************************************/
/***********************************************************************************/

template class AcousticStructureCouplingCondition<2, false>;
template class AcousticStructureCouplingCondition<3, false>;
template class AcousticStructureCouplingCondition<2, true>;
template class AcousticStructureCouplingCondition<3, true>;

} // Namespace Kratos


