// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/geometrical_projection_utilities.h"
/* Mortar includes */
#include "custom_conditions/mesh_tying_mortar_condition.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster> >( NewId, pGeom, pProperties, pMasterGeom);
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::~MeshTyingMortarCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    // We get the unkown variable
    const std::string variable_name = GetProperties().Has(TYING_VARIABLE) ? GetProperties().GetValue(TYING_VARIABLE) : "DISPLACEMENT";
    if (KratosComponents<Variable<double>>::Has(variable_name)) {
        mDoubleVariables.push_back(KratosComponents<Variable<double>>::Get(variable_name));
    } else if (KratosComponents<Array1DComponentsType>::Has(variable_name)) {
        mArray1DVariables.push_back(KratosComponents<Array1DComponentsType>::Get(variable_name));
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name)) {
        mArray1DVariables.push_back(KratosComponents<Array1DComponentsType>::Get(variable_name+"_X"));
        mArray1DVariables.push_back(KratosComponents<Array1DComponentsType>::Get(variable_name+"_Y"));
        if (TDim == 3) mArray1DVariables.push_back(KratosComponents<Array1DComponentsType>::Get(variable_name+"_Z"));
    } else {
        KRATOS_ERROR << "Compatible variables are: double, array_1d<double, 3> or components of the former" << std::endl;
    }

    // We define the integration method
    mIntegrationOrder = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;

    // The slave geometry
    GeometryType& r_slave_geometry = this->GetGeometry();
    const array_1d<double, 3>& r_normal_slave = this->GetValue(NORMAL);

    // Create and initialize condition variables:
    GeneralVariables rVariables;

    // Create Ae matrix
    MatrixDualLM Ae;

    // The master geometry
    GeometryType& r_master_geometry = this->GetPairedGeometry();
    const array_1d<double, 3>& r_normal_master = this->GetValue(PAIRED_NORMAL);
    // Initialize general variables for the current master element
    rVariables.Initialize();

    // Initialize the mortar operators
    mrThisMortarConditionMatrices.Initialize();

    // We call the exact integration utility
    IntegrationUtility integration_utility = IntegrationUtility (mIntegrationOrder);

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(r_slave_geometry, r_normal_slave, r_master_geometry, r_normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    if ((is_inside == true) && ((integration_area/r_slave_geometry.Area()) > 1.0e-3 * r_slave_geometry.Area())) {
        IntegrationMethod this_integration_method = GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Initialize the mortar operators
        mrThisMortarConditionMatrices.Initialize();

        const bool dual_LM = CalculateAe(r_normal_master, Ae, rVariables, conditions_points_slave, this_integration_method);

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates
                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    this->CalculateKinematics( rVariables, Ae, r_normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double integration_weight = integration_points_slave[point_number].Weight();

                    mrThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }
    } else { // We deactivate
        this->Set(ACTIVE, false);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // TODO: Add things if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Compute the matrix size
    const TensorValue tensor_value = (mDoubleVariables.size() == 1) ? ScalarValue : static_cast<TensorValue>(TDim);
    const SizeType matrix_size = tensor_value * (2 * NumNodes + NumNodesMaster);

    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size )
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );

    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != matrix_size )
        rRightHandSideVector.resize( matrix_size, false );

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // Compute the matrix size
    const TensorValue tensor_value = (mDoubleVariables.size() == 1) ? ScalarValue : static_cast<TensorValue>(TDim);
    const SizeType matrix_size = tensor_value * (2 * NumNodes + NumNodesMaster);

    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size )
        rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );

    // Creating an auxiliar vector
    VectorType aux_right_hand_side_vector = Vector();

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo, true, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // Creating an auxiliar matrix
    MatrixType aux_left_hand_side_matrix = Matrix();

    // Compute the matrix size
    const TensorValue tensor_value = (mDoubleVariables.size() == 1) ? ScalarValue : static_cast<TensorValue>(TDim);
    const SizeType matrix_size = tensor_value * (2 * NumNodes + NumNodesMaster);

    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != matrix_size )
        rRightHandSideVector.resize( matrix_size, false );

    // Calculate condition system
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim, TNumNodesElem, TNumNodesElemMaster>::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY;

    // Create the current DoF data
    const TensorValue tensor_value = (mDoubleVariables.size() == 1) ? ScalarValue : static_cast<TensorValue>(TDim);

    if (tensor_value ==  ScalarValue) {
        DofData<ScalarValue> dof_data;

        // Initialize the DoF data
        this->InitializeDofData<ScalarValue>(dof_data);

        // Update slave element info
        dof_data.UpdateMasterPair(this->GetPairedGeometry(), mDoubleVariables, mArray1DVariables);

        // Assemble of the matrix is required
        if ( ComputeLHS ) {
            // Calculate the local contribution
            this->CalculateLocalLHS<ScalarValue>(rLeftHandSideMatrix, mrThisMortarConditionMatrices, dof_data);
        }

        // Assemble of the vector is required
        if ( ComputeRHS) {
            // Calculate the local contribution
            this->CalculateLocalRHS<ScalarValue>( rRightHandSideVector, mrThisMortarConditionMatrices, dof_data);
        }
    } else {
        DofData<static_cast<TensorValue>(TDim)> dof_data;

        // Initialize the DoF data
        this->InitializeDofData<static_cast<TensorValue>(TDim)>(dof_data);

        // Update slave element info
        dof_data.UpdateMasterPair(this->GetPairedGeometry(), mDoubleVariables, mArray1DVariables);

        // Assemble of the matrix is required
        if ( ComputeLHS ) {
            // Calculate the local contribution
            this->CalculateLocalLHS<static_cast<TensorValue>(TDim)>(rLeftHandSideMatrix, mrThisMortarConditionMatrices, dof_data);
        }

        // Assemble of the vector is required
        if ( ComputeRHS) {
            // Calculate the local contribution
            this->CalculateLocalRHS<static_cast<TensorValue>(TDim)>( rRightHandSideVector, mrThisMortarConditionMatrices, dof_data);
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
bool MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateAe(
    const array_1d<double, 3>& rNormalMaster,
    MatrixDualLM& rAe,
    GeneralVariables& rVariables,
    ConditionArrayListType& rConditionsPointsSlave,
    IntegrationMethod ThisIntegrationMethod
    )
{
    // We initilize the Ae components
    AeData rAeData;
    rAeData.Initialize();

    rAe = ZeroMatrix(NumNodes, NumNodes);

    // The slave geometry
    GeometryType& r_slave_geometry = this->GetGeometry();

    // Initialize general variables for the current master element
    rVariables.Initialize();

    // Calculating the proportion between the integrated area and segment area
    for (IndexType i_geom = 0; i_geom < rConditionsPointsSlave.size(); ++i_geom) {
        PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            r_slave_geometry.GlobalCoordinates(global_point, rConditionsPointsSlave[i_geom][i_node]);
            points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
        }

        DecompositionType decomp_geom( points_array );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

        if (bad_shape == false) {
            const GeometryType::IntegrationPointsArrayType integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

            // Integrating the mortar operators
            for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                // We compute the local coordinates
                const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                PointType local_point_parent;
                PointType gp_global;
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                // Calculate the kinematic variables
                // We compute the current configuration
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                this->CalculateKinematics( rVariables, rAe, rNormalMaster, local_point_decomp, local_point_parent, decomp_geom, false);

                // Integrate
                const double integration_weight = integration_points_slave[point_number].Weight();

                rAeData.CalculateAeComponents(rVariables, integration_weight);
            }
        }
    }

    return rAeData.CalculateAe(rAe);
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateKinematics(
    GeneralVariables& rVariables,
    const MatrixDualLM& rAe,
    const array_1d<double, 3>& rNormalMaster,
    const PointType& rLocalPointDecomp,
    const PointType& rLocalPointParent,
    GeometryPointType& rGeometryDecomp,
    const bool DualLM
    )
{
    /// SLAVE CONDITION ///
    /* SHAPE FUNCTIONS */
    GetGeometry().ShapeFunctionsValues( rVariables.NSlave, rLocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = (DualLM == true) ? prod(rAe, rVariables.NSlave) : rVariables.NSlave;

    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.DetjSlave = rGeometryDecomp.DeterminantOfJacobian( rLocalPointDecomp );

    KRATOS_ERROR_IF(rVariables.DetjSlave < 0.0) << "WARNING:: CONDITION ID: " << this->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;

    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, rNormalMaster, rLocalPointParent);
}

/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const array_1d<double, 3>& rNormalMaster,
    const PointType& rLocalPoint
    )
{
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    PointType projected_gp_global;
    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, GetGeometry());

    GeometryType::CoordinatesArrayType slave_gp_global;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, rLocalPoint );
    GeometricalProjectionUtilities::FastProjectDirection( r_master_geometry, slave_gp_global, projected_gp_global, rNormalMaster, -gp_normal ); // The opposite direction

    GeometryType::CoordinatesArrayType projected_gp_local;

    r_master_geometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;

    // SHAPE FUNCTIONS
    r_master_geometry.ShapeFunctionsValues( rVariables.NMaster, projected_gp_local );
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_lhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_rhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Compute the matrix size
    const TensorValue tensor_value = (mDoubleVariables.size() == 1) ? ScalarValue : static_cast<TensorValue>(TDim);
    const SizeType matrix_size = tensor_value * (2 * NumNodes + NumNodesMaster);

    if (rResult.size() != matrix_size) {
        rResult.resize( matrix_size, false );
    }

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE, LM ] */
    // Master Nodes DoF Equation IDs
    GeometryType& r_current_master = this->GetPairedGeometry();

    if (tensor_value == ScalarValue) {
        for ( IndexType i_master = 0; i_master < NumNodesMaster; ++i_master ) {
            rResult[index++] = r_current_master[i_master].GetDof( mDoubleVariables[0] ).EquationId( );
        }
    } else {
        for ( IndexType i_master = 0; i_master < NumNodesMaster; ++i_master ) {
            NodeType& r_master_node = r_current_master[i_master];
            for (IndexType i_dof = 0; i_dof < mArray1DVariables.size(); i_dof++)
                rResult[index++] = r_master_node.GetDof( mArray1DVariables[i_dof] ).EquationId( );
        }
    }

    // Slave Nodes DoF Equation IDs
    if (tensor_value == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            rResult[index++] = this->GetGeometry()[i_slave].GetDof( mDoubleVariables[0] ).EquationId( );
        }
    } else {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            for (IndexType i_dof = 0; i_dof < mArray1DVariables.size(); i_dof++)
                rResult[index++] = slave_node.GetDof( mArray1DVariables[i_dof] ).EquationId( );
        }
    }

    // Slave Nodes LM Equation IDs
    if (tensor_value == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave )  {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rResult[index++] = slave_node.GetDof( SCALAR_LAGRANGE_MULTIPLIER ).EquationId( );
        }
    } else {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
            if (mArray1DVariables.size() > 1) rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
            if (mArray1DVariables.size() > 2) rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim, TNumNodesElem, TNumNodesElemMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    const TensorValue tensor_value = (mDoubleVariables.size() == 1) ? ScalarValue : static_cast<TensorValue>(TDim);
    const SizeType matrix_size = tensor_value * (2 * NumNodes + NumNodesMaster);

    if (rConditionalDofList.size() != matrix_size) {
        rConditionalDofList.resize( matrix_size );
    }

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE, LM ] */
    // Master Nodes DoF Equation IDs
    GeometryType& current_master = this->GetPairedGeometry();

    if (tensor_value == ScalarValue) {
        for ( IndexType i_master = 0; i_master < NumNodesMaster; ++i_master )  {
            rConditionalDofList[index++] = current_master[i_master].pGetDof( mDoubleVariables[0] );
        }
    } else {
        for ( IndexType i_master = 0; i_master < NumNodesMaster; ++i_master ) {
            NodeType& master_node = current_master[i_master];
            for (IndexType i_dof = 0; i_dof < mArray1DVariables.size(); i_dof++)
                rConditionalDofList[index++] = master_node.pGetDof( mArray1DVariables[i_dof] );
        }
    }

    // Slave Nodes DoF Equation IDs
    if (tensor_value == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            rConditionalDofList[index++] = this->GetGeometry()[i_slave].pGetDof( mDoubleVariables[0] );
        }
    } else {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            for (IndexType i_dof = 0; i_dof < mArray1DVariables.size(); i_dof++)
                rConditionalDofList[index++] = slave_node.pGetDof( mArray1DVariables[i_dof] );
        }
    }

    // Slave Nodes LM Equation IDs
    if (tensor_value == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rConditionalDofList[index++] = slave_node.pGetDof( SCALAR_LAGRANGE_MULTIPLIER );
        }
    } else {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
            if (mArray1DVariables.size() > 1) rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
            if (mArray1DVariables.size() > 2) rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
        }
    }

    KRATOS_CATCH( "" );
}


//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints();

    if ( rOutput.size() != r_integration_points.size() ) {
        rOutput.resize( r_integration_points.size() );
    }

    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        rOutput[point_number] = 0.0;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints();

    if ( rOutput.size() != r_integration_points.size() ) {
        rOutput.resize( r_integration_points.size() );
    }

    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        rOutput[point_number] = ZeroVector(3);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
void MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints();

    if ( rOutput.size() != r_integration_points.size() ) {
        rOutput.resize( r_integration_points.size() );
    }

    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        rOutput[point_number] = ZeroVector(3);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodesElem, SizeType TNumNodesElemMaster>
int MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Condition::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    KRATOS_ERROR_IF(BaseType::mpPairedGeometry == nullptr) << "YOU HAVE NOT INITIALIZED THE PAIR GEOMETRY IN THE MeshTyingMortarCondition" << std::endl;

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class MeshTyingMortarCondition<2, 3, 3>; // 2DLine/Triangle
template class MeshTyingMortarCondition<2, 4, 4>; // 2DLine/Quadrilateral
template class MeshTyingMortarCondition<3, 4, 4>; // 3D Triangle/Tetrahedron
template class MeshTyingMortarCondition<3, 8, 8>; // 3D Quadrilateral/Hexahedra
template class MeshTyingMortarCondition<3, 4, 8>; // 3D Triangle/Tetrahedron
template class MeshTyingMortarCondition<3, 8, 4>; // 3D Quadrilateral/Hexahedra

} // Namespace Kratos
