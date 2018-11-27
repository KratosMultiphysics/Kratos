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
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name)) {
        mArray1DVariables.push_back(KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name));
    } else {
        KRATOS_ERROR << "Compatible variables are: double or array_1d<double, 3> " << std::endl;
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

template< >
template< >
void MeshTyingMortarCondition<2,3,3>::CalculateLocalLHS<MeshTyingMortarCondition<2,3,3>::ScalarValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=clhs0;
    rLocalLHS(0,5)=clhs1;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=clhs2;
    rLocalLHS(1,5)=clhs3;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=DOperator(0,0);
    rLocalLHS(2,5)=DOperator(1,0);
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=DOperator(0,1);
    rLocalLHS(3,5)=DOperator(1,1);
    rLocalLHS(4,0)=clhs0;
    rLocalLHS(4,1)=clhs2;
    rLocalLHS(4,2)=DOperator(0,0);
    rLocalLHS(4,3)=DOperator(0,1);
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(5,0)=clhs1;
    rLocalLHS(5,1)=clhs3;
    rLocalLHS(5,2)=DOperator(1,0);
    rLocalLHS(5,3)=DOperator(1,1);
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<2,3,3>::CalculateLocalLHS<MeshTyingMortarCondition<2,3,3>::Vector2DValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector2DValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=0;
    rLocalLHS(0,8)=clhs0;
    rLocalLHS(0,9)=0;
    rLocalLHS(0,10)=clhs1;
    rLocalLHS(0,11)=0;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=0;
    rLocalLHS(1,8)=0;
    rLocalLHS(1,9)=clhs0;
    rLocalLHS(1,10)=0;
    rLocalLHS(1,11)=clhs1;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=0;
    rLocalLHS(2,8)=clhs2;
    rLocalLHS(2,9)=0;
    rLocalLHS(2,10)=clhs3;
    rLocalLHS(2,11)=0;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=0;
    rLocalLHS(3,8)=0;
    rLocalLHS(3,9)=clhs2;
    rLocalLHS(3,10)=0;
    rLocalLHS(3,11)=clhs3;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=0;
    rLocalLHS(4,8)=DOperator(0,0);
    rLocalLHS(4,9)=0;
    rLocalLHS(4,10)=DOperator(1,0);
    rLocalLHS(4,11)=0;
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=0;
    rLocalLHS(5,8)=0;
    rLocalLHS(5,9)=DOperator(0,0);
    rLocalLHS(5,10)=0;
    rLocalLHS(5,11)=DOperator(1,0);
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=DOperator(0,1);
    rLocalLHS(6,9)=0;
    rLocalLHS(6,10)=DOperator(1,1);
    rLocalLHS(6,11)=0;
    rLocalLHS(7,0)=0;
    rLocalLHS(7,1)=0;
    rLocalLHS(7,2)=0;
    rLocalLHS(7,3)=0;
    rLocalLHS(7,4)=0;
    rLocalLHS(7,5)=0;
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=DOperator(0,1);
    rLocalLHS(7,10)=0;
    rLocalLHS(7,11)=DOperator(1,1);
    rLocalLHS(8,0)=clhs0;
    rLocalLHS(8,1)=0;
    rLocalLHS(8,2)=clhs2;
    rLocalLHS(8,3)=0;
    rLocalLHS(8,4)=DOperator(0,0);
    rLocalLHS(8,5)=0;
    rLocalLHS(8,6)=DOperator(0,1);
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(8,11)=0;
    rLocalLHS(9,0)=0;
    rLocalLHS(9,1)=clhs0;
    rLocalLHS(9,2)=0;
    rLocalLHS(9,3)=clhs2;
    rLocalLHS(9,4)=0;
    rLocalLHS(9,5)=DOperator(0,0);
    rLocalLHS(9,6)=0;
    rLocalLHS(9,7)=DOperator(0,1);
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(9,11)=0;
    rLocalLHS(10,0)=clhs1;
    rLocalLHS(10,1)=0;
    rLocalLHS(10,2)=clhs3;
    rLocalLHS(10,3)=0;
    rLocalLHS(10,4)=DOperator(1,0);
    rLocalLHS(10,5)=0;
    rLocalLHS(10,6)=DOperator(1,1);
    rLocalLHS(10,7)=0;
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;
    rLocalLHS(10,11)=0;
    rLocalLHS(11,0)=0;
    rLocalLHS(11,1)=clhs1;
    rLocalLHS(11,2)=0;
    rLocalLHS(11,3)=clhs3;
    rLocalLHS(11,4)=0;
    rLocalLHS(11,5)=DOperator(1,0);
    rLocalLHS(11,6)=0;
    rLocalLHS(11,7)=DOperator(1,1);
    rLocalLHS(11,8)=0;
    rLocalLHS(11,9)=0;
    rLocalLHS(11,10)=0;
    rLocalLHS(11,11)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<2,4,4>::CalculateLocalLHS<MeshTyingMortarCondition<2,4,4>::ScalarValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=clhs0;
    rLocalLHS(0,5)=clhs1;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=clhs2;
    rLocalLHS(1,5)=clhs3;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=DOperator(0,0);
    rLocalLHS(2,5)=DOperator(1,0);
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=DOperator(0,1);
    rLocalLHS(3,5)=DOperator(1,1);
    rLocalLHS(4,0)=clhs0;
    rLocalLHS(4,1)=clhs2;
    rLocalLHS(4,2)=DOperator(0,0);
    rLocalLHS(4,3)=DOperator(0,1);
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(5,0)=clhs1;
    rLocalLHS(5,1)=clhs3;
    rLocalLHS(5,2)=DOperator(1,0);
    rLocalLHS(5,3)=DOperator(1,1);
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<2,4,4>::CalculateLocalLHS<MeshTyingMortarCondition<2,4,4>::Vector2DValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector2DValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=0;
    rLocalLHS(0,8)=clhs0;
    rLocalLHS(0,9)=0;
    rLocalLHS(0,10)=clhs1;
    rLocalLHS(0,11)=0;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=0;
    rLocalLHS(1,8)=0;
    rLocalLHS(1,9)=clhs0;
    rLocalLHS(1,10)=0;
    rLocalLHS(1,11)=clhs1;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=0;
    rLocalLHS(2,8)=clhs2;
    rLocalLHS(2,9)=0;
    rLocalLHS(2,10)=clhs3;
    rLocalLHS(2,11)=0;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=0;
    rLocalLHS(3,8)=0;
    rLocalLHS(3,9)=clhs2;
    rLocalLHS(3,10)=0;
    rLocalLHS(3,11)=clhs3;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=0;
    rLocalLHS(4,8)=DOperator(0,0);
    rLocalLHS(4,9)=0;
    rLocalLHS(4,10)=DOperator(1,0);
    rLocalLHS(4,11)=0;
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=0;
    rLocalLHS(5,8)=0;
    rLocalLHS(5,9)=DOperator(0,0);
    rLocalLHS(5,10)=0;
    rLocalLHS(5,11)=DOperator(1,0);
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=DOperator(0,1);
    rLocalLHS(6,9)=0;
    rLocalLHS(6,10)=DOperator(1,1);
    rLocalLHS(6,11)=0;
    rLocalLHS(7,0)=0;
    rLocalLHS(7,1)=0;
    rLocalLHS(7,2)=0;
    rLocalLHS(7,3)=0;
    rLocalLHS(7,4)=0;
    rLocalLHS(7,5)=0;
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=DOperator(0,1);
    rLocalLHS(7,10)=0;
    rLocalLHS(7,11)=DOperator(1,1);
    rLocalLHS(8,0)=clhs0;
    rLocalLHS(8,1)=0;
    rLocalLHS(8,2)=clhs2;
    rLocalLHS(8,3)=0;
    rLocalLHS(8,4)=DOperator(0,0);
    rLocalLHS(8,5)=0;
    rLocalLHS(8,6)=DOperator(0,1);
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(8,11)=0;
    rLocalLHS(9,0)=0;
    rLocalLHS(9,1)=clhs0;
    rLocalLHS(9,2)=0;
    rLocalLHS(9,3)=clhs2;
    rLocalLHS(9,4)=0;
    rLocalLHS(9,5)=DOperator(0,0);
    rLocalLHS(9,6)=0;
    rLocalLHS(9,7)=DOperator(0,1);
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(9,11)=0;
    rLocalLHS(10,0)=clhs1;
    rLocalLHS(10,1)=0;
    rLocalLHS(10,2)=clhs3;
    rLocalLHS(10,3)=0;
    rLocalLHS(10,4)=DOperator(1,0);
    rLocalLHS(10,5)=0;
    rLocalLHS(10,6)=DOperator(1,1);
    rLocalLHS(10,7)=0;
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;
    rLocalLHS(10,11)=0;
    rLocalLHS(11,0)=0;
    rLocalLHS(11,1)=clhs1;
    rLocalLHS(11,2)=0;
    rLocalLHS(11,3)=clhs3;
    rLocalLHS(11,4)=0;
    rLocalLHS(11,5)=DOperator(1,0);
    rLocalLHS(11,6)=0;
    rLocalLHS(11,7)=DOperator(1,1);
    rLocalLHS(11,8)=0;
    rLocalLHS(11,9)=0;
    rLocalLHS(11,10)=0;
    rLocalLHS(11,11)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,4,4>::CalculateLocalLHS<MeshTyingMortarCondition<3,4,4>::ScalarValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(0,1);
    const double clhs4 =     -MOperator(1,1);
    const double clhs5 =     -MOperator(2,1);
    const double clhs6 =     -MOperator(0,2);
    const double clhs7 =     -MOperator(1,2);
    const double clhs8 =     -MOperator(2,2);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=clhs0;
    rLocalLHS(0,7)=clhs1;
    rLocalLHS(0,8)=clhs2;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=clhs3;
    rLocalLHS(1,7)=clhs4;
    rLocalLHS(1,8)=clhs5;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=clhs6;
    rLocalLHS(2,7)=clhs7;
    rLocalLHS(2,8)=clhs8;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=DOperator(0,0);
    rLocalLHS(3,7)=DOperator(1,0);
    rLocalLHS(3,8)=DOperator(2,0);
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=DOperator(0,1);
    rLocalLHS(4,7)=DOperator(1,1);
    rLocalLHS(4,8)=DOperator(2,1);
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=DOperator(0,2);
    rLocalLHS(5,7)=DOperator(1,2);
    rLocalLHS(5,8)=DOperator(2,2);
    rLocalLHS(6,0)=clhs0;
    rLocalLHS(6,1)=clhs3;
    rLocalLHS(6,2)=clhs6;
    rLocalLHS(6,3)=DOperator(0,0);
    rLocalLHS(6,4)=DOperator(0,1);
    rLocalLHS(6,5)=DOperator(0,2);
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=0;
    rLocalLHS(7,0)=clhs1;
    rLocalLHS(7,1)=clhs4;
    rLocalLHS(7,2)=clhs7;
    rLocalLHS(7,3)=DOperator(1,0);
    rLocalLHS(7,4)=DOperator(1,1);
    rLocalLHS(7,5)=DOperator(1,2);
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(8,0)=clhs2;
    rLocalLHS(8,1)=clhs5;
    rLocalLHS(8,2)=clhs8;
    rLocalLHS(8,3)=DOperator(2,0);
    rLocalLHS(8,4)=DOperator(2,1);
    rLocalLHS(8,5)=DOperator(2,2);
    rLocalLHS(8,6)=0;
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,4,4>::CalculateLocalLHS<MeshTyingMortarCondition<3,4,4>::Vector3DValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(0,1);
    const double clhs4 =     -MOperator(1,1);
    const double clhs5 =     -MOperator(2,1);
    const double clhs6 =     -MOperator(0,2);
    const double clhs7 =     -MOperator(1,2);
    const double clhs8 =     -MOperator(2,2);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=0;
    rLocalLHS(0,8)=0;
    rLocalLHS(0,9)=0;
    rLocalLHS(0,10)=0;
    rLocalLHS(0,11)=0;
    rLocalLHS(0,12)=0;
    rLocalLHS(0,13)=0;
    rLocalLHS(0,14)=0;
    rLocalLHS(0,15)=0;
    rLocalLHS(0,16)=0;
    rLocalLHS(0,17)=0;
    rLocalLHS(0,18)=clhs0;
    rLocalLHS(0,19)=0;
    rLocalLHS(0,20)=0;
    rLocalLHS(0,21)=clhs1;
    rLocalLHS(0,22)=0;
    rLocalLHS(0,23)=0;
    rLocalLHS(0,24)=clhs2;
    rLocalLHS(0,25)=0;
    rLocalLHS(0,26)=0;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=0;
    rLocalLHS(1,8)=0;
    rLocalLHS(1,9)=0;
    rLocalLHS(1,10)=0;
    rLocalLHS(1,11)=0;
    rLocalLHS(1,12)=0;
    rLocalLHS(1,13)=0;
    rLocalLHS(1,14)=0;
    rLocalLHS(1,15)=0;
    rLocalLHS(1,16)=0;
    rLocalLHS(1,17)=0;
    rLocalLHS(1,18)=0;
    rLocalLHS(1,19)=clhs0;
    rLocalLHS(1,20)=0;
    rLocalLHS(1,21)=0;
    rLocalLHS(1,22)=clhs1;
    rLocalLHS(1,23)=0;
    rLocalLHS(1,24)=0;
    rLocalLHS(1,25)=clhs2;
    rLocalLHS(1,26)=0;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=0;
    rLocalLHS(2,8)=0;
    rLocalLHS(2,9)=0;
    rLocalLHS(2,10)=0;
    rLocalLHS(2,11)=0;
    rLocalLHS(2,12)=0;
    rLocalLHS(2,13)=0;
    rLocalLHS(2,14)=0;
    rLocalLHS(2,15)=0;
    rLocalLHS(2,16)=0;
    rLocalLHS(2,17)=0;
    rLocalLHS(2,18)=0;
    rLocalLHS(2,19)=0;
    rLocalLHS(2,20)=clhs0;
    rLocalLHS(2,21)=0;
    rLocalLHS(2,22)=0;
    rLocalLHS(2,23)=clhs1;
    rLocalLHS(2,24)=0;
    rLocalLHS(2,25)=0;
    rLocalLHS(2,26)=clhs2;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=0;
    rLocalLHS(3,8)=0;
    rLocalLHS(3,9)=0;
    rLocalLHS(3,10)=0;
    rLocalLHS(3,11)=0;
    rLocalLHS(3,12)=0;
    rLocalLHS(3,13)=0;
    rLocalLHS(3,14)=0;
    rLocalLHS(3,15)=0;
    rLocalLHS(3,16)=0;
    rLocalLHS(3,17)=0;
    rLocalLHS(3,18)=clhs3;
    rLocalLHS(3,19)=0;
    rLocalLHS(3,20)=0;
    rLocalLHS(3,21)=clhs4;
    rLocalLHS(3,22)=0;
    rLocalLHS(3,23)=0;
    rLocalLHS(3,24)=clhs5;
    rLocalLHS(3,25)=0;
    rLocalLHS(3,26)=0;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=0;
    rLocalLHS(4,8)=0;
    rLocalLHS(4,9)=0;
    rLocalLHS(4,10)=0;
    rLocalLHS(4,11)=0;
    rLocalLHS(4,12)=0;
    rLocalLHS(4,13)=0;
    rLocalLHS(4,14)=0;
    rLocalLHS(4,15)=0;
    rLocalLHS(4,16)=0;
    rLocalLHS(4,17)=0;
    rLocalLHS(4,18)=0;
    rLocalLHS(4,19)=clhs3;
    rLocalLHS(4,20)=0;
    rLocalLHS(4,21)=0;
    rLocalLHS(4,22)=clhs4;
    rLocalLHS(4,23)=0;
    rLocalLHS(4,24)=0;
    rLocalLHS(4,25)=clhs5;
    rLocalLHS(4,26)=0;
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=0;
    rLocalLHS(5,8)=0;
    rLocalLHS(5,9)=0;
    rLocalLHS(5,10)=0;
    rLocalLHS(5,11)=0;
    rLocalLHS(5,12)=0;
    rLocalLHS(5,13)=0;
    rLocalLHS(5,14)=0;
    rLocalLHS(5,15)=0;
    rLocalLHS(5,16)=0;
    rLocalLHS(5,17)=0;
    rLocalLHS(5,18)=0;
    rLocalLHS(5,19)=0;
    rLocalLHS(5,20)=clhs3;
    rLocalLHS(5,21)=0;
    rLocalLHS(5,22)=0;
    rLocalLHS(5,23)=clhs4;
    rLocalLHS(5,24)=0;
    rLocalLHS(5,25)=0;
    rLocalLHS(5,26)=clhs5;
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=0;
    rLocalLHS(6,9)=0;
    rLocalLHS(6,10)=0;
    rLocalLHS(6,11)=0;
    rLocalLHS(6,12)=0;
    rLocalLHS(6,13)=0;
    rLocalLHS(6,14)=0;
    rLocalLHS(6,15)=0;
    rLocalLHS(6,16)=0;
    rLocalLHS(6,17)=0;
    rLocalLHS(6,18)=clhs6;
    rLocalLHS(6,19)=0;
    rLocalLHS(6,20)=0;
    rLocalLHS(6,21)=clhs7;
    rLocalLHS(6,22)=0;
    rLocalLHS(6,23)=0;
    rLocalLHS(6,24)=clhs8;
    rLocalLHS(6,25)=0;
    rLocalLHS(6,26)=0;
    rLocalLHS(7,0)=0;
    rLocalLHS(7,1)=0;
    rLocalLHS(7,2)=0;
    rLocalLHS(7,3)=0;
    rLocalLHS(7,4)=0;
    rLocalLHS(7,5)=0;
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=0;
    rLocalLHS(7,10)=0;
    rLocalLHS(7,11)=0;
    rLocalLHS(7,12)=0;
    rLocalLHS(7,13)=0;
    rLocalLHS(7,14)=0;
    rLocalLHS(7,15)=0;
    rLocalLHS(7,16)=0;
    rLocalLHS(7,17)=0;
    rLocalLHS(7,18)=0;
    rLocalLHS(7,19)=clhs6;
    rLocalLHS(7,20)=0;
    rLocalLHS(7,21)=0;
    rLocalLHS(7,22)=clhs7;
    rLocalLHS(7,23)=0;
    rLocalLHS(7,24)=0;
    rLocalLHS(7,25)=clhs8;
    rLocalLHS(7,26)=0;
    rLocalLHS(8,0)=0;
    rLocalLHS(8,1)=0;
    rLocalLHS(8,2)=0;
    rLocalLHS(8,3)=0;
    rLocalLHS(8,4)=0;
    rLocalLHS(8,5)=0;
    rLocalLHS(8,6)=0;
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(8,11)=0;
    rLocalLHS(8,12)=0;
    rLocalLHS(8,13)=0;
    rLocalLHS(8,14)=0;
    rLocalLHS(8,15)=0;
    rLocalLHS(8,16)=0;
    rLocalLHS(8,17)=0;
    rLocalLHS(8,18)=0;
    rLocalLHS(8,19)=0;
    rLocalLHS(8,20)=clhs6;
    rLocalLHS(8,21)=0;
    rLocalLHS(8,22)=0;
    rLocalLHS(8,23)=clhs7;
    rLocalLHS(8,24)=0;
    rLocalLHS(8,25)=0;
    rLocalLHS(8,26)=clhs8;
    rLocalLHS(9,0)=0;
    rLocalLHS(9,1)=0;
    rLocalLHS(9,2)=0;
    rLocalLHS(9,3)=0;
    rLocalLHS(9,4)=0;
    rLocalLHS(9,5)=0;
    rLocalLHS(9,6)=0;
    rLocalLHS(9,7)=0;
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(9,11)=0;
    rLocalLHS(9,12)=0;
    rLocalLHS(9,13)=0;
    rLocalLHS(9,14)=0;
    rLocalLHS(9,15)=0;
    rLocalLHS(9,16)=0;
    rLocalLHS(9,17)=0;
    rLocalLHS(9,18)=DOperator(0,0);
    rLocalLHS(9,19)=0;
    rLocalLHS(9,20)=0;
    rLocalLHS(9,21)=DOperator(1,0);
    rLocalLHS(9,22)=0;
    rLocalLHS(9,23)=0;
    rLocalLHS(9,24)=DOperator(2,0);
    rLocalLHS(9,25)=0;
    rLocalLHS(9,26)=0;
    rLocalLHS(10,0)=0;
    rLocalLHS(10,1)=0;
    rLocalLHS(10,2)=0;
    rLocalLHS(10,3)=0;
    rLocalLHS(10,4)=0;
    rLocalLHS(10,5)=0;
    rLocalLHS(10,6)=0;
    rLocalLHS(10,7)=0;
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;
    rLocalLHS(10,11)=0;
    rLocalLHS(10,12)=0;
    rLocalLHS(10,13)=0;
    rLocalLHS(10,14)=0;
    rLocalLHS(10,15)=0;
    rLocalLHS(10,16)=0;
    rLocalLHS(10,17)=0;
    rLocalLHS(10,18)=0;
    rLocalLHS(10,19)=DOperator(0,0);
    rLocalLHS(10,20)=0;
    rLocalLHS(10,21)=0;
    rLocalLHS(10,22)=DOperator(1,0);
    rLocalLHS(10,23)=0;
    rLocalLHS(10,24)=0;
    rLocalLHS(10,25)=DOperator(2,0);
    rLocalLHS(10,26)=0;
    rLocalLHS(11,0)=0;
    rLocalLHS(11,1)=0;
    rLocalLHS(11,2)=0;
    rLocalLHS(11,3)=0;
    rLocalLHS(11,4)=0;
    rLocalLHS(11,5)=0;
    rLocalLHS(11,6)=0;
    rLocalLHS(11,7)=0;
    rLocalLHS(11,8)=0;
    rLocalLHS(11,9)=0;
    rLocalLHS(11,10)=0;
    rLocalLHS(11,11)=0;
    rLocalLHS(11,12)=0;
    rLocalLHS(11,13)=0;
    rLocalLHS(11,14)=0;
    rLocalLHS(11,15)=0;
    rLocalLHS(11,16)=0;
    rLocalLHS(11,17)=0;
    rLocalLHS(11,18)=0;
    rLocalLHS(11,19)=0;
    rLocalLHS(11,20)=DOperator(0,0);
    rLocalLHS(11,21)=0;
    rLocalLHS(11,22)=0;
    rLocalLHS(11,23)=DOperator(1,0);
    rLocalLHS(11,24)=0;
    rLocalLHS(11,25)=0;
    rLocalLHS(11,26)=DOperator(2,0);
    rLocalLHS(12,0)=0;
    rLocalLHS(12,1)=0;
    rLocalLHS(12,2)=0;
    rLocalLHS(12,3)=0;
    rLocalLHS(12,4)=0;
    rLocalLHS(12,5)=0;
    rLocalLHS(12,6)=0;
    rLocalLHS(12,7)=0;
    rLocalLHS(12,8)=0;
    rLocalLHS(12,9)=0;
    rLocalLHS(12,10)=0;
    rLocalLHS(12,11)=0;
    rLocalLHS(12,12)=0;
    rLocalLHS(12,13)=0;
    rLocalLHS(12,14)=0;
    rLocalLHS(12,15)=0;
    rLocalLHS(12,16)=0;
    rLocalLHS(12,17)=0;
    rLocalLHS(12,18)=DOperator(0,1);
    rLocalLHS(12,19)=0;
    rLocalLHS(12,20)=0;
    rLocalLHS(12,21)=DOperator(1,1);
    rLocalLHS(12,22)=0;
    rLocalLHS(12,23)=0;
    rLocalLHS(12,24)=DOperator(2,1);
    rLocalLHS(12,25)=0;
    rLocalLHS(12,26)=0;
    rLocalLHS(13,0)=0;
    rLocalLHS(13,1)=0;
    rLocalLHS(13,2)=0;
    rLocalLHS(13,3)=0;
    rLocalLHS(13,4)=0;
    rLocalLHS(13,5)=0;
    rLocalLHS(13,6)=0;
    rLocalLHS(13,7)=0;
    rLocalLHS(13,8)=0;
    rLocalLHS(13,9)=0;
    rLocalLHS(13,10)=0;
    rLocalLHS(13,11)=0;
    rLocalLHS(13,12)=0;
    rLocalLHS(13,13)=0;
    rLocalLHS(13,14)=0;
    rLocalLHS(13,15)=0;
    rLocalLHS(13,16)=0;
    rLocalLHS(13,17)=0;
    rLocalLHS(13,18)=0;
    rLocalLHS(13,19)=DOperator(0,1);
    rLocalLHS(13,20)=0;
    rLocalLHS(13,21)=0;
    rLocalLHS(13,22)=DOperator(1,1);
    rLocalLHS(13,23)=0;
    rLocalLHS(13,24)=0;
    rLocalLHS(13,25)=DOperator(2,1);
    rLocalLHS(13,26)=0;
    rLocalLHS(14,0)=0;
    rLocalLHS(14,1)=0;
    rLocalLHS(14,2)=0;
    rLocalLHS(14,3)=0;
    rLocalLHS(14,4)=0;
    rLocalLHS(14,5)=0;
    rLocalLHS(14,6)=0;
    rLocalLHS(14,7)=0;
    rLocalLHS(14,8)=0;
    rLocalLHS(14,9)=0;
    rLocalLHS(14,10)=0;
    rLocalLHS(14,11)=0;
    rLocalLHS(14,12)=0;
    rLocalLHS(14,13)=0;
    rLocalLHS(14,14)=0;
    rLocalLHS(14,15)=0;
    rLocalLHS(14,16)=0;
    rLocalLHS(14,17)=0;
    rLocalLHS(14,18)=0;
    rLocalLHS(14,19)=0;
    rLocalLHS(14,20)=DOperator(0,1);
    rLocalLHS(14,21)=0;
    rLocalLHS(14,22)=0;
    rLocalLHS(14,23)=DOperator(1,1);
    rLocalLHS(14,24)=0;
    rLocalLHS(14,25)=0;
    rLocalLHS(14,26)=DOperator(2,1);
    rLocalLHS(15,0)=0;
    rLocalLHS(15,1)=0;
    rLocalLHS(15,2)=0;
    rLocalLHS(15,3)=0;
    rLocalLHS(15,4)=0;
    rLocalLHS(15,5)=0;
    rLocalLHS(15,6)=0;
    rLocalLHS(15,7)=0;
    rLocalLHS(15,8)=0;
    rLocalLHS(15,9)=0;
    rLocalLHS(15,10)=0;
    rLocalLHS(15,11)=0;
    rLocalLHS(15,12)=0;
    rLocalLHS(15,13)=0;
    rLocalLHS(15,14)=0;
    rLocalLHS(15,15)=0;
    rLocalLHS(15,16)=0;
    rLocalLHS(15,17)=0;
    rLocalLHS(15,18)=DOperator(0,2);
    rLocalLHS(15,19)=0;
    rLocalLHS(15,20)=0;
    rLocalLHS(15,21)=DOperator(1,2);
    rLocalLHS(15,22)=0;
    rLocalLHS(15,23)=0;
    rLocalLHS(15,24)=DOperator(2,2);
    rLocalLHS(15,25)=0;
    rLocalLHS(15,26)=0;
    rLocalLHS(16,0)=0;
    rLocalLHS(16,1)=0;
    rLocalLHS(16,2)=0;
    rLocalLHS(16,3)=0;
    rLocalLHS(16,4)=0;
    rLocalLHS(16,5)=0;
    rLocalLHS(16,6)=0;
    rLocalLHS(16,7)=0;
    rLocalLHS(16,8)=0;
    rLocalLHS(16,9)=0;
    rLocalLHS(16,10)=0;
    rLocalLHS(16,11)=0;
    rLocalLHS(16,12)=0;
    rLocalLHS(16,13)=0;
    rLocalLHS(16,14)=0;
    rLocalLHS(16,15)=0;
    rLocalLHS(16,16)=0;
    rLocalLHS(16,17)=0;
    rLocalLHS(16,18)=0;
    rLocalLHS(16,19)=DOperator(0,2);
    rLocalLHS(16,20)=0;
    rLocalLHS(16,21)=0;
    rLocalLHS(16,22)=DOperator(1,2);
    rLocalLHS(16,23)=0;
    rLocalLHS(16,24)=0;
    rLocalLHS(16,25)=DOperator(2,2);
    rLocalLHS(16,26)=0;
    rLocalLHS(17,0)=0;
    rLocalLHS(17,1)=0;
    rLocalLHS(17,2)=0;
    rLocalLHS(17,3)=0;
    rLocalLHS(17,4)=0;
    rLocalLHS(17,5)=0;
    rLocalLHS(17,6)=0;
    rLocalLHS(17,7)=0;
    rLocalLHS(17,8)=0;
    rLocalLHS(17,9)=0;
    rLocalLHS(17,10)=0;
    rLocalLHS(17,11)=0;
    rLocalLHS(17,12)=0;
    rLocalLHS(17,13)=0;
    rLocalLHS(17,14)=0;
    rLocalLHS(17,15)=0;
    rLocalLHS(17,16)=0;
    rLocalLHS(17,17)=0;
    rLocalLHS(17,18)=0;
    rLocalLHS(17,19)=0;
    rLocalLHS(17,20)=DOperator(0,2);
    rLocalLHS(17,21)=0;
    rLocalLHS(17,22)=0;
    rLocalLHS(17,23)=DOperator(1,2);
    rLocalLHS(17,24)=0;
    rLocalLHS(17,25)=0;
    rLocalLHS(17,26)=DOperator(2,2);
    rLocalLHS(18,0)=clhs0;
    rLocalLHS(18,1)=0;
    rLocalLHS(18,2)=0;
    rLocalLHS(18,3)=clhs3;
    rLocalLHS(18,4)=0;
    rLocalLHS(18,5)=0;
    rLocalLHS(18,6)=clhs6;
    rLocalLHS(18,7)=0;
    rLocalLHS(18,8)=0;
    rLocalLHS(18,9)=DOperator(0,0);
    rLocalLHS(18,10)=0;
    rLocalLHS(18,11)=0;
    rLocalLHS(18,12)=DOperator(0,1);
    rLocalLHS(18,13)=0;
    rLocalLHS(18,14)=0;
    rLocalLHS(18,15)=DOperator(0,2);
    rLocalLHS(18,16)=0;
    rLocalLHS(18,17)=0;
    rLocalLHS(18,18)=0;
    rLocalLHS(18,19)=0;
    rLocalLHS(18,20)=0;
    rLocalLHS(18,21)=0;
    rLocalLHS(18,22)=0;
    rLocalLHS(18,23)=0;
    rLocalLHS(18,24)=0;
    rLocalLHS(18,25)=0;
    rLocalLHS(18,26)=0;
    rLocalLHS(19,0)=0;
    rLocalLHS(19,1)=clhs0;
    rLocalLHS(19,2)=0;
    rLocalLHS(19,3)=0;
    rLocalLHS(19,4)=clhs3;
    rLocalLHS(19,5)=0;
    rLocalLHS(19,6)=0;
    rLocalLHS(19,7)=clhs6;
    rLocalLHS(19,8)=0;
    rLocalLHS(19,9)=0;
    rLocalLHS(19,10)=DOperator(0,0);
    rLocalLHS(19,11)=0;
    rLocalLHS(19,12)=0;
    rLocalLHS(19,13)=DOperator(0,1);
    rLocalLHS(19,14)=0;
    rLocalLHS(19,15)=0;
    rLocalLHS(19,16)=DOperator(0,2);
    rLocalLHS(19,17)=0;
    rLocalLHS(19,18)=0;
    rLocalLHS(19,19)=0;
    rLocalLHS(19,20)=0;
    rLocalLHS(19,21)=0;
    rLocalLHS(19,22)=0;
    rLocalLHS(19,23)=0;
    rLocalLHS(19,24)=0;
    rLocalLHS(19,25)=0;
    rLocalLHS(19,26)=0;
    rLocalLHS(20,0)=0;
    rLocalLHS(20,1)=0;
    rLocalLHS(20,2)=clhs0;
    rLocalLHS(20,3)=0;
    rLocalLHS(20,4)=0;
    rLocalLHS(20,5)=clhs3;
    rLocalLHS(20,6)=0;
    rLocalLHS(20,7)=0;
    rLocalLHS(20,8)=clhs6;
    rLocalLHS(20,9)=0;
    rLocalLHS(20,10)=0;
    rLocalLHS(20,11)=DOperator(0,0);
    rLocalLHS(20,12)=0;
    rLocalLHS(20,13)=0;
    rLocalLHS(20,14)=DOperator(0,1);
    rLocalLHS(20,15)=0;
    rLocalLHS(20,16)=0;
    rLocalLHS(20,17)=DOperator(0,2);
    rLocalLHS(20,18)=0;
    rLocalLHS(20,19)=0;
    rLocalLHS(20,20)=0;
    rLocalLHS(20,21)=0;
    rLocalLHS(20,22)=0;
    rLocalLHS(20,23)=0;
    rLocalLHS(20,24)=0;
    rLocalLHS(20,25)=0;
    rLocalLHS(20,26)=0;
    rLocalLHS(21,0)=clhs1;
    rLocalLHS(21,1)=0;
    rLocalLHS(21,2)=0;
    rLocalLHS(21,3)=clhs4;
    rLocalLHS(21,4)=0;
    rLocalLHS(21,5)=0;
    rLocalLHS(21,6)=clhs7;
    rLocalLHS(21,7)=0;
    rLocalLHS(21,8)=0;
    rLocalLHS(21,9)=DOperator(1,0);
    rLocalLHS(21,10)=0;
    rLocalLHS(21,11)=0;
    rLocalLHS(21,12)=DOperator(1,1);
    rLocalLHS(21,13)=0;
    rLocalLHS(21,14)=0;
    rLocalLHS(21,15)=DOperator(1,2);
    rLocalLHS(21,16)=0;
    rLocalLHS(21,17)=0;
    rLocalLHS(21,18)=0;
    rLocalLHS(21,19)=0;
    rLocalLHS(21,20)=0;
    rLocalLHS(21,21)=0;
    rLocalLHS(21,22)=0;
    rLocalLHS(21,23)=0;
    rLocalLHS(21,24)=0;
    rLocalLHS(21,25)=0;
    rLocalLHS(21,26)=0;
    rLocalLHS(22,0)=0;
    rLocalLHS(22,1)=clhs1;
    rLocalLHS(22,2)=0;
    rLocalLHS(22,3)=0;
    rLocalLHS(22,4)=clhs4;
    rLocalLHS(22,5)=0;
    rLocalLHS(22,6)=0;
    rLocalLHS(22,7)=clhs7;
    rLocalLHS(22,8)=0;
    rLocalLHS(22,9)=0;
    rLocalLHS(22,10)=DOperator(1,0);
    rLocalLHS(22,11)=0;
    rLocalLHS(22,12)=0;
    rLocalLHS(22,13)=DOperator(1,1);
    rLocalLHS(22,14)=0;
    rLocalLHS(22,15)=0;
    rLocalLHS(22,16)=DOperator(1,2);
    rLocalLHS(22,17)=0;
    rLocalLHS(22,18)=0;
    rLocalLHS(22,19)=0;
    rLocalLHS(22,20)=0;
    rLocalLHS(22,21)=0;
    rLocalLHS(22,22)=0;
    rLocalLHS(22,23)=0;
    rLocalLHS(22,24)=0;
    rLocalLHS(22,25)=0;
    rLocalLHS(22,26)=0;
    rLocalLHS(23,0)=0;
    rLocalLHS(23,1)=0;
    rLocalLHS(23,2)=clhs1;
    rLocalLHS(23,3)=0;
    rLocalLHS(23,4)=0;
    rLocalLHS(23,5)=clhs4;
    rLocalLHS(23,6)=0;
    rLocalLHS(23,7)=0;
    rLocalLHS(23,8)=clhs7;
    rLocalLHS(23,9)=0;
    rLocalLHS(23,10)=0;
    rLocalLHS(23,11)=DOperator(1,0);
    rLocalLHS(23,12)=0;
    rLocalLHS(23,13)=0;
    rLocalLHS(23,14)=DOperator(1,1);
    rLocalLHS(23,15)=0;
    rLocalLHS(23,16)=0;
    rLocalLHS(23,17)=DOperator(1,2);
    rLocalLHS(23,18)=0;
    rLocalLHS(23,19)=0;
    rLocalLHS(23,20)=0;
    rLocalLHS(23,21)=0;
    rLocalLHS(23,22)=0;
    rLocalLHS(23,23)=0;
    rLocalLHS(23,24)=0;
    rLocalLHS(23,25)=0;
    rLocalLHS(23,26)=0;
    rLocalLHS(24,0)=clhs2;
    rLocalLHS(24,1)=0;
    rLocalLHS(24,2)=0;
    rLocalLHS(24,3)=clhs5;
    rLocalLHS(24,4)=0;
    rLocalLHS(24,5)=0;
    rLocalLHS(24,6)=clhs8;
    rLocalLHS(24,7)=0;
    rLocalLHS(24,8)=0;
    rLocalLHS(24,9)=DOperator(2,0);
    rLocalLHS(24,10)=0;
    rLocalLHS(24,11)=0;
    rLocalLHS(24,12)=DOperator(2,1);
    rLocalLHS(24,13)=0;
    rLocalLHS(24,14)=0;
    rLocalLHS(24,15)=DOperator(2,2);
    rLocalLHS(24,16)=0;
    rLocalLHS(24,17)=0;
    rLocalLHS(24,18)=0;
    rLocalLHS(24,19)=0;
    rLocalLHS(24,20)=0;
    rLocalLHS(24,21)=0;
    rLocalLHS(24,22)=0;
    rLocalLHS(24,23)=0;
    rLocalLHS(24,24)=0;
    rLocalLHS(24,25)=0;
    rLocalLHS(24,26)=0;
    rLocalLHS(25,0)=0;
    rLocalLHS(25,1)=clhs2;
    rLocalLHS(25,2)=0;
    rLocalLHS(25,3)=0;
    rLocalLHS(25,4)=clhs5;
    rLocalLHS(25,5)=0;
    rLocalLHS(25,6)=0;
    rLocalLHS(25,7)=clhs8;
    rLocalLHS(25,8)=0;
    rLocalLHS(25,9)=0;
    rLocalLHS(25,10)=DOperator(2,0);
    rLocalLHS(25,11)=0;
    rLocalLHS(25,12)=0;
    rLocalLHS(25,13)=DOperator(2,1);
    rLocalLHS(25,14)=0;
    rLocalLHS(25,15)=0;
    rLocalLHS(25,16)=DOperator(2,2);
    rLocalLHS(25,17)=0;
    rLocalLHS(25,18)=0;
    rLocalLHS(25,19)=0;
    rLocalLHS(25,20)=0;
    rLocalLHS(25,21)=0;
    rLocalLHS(25,22)=0;
    rLocalLHS(25,23)=0;
    rLocalLHS(25,24)=0;
    rLocalLHS(25,25)=0;
    rLocalLHS(25,26)=0;
    rLocalLHS(26,0)=0;
    rLocalLHS(26,1)=0;
    rLocalLHS(26,2)=clhs2;
    rLocalLHS(26,3)=0;
    rLocalLHS(26,4)=0;
    rLocalLHS(26,5)=clhs5;
    rLocalLHS(26,6)=0;
    rLocalLHS(26,7)=0;
    rLocalLHS(26,8)=clhs8;
    rLocalLHS(26,9)=0;
    rLocalLHS(26,10)=0;
    rLocalLHS(26,11)=DOperator(2,0);
    rLocalLHS(26,12)=0;
    rLocalLHS(26,13)=0;
    rLocalLHS(26,14)=DOperator(2,1);
    rLocalLHS(26,15)=0;
    rLocalLHS(26,16)=0;
    rLocalLHS(26,17)=DOperator(2,2);
    rLocalLHS(26,18)=0;
    rLocalLHS(26,19)=0;
    rLocalLHS(26,20)=0;
    rLocalLHS(26,21)=0;
    rLocalLHS(26,22)=0;
    rLocalLHS(26,23)=0;
    rLocalLHS(26,24)=0;
    rLocalLHS(26,25)=0;
    rLocalLHS(26,26)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,8,8>::CalculateLocalLHS<MeshTyingMortarCondition<3,8,8>::ScalarValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(3,0);
    const double clhs4 =     -MOperator(0,1);
    const double clhs5 =     -MOperator(1,1);
    const double clhs6 =     -MOperator(2,1);
    const double clhs7 =     -MOperator(3,1);
    const double clhs8 =     -MOperator(0,2);
    const double clhs9 =     -MOperator(1,2);
    const double clhs10 =     -MOperator(2,2);
    const double clhs11 =     -MOperator(3,2);
    const double clhs12 =     -MOperator(0,3);
    const double clhs13 =     -MOperator(1,3);
    const double clhs14 =     -MOperator(2,3);
    const double clhs15 =     -MOperator(3,3);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=0;
    rLocalLHS(0,8)=clhs0;
    rLocalLHS(0,9)=clhs1;
    rLocalLHS(0,10)=clhs2;
    rLocalLHS(0,11)=clhs3;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=0;
    rLocalLHS(1,8)=clhs4;
    rLocalLHS(1,9)=clhs5;
    rLocalLHS(1,10)=clhs6;
    rLocalLHS(1,11)=clhs7;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=0;
    rLocalLHS(2,8)=clhs8;
    rLocalLHS(2,9)=clhs9;
    rLocalLHS(2,10)=clhs10;
    rLocalLHS(2,11)=clhs11;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=0;
    rLocalLHS(3,8)=clhs12;
    rLocalLHS(3,9)=clhs13;
    rLocalLHS(3,10)=clhs14;
    rLocalLHS(3,11)=clhs15;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=0;
    rLocalLHS(4,8)=DOperator(0,0);
    rLocalLHS(4,9)=DOperator(1,0);
    rLocalLHS(4,10)=DOperator(2,0);
    rLocalLHS(4,11)=DOperator(3,0);
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=0;
    rLocalLHS(5,8)=DOperator(0,1);
    rLocalLHS(5,9)=DOperator(1,1);
    rLocalLHS(5,10)=DOperator(2,1);
    rLocalLHS(5,11)=DOperator(3,1);
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=DOperator(0,2);
    rLocalLHS(6,9)=DOperator(1,2);
    rLocalLHS(6,10)=DOperator(2,2);
    rLocalLHS(6,11)=DOperator(3,2);
    rLocalLHS(7,0)=0;
    rLocalLHS(7,1)=0;
    rLocalLHS(7,2)=0;
    rLocalLHS(7,3)=0;
    rLocalLHS(7,4)=0;
    rLocalLHS(7,5)=0;
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=DOperator(0,3);
    rLocalLHS(7,9)=DOperator(1,3);
    rLocalLHS(7,10)=DOperator(2,3);
    rLocalLHS(7,11)=DOperator(3,3);
    rLocalLHS(8,0)=clhs0;
    rLocalLHS(8,1)=clhs4;
    rLocalLHS(8,2)=clhs8;
    rLocalLHS(8,3)=clhs12;
    rLocalLHS(8,4)=DOperator(0,0);
    rLocalLHS(8,5)=DOperator(0,1);
    rLocalLHS(8,6)=DOperator(0,2);
    rLocalLHS(8,7)=DOperator(0,3);
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(8,11)=0;
    rLocalLHS(9,0)=clhs1;
    rLocalLHS(9,1)=clhs5;
    rLocalLHS(9,2)=clhs9;
    rLocalLHS(9,3)=clhs13;
    rLocalLHS(9,4)=DOperator(1,0);
    rLocalLHS(9,5)=DOperator(1,1);
    rLocalLHS(9,6)=DOperator(1,2);
    rLocalLHS(9,7)=DOperator(1,3);
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(9,11)=0;
    rLocalLHS(10,0)=clhs2;
    rLocalLHS(10,1)=clhs6;
    rLocalLHS(10,2)=clhs10;
    rLocalLHS(10,3)=clhs14;
    rLocalLHS(10,4)=DOperator(2,0);
    rLocalLHS(10,5)=DOperator(2,1);
    rLocalLHS(10,6)=DOperator(2,2);
    rLocalLHS(10,7)=DOperator(2,3);
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;
    rLocalLHS(10,11)=0;
    rLocalLHS(11,0)=clhs3;
    rLocalLHS(11,1)=clhs7;
    rLocalLHS(11,2)=clhs11;
    rLocalLHS(11,3)=clhs15;
    rLocalLHS(11,4)=DOperator(3,0);
    rLocalLHS(11,5)=DOperator(3,1);
    rLocalLHS(11,6)=DOperator(3,2);
    rLocalLHS(11,7)=DOperator(3,3);
    rLocalLHS(11,8)=0;
    rLocalLHS(11,9)=0;
    rLocalLHS(11,10)=0;
    rLocalLHS(11,11)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,8,8>::CalculateLocalLHS<MeshTyingMortarCondition<3,8,8>::Vector3DValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(3,0);
    const double clhs4 =     -MOperator(0,1);
    const double clhs5 =     -MOperator(1,1);
    const double clhs6 =     -MOperator(2,1);
    const double clhs7 =     -MOperator(3,1);
    const double clhs8 =     -MOperator(0,2);
    const double clhs9 =     -MOperator(1,2);
    const double clhs10 =     -MOperator(2,2);
    const double clhs11 =     -MOperator(3,2);
    const double clhs12 =     -MOperator(0,3);
    const double clhs13 =     -MOperator(1,3);
    const double clhs14 =     -MOperator(2,3);
    const double clhs15 =     -MOperator(3,3);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=0;
    rLocalLHS(0,8)=0;
    rLocalLHS(0,9)=0;
    rLocalLHS(0,10)=0;
    rLocalLHS(0,11)=0;
    rLocalLHS(0,12)=0;
    rLocalLHS(0,13)=0;
    rLocalLHS(0,14)=0;
    rLocalLHS(0,15)=0;
    rLocalLHS(0,16)=0;
    rLocalLHS(0,17)=0;
    rLocalLHS(0,18)=0;
    rLocalLHS(0,19)=0;
    rLocalLHS(0,20)=0;
    rLocalLHS(0,21)=0;
    rLocalLHS(0,22)=0;
    rLocalLHS(0,23)=0;
    rLocalLHS(0,24)=clhs0;
    rLocalLHS(0,25)=0;
    rLocalLHS(0,26)=0;
    rLocalLHS(0,27)=clhs1;
    rLocalLHS(0,28)=0;
    rLocalLHS(0,29)=0;
    rLocalLHS(0,30)=clhs2;
    rLocalLHS(0,31)=0;
    rLocalLHS(0,32)=0;
    rLocalLHS(0,33)=clhs3;
    rLocalLHS(0,34)=0;
    rLocalLHS(0,35)=0;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=0;
    rLocalLHS(1,8)=0;
    rLocalLHS(1,9)=0;
    rLocalLHS(1,10)=0;
    rLocalLHS(1,11)=0;
    rLocalLHS(1,12)=0;
    rLocalLHS(1,13)=0;
    rLocalLHS(1,14)=0;
    rLocalLHS(1,15)=0;
    rLocalLHS(1,16)=0;
    rLocalLHS(1,17)=0;
    rLocalLHS(1,18)=0;
    rLocalLHS(1,19)=0;
    rLocalLHS(1,20)=0;
    rLocalLHS(1,21)=0;
    rLocalLHS(1,22)=0;
    rLocalLHS(1,23)=0;
    rLocalLHS(1,24)=0;
    rLocalLHS(1,25)=clhs0;
    rLocalLHS(1,26)=0;
    rLocalLHS(1,27)=0;
    rLocalLHS(1,28)=clhs1;
    rLocalLHS(1,29)=0;
    rLocalLHS(1,30)=0;
    rLocalLHS(1,31)=clhs2;
    rLocalLHS(1,32)=0;
    rLocalLHS(1,33)=0;
    rLocalLHS(1,34)=clhs3;
    rLocalLHS(1,35)=0;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=0;
    rLocalLHS(2,8)=0;
    rLocalLHS(2,9)=0;
    rLocalLHS(2,10)=0;
    rLocalLHS(2,11)=0;
    rLocalLHS(2,12)=0;
    rLocalLHS(2,13)=0;
    rLocalLHS(2,14)=0;
    rLocalLHS(2,15)=0;
    rLocalLHS(2,16)=0;
    rLocalLHS(2,17)=0;
    rLocalLHS(2,18)=0;
    rLocalLHS(2,19)=0;
    rLocalLHS(2,20)=0;
    rLocalLHS(2,21)=0;
    rLocalLHS(2,22)=0;
    rLocalLHS(2,23)=0;
    rLocalLHS(2,24)=0;
    rLocalLHS(2,25)=0;
    rLocalLHS(2,26)=clhs0;
    rLocalLHS(2,27)=0;
    rLocalLHS(2,28)=0;
    rLocalLHS(2,29)=clhs1;
    rLocalLHS(2,30)=0;
    rLocalLHS(2,31)=0;
    rLocalLHS(2,32)=clhs2;
    rLocalLHS(2,33)=0;
    rLocalLHS(2,34)=0;
    rLocalLHS(2,35)=clhs3;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=0;
    rLocalLHS(3,8)=0;
    rLocalLHS(3,9)=0;
    rLocalLHS(3,10)=0;
    rLocalLHS(3,11)=0;
    rLocalLHS(3,12)=0;
    rLocalLHS(3,13)=0;
    rLocalLHS(3,14)=0;
    rLocalLHS(3,15)=0;
    rLocalLHS(3,16)=0;
    rLocalLHS(3,17)=0;
    rLocalLHS(3,18)=0;
    rLocalLHS(3,19)=0;
    rLocalLHS(3,20)=0;
    rLocalLHS(3,21)=0;
    rLocalLHS(3,22)=0;
    rLocalLHS(3,23)=0;
    rLocalLHS(3,24)=clhs4;
    rLocalLHS(3,25)=0;
    rLocalLHS(3,26)=0;
    rLocalLHS(3,27)=clhs5;
    rLocalLHS(3,28)=0;
    rLocalLHS(3,29)=0;
    rLocalLHS(3,30)=clhs6;
    rLocalLHS(3,31)=0;
    rLocalLHS(3,32)=0;
    rLocalLHS(3,33)=clhs7;
    rLocalLHS(3,34)=0;
    rLocalLHS(3,35)=0;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=0;
    rLocalLHS(4,8)=0;
    rLocalLHS(4,9)=0;
    rLocalLHS(4,10)=0;
    rLocalLHS(4,11)=0;
    rLocalLHS(4,12)=0;
    rLocalLHS(4,13)=0;
    rLocalLHS(4,14)=0;
    rLocalLHS(4,15)=0;
    rLocalLHS(4,16)=0;
    rLocalLHS(4,17)=0;
    rLocalLHS(4,18)=0;
    rLocalLHS(4,19)=0;
    rLocalLHS(4,20)=0;
    rLocalLHS(4,21)=0;
    rLocalLHS(4,22)=0;
    rLocalLHS(4,23)=0;
    rLocalLHS(4,24)=0;
    rLocalLHS(4,25)=clhs4;
    rLocalLHS(4,26)=0;
    rLocalLHS(4,27)=0;
    rLocalLHS(4,28)=clhs5;
    rLocalLHS(4,29)=0;
    rLocalLHS(4,30)=0;
    rLocalLHS(4,31)=clhs6;
    rLocalLHS(4,32)=0;
    rLocalLHS(4,33)=0;
    rLocalLHS(4,34)=clhs7;
    rLocalLHS(4,35)=0;
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=0;
    rLocalLHS(5,8)=0;
    rLocalLHS(5,9)=0;
    rLocalLHS(5,10)=0;
    rLocalLHS(5,11)=0;
    rLocalLHS(5,12)=0;
    rLocalLHS(5,13)=0;
    rLocalLHS(5,14)=0;
    rLocalLHS(5,15)=0;
    rLocalLHS(5,16)=0;
    rLocalLHS(5,17)=0;
    rLocalLHS(5,18)=0;
    rLocalLHS(5,19)=0;
    rLocalLHS(5,20)=0;
    rLocalLHS(5,21)=0;
    rLocalLHS(5,22)=0;
    rLocalLHS(5,23)=0;
    rLocalLHS(5,24)=0;
    rLocalLHS(5,25)=0;
    rLocalLHS(5,26)=clhs4;
    rLocalLHS(5,27)=0;
    rLocalLHS(5,28)=0;
    rLocalLHS(5,29)=clhs5;
    rLocalLHS(5,30)=0;
    rLocalLHS(5,31)=0;
    rLocalLHS(5,32)=clhs6;
    rLocalLHS(5,33)=0;
    rLocalLHS(5,34)=0;
    rLocalLHS(5,35)=clhs7;
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=0;
    rLocalLHS(6,9)=0;
    rLocalLHS(6,10)=0;
    rLocalLHS(6,11)=0;
    rLocalLHS(6,12)=0;
    rLocalLHS(6,13)=0;
    rLocalLHS(6,14)=0;
    rLocalLHS(6,15)=0;
    rLocalLHS(6,16)=0;
    rLocalLHS(6,17)=0;
    rLocalLHS(6,18)=0;
    rLocalLHS(6,19)=0;
    rLocalLHS(6,20)=0;
    rLocalLHS(6,21)=0;
    rLocalLHS(6,22)=0;
    rLocalLHS(6,23)=0;
    rLocalLHS(6,24)=clhs8;
    rLocalLHS(6,25)=0;
    rLocalLHS(6,26)=0;
    rLocalLHS(6,27)=clhs9;
    rLocalLHS(6,28)=0;
    rLocalLHS(6,29)=0;
    rLocalLHS(6,30)=clhs10;
    rLocalLHS(6,31)=0;
    rLocalLHS(6,32)=0;
    rLocalLHS(6,33)=clhs11;
    rLocalLHS(6,34)=0;
    rLocalLHS(6,35)=0;
    rLocalLHS(7,0)=0;
    rLocalLHS(7,1)=0;
    rLocalLHS(7,2)=0;
    rLocalLHS(7,3)=0;
    rLocalLHS(7,4)=0;
    rLocalLHS(7,5)=0;
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=0;
    rLocalLHS(7,10)=0;
    rLocalLHS(7,11)=0;
    rLocalLHS(7,12)=0;
    rLocalLHS(7,13)=0;
    rLocalLHS(7,14)=0;
    rLocalLHS(7,15)=0;
    rLocalLHS(7,16)=0;
    rLocalLHS(7,17)=0;
    rLocalLHS(7,18)=0;
    rLocalLHS(7,19)=0;
    rLocalLHS(7,20)=0;
    rLocalLHS(7,21)=0;
    rLocalLHS(7,22)=0;
    rLocalLHS(7,23)=0;
    rLocalLHS(7,24)=0;
    rLocalLHS(7,25)=clhs8;
    rLocalLHS(7,26)=0;
    rLocalLHS(7,27)=0;
    rLocalLHS(7,28)=clhs9;
    rLocalLHS(7,29)=0;
    rLocalLHS(7,30)=0;
    rLocalLHS(7,31)=clhs10;
    rLocalLHS(7,32)=0;
    rLocalLHS(7,33)=0;
    rLocalLHS(7,34)=clhs11;
    rLocalLHS(7,35)=0;
    rLocalLHS(8,0)=0;
    rLocalLHS(8,1)=0;
    rLocalLHS(8,2)=0;
    rLocalLHS(8,3)=0;
    rLocalLHS(8,4)=0;
    rLocalLHS(8,5)=0;
    rLocalLHS(8,6)=0;
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(8,11)=0;
    rLocalLHS(8,12)=0;
    rLocalLHS(8,13)=0;
    rLocalLHS(8,14)=0;
    rLocalLHS(8,15)=0;
    rLocalLHS(8,16)=0;
    rLocalLHS(8,17)=0;
    rLocalLHS(8,18)=0;
    rLocalLHS(8,19)=0;
    rLocalLHS(8,20)=0;
    rLocalLHS(8,21)=0;
    rLocalLHS(8,22)=0;
    rLocalLHS(8,23)=0;
    rLocalLHS(8,24)=0;
    rLocalLHS(8,25)=0;
    rLocalLHS(8,26)=clhs8;
    rLocalLHS(8,27)=0;
    rLocalLHS(8,28)=0;
    rLocalLHS(8,29)=clhs9;
    rLocalLHS(8,30)=0;
    rLocalLHS(8,31)=0;
    rLocalLHS(8,32)=clhs10;
    rLocalLHS(8,33)=0;
    rLocalLHS(8,34)=0;
    rLocalLHS(8,35)=clhs11;
    rLocalLHS(9,0)=0;
    rLocalLHS(9,1)=0;
    rLocalLHS(9,2)=0;
    rLocalLHS(9,3)=0;
    rLocalLHS(9,4)=0;
    rLocalLHS(9,5)=0;
    rLocalLHS(9,6)=0;
    rLocalLHS(9,7)=0;
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(9,11)=0;
    rLocalLHS(9,12)=0;
    rLocalLHS(9,13)=0;
    rLocalLHS(9,14)=0;
    rLocalLHS(9,15)=0;
    rLocalLHS(9,16)=0;
    rLocalLHS(9,17)=0;
    rLocalLHS(9,18)=0;
    rLocalLHS(9,19)=0;
    rLocalLHS(9,20)=0;
    rLocalLHS(9,21)=0;
    rLocalLHS(9,22)=0;
    rLocalLHS(9,23)=0;
    rLocalLHS(9,24)=clhs12;
    rLocalLHS(9,25)=0;
    rLocalLHS(9,26)=0;
    rLocalLHS(9,27)=clhs13;
    rLocalLHS(9,28)=0;
    rLocalLHS(9,29)=0;
    rLocalLHS(9,30)=clhs14;
    rLocalLHS(9,31)=0;
    rLocalLHS(9,32)=0;
    rLocalLHS(9,33)=clhs15;
    rLocalLHS(9,34)=0;
    rLocalLHS(9,35)=0;
    rLocalLHS(10,0)=0;
    rLocalLHS(10,1)=0;
    rLocalLHS(10,2)=0;
    rLocalLHS(10,3)=0;
    rLocalLHS(10,4)=0;
    rLocalLHS(10,5)=0;
    rLocalLHS(10,6)=0;
    rLocalLHS(10,7)=0;
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;
    rLocalLHS(10,11)=0;
    rLocalLHS(10,12)=0;
    rLocalLHS(10,13)=0;
    rLocalLHS(10,14)=0;
    rLocalLHS(10,15)=0;
    rLocalLHS(10,16)=0;
    rLocalLHS(10,17)=0;
    rLocalLHS(10,18)=0;
    rLocalLHS(10,19)=0;
    rLocalLHS(10,20)=0;
    rLocalLHS(10,21)=0;
    rLocalLHS(10,22)=0;
    rLocalLHS(10,23)=0;
    rLocalLHS(10,24)=0;
    rLocalLHS(10,25)=clhs12;
    rLocalLHS(10,26)=0;
    rLocalLHS(10,27)=0;
    rLocalLHS(10,28)=clhs13;
    rLocalLHS(10,29)=0;
    rLocalLHS(10,30)=0;
    rLocalLHS(10,31)=clhs14;
    rLocalLHS(10,32)=0;
    rLocalLHS(10,33)=0;
    rLocalLHS(10,34)=clhs15;
    rLocalLHS(10,35)=0;
    rLocalLHS(11,0)=0;
    rLocalLHS(11,1)=0;
    rLocalLHS(11,2)=0;
    rLocalLHS(11,3)=0;
    rLocalLHS(11,4)=0;
    rLocalLHS(11,5)=0;
    rLocalLHS(11,6)=0;
    rLocalLHS(11,7)=0;
    rLocalLHS(11,8)=0;
    rLocalLHS(11,9)=0;
    rLocalLHS(11,10)=0;
    rLocalLHS(11,11)=0;
    rLocalLHS(11,12)=0;
    rLocalLHS(11,13)=0;
    rLocalLHS(11,14)=0;
    rLocalLHS(11,15)=0;
    rLocalLHS(11,16)=0;
    rLocalLHS(11,17)=0;
    rLocalLHS(11,18)=0;
    rLocalLHS(11,19)=0;
    rLocalLHS(11,20)=0;
    rLocalLHS(11,21)=0;
    rLocalLHS(11,22)=0;
    rLocalLHS(11,23)=0;
    rLocalLHS(11,24)=0;
    rLocalLHS(11,25)=0;
    rLocalLHS(11,26)=clhs12;
    rLocalLHS(11,27)=0;
    rLocalLHS(11,28)=0;
    rLocalLHS(11,29)=clhs13;
    rLocalLHS(11,30)=0;
    rLocalLHS(11,31)=0;
    rLocalLHS(11,32)=clhs14;
    rLocalLHS(11,33)=0;
    rLocalLHS(11,34)=0;
    rLocalLHS(11,35)=clhs15;
    rLocalLHS(12,0)=0;
    rLocalLHS(12,1)=0;
    rLocalLHS(12,2)=0;
    rLocalLHS(12,3)=0;
    rLocalLHS(12,4)=0;
    rLocalLHS(12,5)=0;
    rLocalLHS(12,6)=0;
    rLocalLHS(12,7)=0;
    rLocalLHS(12,8)=0;
    rLocalLHS(12,9)=0;
    rLocalLHS(12,10)=0;
    rLocalLHS(12,11)=0;
    rLocalLHS(12,12)=0;
    rLocalLHS(12,13)=0;
    rLocalLHS(12,14)=0;
    rLocalLHS(12,15)=0;
    rLocalLHS(12,16)=0;
    rLocalLHS(12,17)=0;
    rLocalLHS(12,18)=0;
    rLocalLHS(12,19)=0;
    rLocalLHS(12,20)=0;
    rLocalLHS(12,21)=0;
    rLocalLHS(12,22)=0;
    rLocalLHS(12,23)=0;
    rLocalLHS(12,24)=DOperator(0,0);
    rLocalLHS(12,25)=0;
    rLocalLHS(12,26)=0;
    rLocalLHS(12,27)=DOperator(1,0);
    rLocalLHS(12,28)=0;
    rLocalLHS(12,29)=0;
    rLocalLHS(12,30)=DOperator(2,0);
    rLocalLHS(12,31)=0;
    rLocalLHS(12,32)=0;
    rLocalLHS(12,33)=DOperator(3,0);
    rLocalLHS(12,34)=0;
    rLocalLHS(12,35)=0;
    rLocalLHS(13,0)=0;
    rLocalLHS(13,1)=0;
    rLocalLHS(13,2)=0;
    rLocalLHS(13,3)=0;
    rLocalLHS(13,4)=0;
    rLocalLHS(13,5)=0;
    rLocalLHS(13,6)=0;
    rLocalLHS(13,7)=0;
    rLocalLHS(13,8)=0;
    rLocalLHS(13,9)=0;
    rLocalLHS(13,10)=0;
    rLocalLHS(13,11)=0;
    rLocalLHS(13,12)=0;
    rLocalLHS(13,13)=0;
    rLocalLHS(13,14)=0;
    rLocalLHS(13,15)=0;
    rLocalLHS(13,16)=0;
    rLocalLHS(13,17)=0;
    rLocalLHS(13,18)=0;
    rLocalLHS(13,19)=0;
    rLocalLHS(13,20)=0;
    rLocalLHS(13,21)=0;
    rLocalLHS(13,22)=0;
    rLocalLHS(13,23)=0;
    rLocalLHS(13,24)=0;
    rLocalLHS(13,25)=DOperator(0,0);
    rLocalLHS(13,26)=0;
    rLocalLHS(13,27)=0;
    rLocalLHS(13,28)=DOperator(1,0);
    rLocalLHS(13,29)=0;
    rLocalLHS(13,30)=0;
    rLocalLHS(13,31)=DOperator(2,0);
    rLocalLHS(13,32)=0;
    rLocalLHS(13,33)=0;
    rLocalLHS(13,34)=DOperator(3,0);
    rLocalLHS(13,35)=0;
    rLocalLHS(14,0)=0;
    rLocalLHS(14,1)=0;
    rLocalLHS(14,2)=0;
    rLocalLHS(14,3)=0;
    rLocalLHS(14,4)=0;
    rLocalLHS(14,5)=0;
    rLocalLHS(14,6)=0;
    rLocalLHS(14,7)=0;
    rLocalLHS(14,8)=0;
    rLocalLHS(14,9)=0;
    rLocalLHS(14,10)=0;
    rLocalLHS(14,11)=0;
    rLocalLHS(14,12)=0;
    rLocalLHS(14,13)=0;
    rLocalLHS(14,14)=0;
    rLocalLHS(14,15)=0;
    rLocalLHS(14,16)=0;
    rLocalLHS(14,17)=0;
    rLocalLHS(14,18)=0;
    rLocalLHS(14,19)=0;
    rLocalLHS(14,20)=0;
    rLocalLHS(14,21)=0;
    rLocalLHS(14,22)=0;
    rLocalLHS(14,23)=0;
    rLocalLHS(14,24)=0;
    rLocalLHS(14,25)=0;
    rLocalLHS(14,26)=DOperator(0,0);
    rLocalLHS(14,27)=0;
    rLocalLHS(14,28)=0;
    rLocalLHS(14,29)=DOperator(1,0);
    rLocalLHS(14,30)=0;
    rLocalLHS(14,31)=0;
    rLocalLHS(14,32)=DOperator(2,0);
    rLocalLHS(14,33)=0;
    rLocalLHS(14,34)=0;
    rLocalLHS(14,35)=DOperator(3,0);
    rLocalLHS(15,0)=0;
    rLocalLHS(15,1)=0;
    rLocalLHS(15,2)=0;
    rLocalLHS(15,3)=0;
    rLocalLHS(15,4)=0;
    rLocalLHS(15,5)=0;
    rLocalLHS(15,6)=0;
    rLocalLHS(15,7)=0;
    rLocalLHS(15,8)=0;
    rLocalLHS(15,9)=0;
    rLocalLHS(15,10)=0;
    rLocalLHS(15,11)=0;
    rLocalLHS(15,12)=0;
    rLocalLHS(15,13)=0;
    rLocalLHS(15,14)=0;
    rLocalLHS(15,15)=0;
    rLocalLHS(15,16)=0;
    rLocalLHS(15,17)=0;
    rLocalLHS(15,18)=0;
    rLocalLHS(15,19)=0;
    rLocalLHS(15,20)=0;
    rLocalLHS(15,21)=0;
    rLocalLHS(15,22)=0;
    rLocalLHS(15,23)=0;
    rLocalLHS(15,24)=DOperator(0,1);
    rLocalLHS(15,25)=0;
    rLocalLHS(15,26)=0;
    rLocalLHS(15,27)=DOperator(1,1);
    rLocalLHS(15,28)=0;
    rLocalLHS(15,29)=0;
    rLocalLHS(15,30)=DOperator(2,1);
    rLocalLHS(15,31)=0;
    rLocalLHS(15,32)=0;
    rLocalLHS(15,33)=DOperator(3,1);
    rLocalLHS(15,34)=0;
    rLocalLHS(15,35)=0;
    rLocalLHS(16,0)=0;
    rLocalLHS(16,1)=0;
    rLocalLHS(16,2)=0;
    rLocalLHS(16,3)=0;
    rLocalLHS(16,4)=0;
    rLocalLHS(16,5)=0;
    rLocalLHS(16,6)=0;
    rLocalLHS(16,7)=0;
    rLocalLHS(16,8)=0;
    rLocalLHS(16,9)=0;
    rLocalLHS(16,10)=0;
    rLocalLHS(16,11)=0;
    rLocalLHS(16,12)=0;
    rLocalLHS(16,13)=0;
    rLocalLHS(16,14)=0;
    rLocalLHS(16,15)=0;
    rLocalLHS(16,16)=0;
    rLocalLHS(16,17)=0;
    rLocalLHS(16,18)=0;
    rLocalLHS(16,19)=0;
    rLocalLHS(16,20)=0;
    rLocalLHS(16,21)=0;
    rLocalLHS(16,22)=0;
    rLocalLHS(16,23)=0;
    rLocalLHS(16,24)=0;
    rLocalLHS(16,25)=DOperator(0,1);
    rLocalLHS(16,26)=0;
    rLocalLHS(16,27)=0;
    rLocalLHS(16,28)=DOperator(1,1);
    rLocalLHS(16,29)=0;
    rLocalLHS(16,30)=0;
    rLocalLHS(16,31)=DOperator(2,1);
    rLocalLHS(16,32)=0;
    rLocalLHS(16,33)=0;
    rLocalLHS(16,34)=DOperator(3,1);
    rLocalLHS(16,35)=0;
    rLocalLHS(17,0)=0;
    rLocalLHS(17,1)=0;
    rLocalLHS(17,2)=0;
    rLocalLHS(17,3)=0;
    rLocalLHS(17,4)=0;
    rLocalLHS(17,5)=0;
    rLocalLHS(17,6)=0;
    rLocalLHS(17,7)=0;
    rLocalLHS(17,8)=0;
    rLocalLHS(17,9)=0;
    rLocalLHS(17,10)=0;
    rLocalLHS(17,11)=0;
    rLocalLHS(17,12)=0;
    rLocalLHS(17,13)=0;
    rLocalLHS(17,14)=0;
    rLocalLHS(17,15)=0;
    rLocalLHS(17,16)=0;
    rLocalLHS(17,17)=0;
    rLocalLHS(17,18)=0;
    rLocalLHS(17,19)=0;
    rLocalLHS(17,20)=0;
    rLocalLHS(17,21)=0;
    rLocalLHS(17,22)=0;
    rLocalLHS(17,23)=0;
    rLocalLHS(17,24)=0;
    rLocalLHS(17,25)=0;
    rLocalLHS(17,26)=DOperator(0,1);
    rLocalLHS(17,27)=0;
    rLocalLHS(17,28)=0;
    rLocalLHS(17,29)=DOperator(1,1);
    rLocalLHS(17,30)=0;
    rLocalLHS(17,31)=0;
    rLocalLHS(17,32)=DOperator(2,1);
    rLocalLHS(17,33)=0;
    rLocalLHS(17,34)=0;
    rLocalLHS(17,35)=DOperator(3,1);
    rLocalLHS(18,0)=0;
    rLocalLHS(18,1)=0;
    rLocalLHS(18,2)=0;
    rLocalLHS(18,3)=0;
    rLocalLHS(18,4)=0;
    rLocalLHS(18,5)=0;
    rLocalLHS(18,6)=0;
    rLocalLHS(18,7)=0;
    rLocalLHS(18,8)=0;
    rLocalLHS(18,9)=0;
    rLocalLHS(18,10)=0;
    rLocalLHS(18,11)=0;
    rLocalLHS(18,12)=0;
    rLocalLHS(18,13)=0;
    rLocalLHS(18,14)=0;
    rLocalLHS(18,15)=0;
    rLocalLHS(18,16)=0;
    rLocalLHS(18,17)=0;
    rLocalLHS(18,18)=0;
    rLocalLHS(18,19)=0;
    rLocalLHS(18,20)=0;
    rLocalLHS(18,21)=0;
    rLocalLHS(18,22)=0;
    rLocalLHS(18,23)=0;
    rLocalLHS(18,24)=DOperator(0,2);
    rLocalLHS(18,25)=0;
    rLocalLHS(18,26)=0;
    rLocalLHS(18,27)=DOperator(1,2);
    rLocalLHS(18,28)=0;
    rLocalLHS(18,29)=0;
    rLocalLHS(18,30)=DOperator(2,2);
    rLocalLHS(18,31)=0;
    rLocalLHS(18,32)=0;
    rLocalLHS(18,33)=DOperator(3,2);
    rLocalLHS(18,34)=0;
    rLocalLHS(18,35)=0;
    rLocalLHS(19,0)=0;
    rLocalLHS(19,1)=0;
    rLocalLHS(19,2)=0;
    rLocalLHS(19,3)=0;
    rLocalLHS(19,4)=0;
    rLocalLHS(19,5)=0;
    rLocalLHS(19,6)=0;
    rLocalLHS(19,7)=0;
    rLocalLHS(19,8)=0;
    rLocalLHS(19,9)=0;
    rLocalLHS(19,10)=0;
    rLocalLHS(19,11)=0;
    rLocalLHS(19,12)=0;
    rLocalLHS(19,13)=0;
    rLocalLHS(19,14)=0;
    rLocalLHS(19,15)=0;
    rLocalLHS(19,16)=0;
    rLocalLHS(19,17)=0;
    rLocalLHS(19,18)=0;
    rLocalLHS(19,19)=0;
    rLocalLHS(19,20)=0;
    rLocalLHS(19,21)=0;
    rLocalLHS(19,22)=0;
    rLocalLHS(19,23)=0;
    rLocalLHS(19,24)=0;
    rLocalLHS(19,25)=DOperator(0,2);
    rLocalLHS(19,26)=0;
    rLocalLHS(19,27)=0;
    rLocalLHS(19,28)=DOperator(1,2);
    rLocalLHS(19,29)=0;
    rLocalLHS(19,30)=0;
    rLocalLHS(19,31)=DOperator(2,2);
    rLocalLHS(19,32)=0;
    rLocalLHS(19,33)=0;
    rLocalLHS(19,34)=DOperator(3,2);
    rLocalLHS(19,35)=0;
    rLocalLHS(20,0)=0;
    rLocalLHS(20,1)=0;
    rLocalLHS(20,2)=0;
    rLocalLHS(20,3)=0;
    rLocalLHS(20,4)=0;
    rLocalLHS(20,5)=0;
    rLocalLHS(20,6)=0;
    rLocalLHS(20,7)=0;
    rLocalLHS(20,8)=0;
    rLocalLHS(20,9)=0;
    rLocalLHS(20,10)=0;
    rLocalLHS(20,11)=0;
    rLocalLHS(20,12)=0;
    rLocalLHS(20,13)=0;
    rLocalLHS(20,14)=0;
    rLocalLHS(20,15)=0;
    rLocalLHS(20,16)=0;
    rLocalLHS(20,17)=0;
    rLocalLHS(20,18)=0;
    rLocalLHS(20,19)=0;
    rLocalLHS(20,20)=0;
    rLocalLHS(20,21)=0;
    rLocalLHS(20,22)=0;
    rLocalLHS(20,23)=0;
    rLocalLHS(20,24)=0;
    rLocalLHS(20,25)=0;
    rLocalLHS(20,26)=DOperator(0,2);
    rLocalLHS(20,27)=0;
    rLocalLHS(20,28)=0;
    rLocalLHS(20,29)=DOperator(1,2);
    rLocalLHS(20,30)=0;
    rLocalLHS(20,31)=0;
    rLocalLHS(20,32)=DOperator(2,2);
    rLocalLHS(20,33)=0;
    rLocalLHS(20,34)=0;
    rLocalLHS(20,35)=DOperator(3,2);
    rLocalLHS(21,0)=0;
    rLocalLHS(21,1)=0;
    rLocalLHS(21,2)=0;
    rLocalLHS(21,3)=0;
    rLocalLHS(21,4)=0;
    rLocalLHS(21,5)=0;
    rLocalLHS(21,6)=0;
    rLocalLHS(21,7)=0;
    rLocalLHS(21,8)=0;
    rLocalLHS(21,9)=0;
    rLocalLHS(21,10)=0;
    rLocalLHS(21,11)=0;
    rLocalLHS(21,12)=0;
    rLocalLHS(21,13)=0;
    rLocalLHS(21,14)=0;
    rLocalLHS(21,15)=0;
    rLocalLHS(21,16)=0;
    rLocalLHS(21,17)=0;
    rLocalLHS(21,18)=0;
    rLocalLHS(21,19)=0;
    rLocalLHS(21,20)=0;
    rLocalLHS(21,21)=0;
    rLocalLHS(21,22)=0;
    rLocalLHS(21,23)=0;
    rLocalLHS(21,24)=DOperator(0,3);
    rLocalLHS(21,25)=0;
    rLocalLHS(21,26)=0;
    rLocalLHS(21,27)=DOperator(1,3);
    rLocalLHS(21,28)=0;
    rLocalLHS(21,29)=0;
    rLocalLHS(21,30)=DOperator(2,3);
    rLocalLHS(21,31)=0;
    rLocalLHS(21,32)=0;
    rLocalLHS(21,33)=DOperator(3,3);
    rLocalLHS(21,34)=0;
    rLocalLHS(21,35)=0;
    rLocalLHS(22,0)=0;
    rLocalLHS(22,1)=0;
    rLocalLHS(22,2)=0;
    rLocalLHS(22,3)=0;
    rLocalLHS(22,4)=0;
    rLocalLHS(22,5)=0;
    rLocalLHS(22,6)=0;
    rLocalLHS(22,7)=0;
    rLocalLHS(22,8)=0;
    rLocalLHS(22,9)=0;
    rLocalLHS(22,10)=0;
    rLocalLHS(22,11)=0;
    rLocalLHS(22,12)=0;
    rLocalLHS(22,13)=0;
    rLocalLHS(22,14)=0;
    rLocalLHS(22,15)=0;
    rLocalLHS(22,16)=0;
    rLocalLHS(22,17)=0;
    rLocalLHS(22,18)=0;
    rLocalLHS(22,19)=0;
    rLocalLHS(22,20)=0;
    rLocalLHS(22,21)=0;
    rLocalLHS(22,22)=0;
    rLocalLHS(22,23)=0;
    rLocalLHS(22,24)=0;
    rLocalLHS(22,25)=DOperator(0,3);
    rLocalLHS(22,26)=0;
    rLocalLHS(22,27)=0;
    rLocalLHS(22,28)=DOperator(1,3);
    rLocalLHS(22,29)=0;
    rLocalLHS(22,30)=0;
    rLocalLHS(22,31)=DOperator(2,3);
    rLocalLHS(22,32)=0;
    rLocalLHS(22,33)=0;
    rLocalLHS(22,34)=DOperator(3,3);
    rLocalLHS(22,35)=0;
    rLocalLHS(23,0)=0;
    rLocalLHS(23,1)=0;
    rLocalLHS(23,2)=0;
    rLocalLHS(23,3)=0;
    rLocalLHS(23,4)=0;
    rLocalLHS(23,5)=0;
    rLocalLHS(23,6)=0;
    rLocalLHS(23,7)=0;
    rLocalLHS(23,8)=0;
    rLocalLHS(23,9)=0;
    rLocalLHS(23,10)=0;
    rLocalLHS(23,11)=0;
    rLocalLHS(23,12)=0;
    rLocalLHS(23,13)=0;
    rLocalLHS(23,14)=0;
    rLocalLHS(23,15)=0;
    rLocalLHS(23,16)=0;
    rLocalLHS(23,17)=0;
    rLocalLHS(23,18)=0;
    rLocalLHS(23,19)=0;
    rLocalLHS(23,20)=0;
    rLocalLHS(23,21)=0;
    rLocalLHS(23,22)=0;
    rLocalLHS(23,23)=0;
    rLocalLHS(23,24)=0;
    rLocalLHS(23,25)=0;
    rLocalLHS(23,26)=DOperator(0,3);
    rLocalLHS(23,27)=0;
    rLocalLHS(23,28)=0;
    rLocalLHS(23,29)=DOperator(1,3);
    rLocalLHS(23,30)=0;
    rLocalLHS(23,31)=0;
    rLocalLHS(23,32)=DOperator(2,3);
    rLocalLHS(23,33)=0;
    rLocalLHS(23,34)=0;
    rLocalLHS(23,35)=DOperator(3,3);
    rLocalLHS(24,0)=clhs0;
    rLocalLHS(24,1)=0;
    rLocalLHS(24,2)=0;
    rLocalLHS(24,3)=clhs4;
    rLocalLHS(24,4)=0;
    rLocalLHS(24,5)=0;
    rLocalLHS(24,6)=clhs8;
    rLocalLHS(24,7)=0;
    rLocalLHS(24,8)=0;
    rLocalLHS(24,9)=clhs12;
    rLocalLHS(24,10)=0;
    rLocalLHS(24,11)=0;
    rLocalLHS(24,12)=DOperator(0,0);
    rLocalLHS(24,13)=0;
    rLocalLHS(24,14)=0;
    rLocalLHS(24,15)=DOperator(0,1);
    rLocalLHS(24,16)=0;
    rLocalLHS(24,17)=0;
    rLocalLHS(24,18)=DOperator(0,2);
    rLocalLHS(24,19)=0;
    rLocalLHS(24,20)=0;
    rLocalLHS(24,21)=DOperator(0,3);
    rLocalLHS(24,22)=0;
    rLocalLHS(24,23)=0;
    rLocalLHS(24,24)=0;
    rLocalLHS(24,25)=0;
    rLocalLHS(24,26)=0;
    rLocalLHS(24,27)=0;
    rLocalLHS(24,28)=0;
    rLocalLHS(24,29)=0;
    rLocalLHS(24,30)=0;
    rLocalLHS(24,31)=0;
    rLocalLHS(24,32)=0;
    rLocalLHS(24,33)=0;
    rLocalLHS(24,34)=0;
    rLocalLHS(24,35)=0;
    rLocalLHS(25,0)=0;
    rLocalLHS(25,1)=clhs0;
    rLocalLHS(25,2)=0;
    rLocalLHS(25,3)=0;
    rLocalLHS(25,4)=clhs4;
    rLocalLHS(25,5)=0;
    rLocalLHS(25,6)=0;
    rLocalLHS(25,7)=clhs8;
    rLocalLHS(25,8)=0;
    rLocalLHS(25,9)=0;
    rLocalLHS(25,10)=clhs12;
    rLocalLHS(25,11)=0;
    rLocalLHS(25,12)=0;
    rLocalLHS(25,13)=DOperator(0,0);
    rLocalLHS(25,14)=0;
    rLocalLHS(25,15)=0;
    rLocalLHS(25,16)=DOperator(0,1);
    rLocalLHS(25,17)=0;
    rLocalLHS(25,18)=0;
    rLocalLHS(25,19)=DOperator(0,2);
    rLocalLHS(25,20)=0;
    rLocalLHS(25,21)=0;
    rLocalLHS(25,22)=DOperator(0,3);
    rLocalLHS(25,23)=0;
    rLocalLHS(25,24)=0;
    rLocalLHS(25,25)=0;
    rLocalLHS(25,26)=0;
    rLocalLHS(25,27)=0;
    rLocalLHS(25,28)=0;
    rLocalLHS(25,29)=0;
    rLocalLHS(25,30)=0;
    rLocalLHS(25,31)=0;
    rLocalLHS(25,32)=0;
    rLocalLHS(25,33)=0;
    rLocalLHS(25,34)=0;
    rLocalLHS(25,35)=0;
    rLocalLHS(26,0)=0;
    rLocalLHS(26,1)=0;
    rLocalLHS(26,2)=clhs0;
    rLocalLHS(26,3)=0;
    rLocalLHS(26,4)=0;
    rLocalLHS(26,5)=clhs4;
    rLocalLHS(26,6)=0;
    rLocalLHS(26,7)=0;
    rLocalLHS(26,8)=clhs8;
    rLocalLHS(26,9)=0;
    rLocalLHS(26,10)=0;
    rLocalLHS(26,11)=clhs12;
    rLocalLHS(26,12)=0;
    rLocalLHS(26,13)=0;
    rLocalLHS(26,14)=DOperator(0,0);
    rLocalLHS(26,15)=0;
    rLocalLHS(26,16)=0;
    rLocalLHS(26,17)=DOperator(0,1);
    rLocalLHS(26,18)=0;
    rLocalLHS(26,19)=0;
    rLocalLHS(26,20)=DOperator(0,2);
    rLocalLHS(26,21)=0;
    rLocalLHS(26,22)=0;
    rLocalLHS(26,23)=DOperator(0,3);
    rLocalLHS(26,24)=0;
    rLocalLHS(26,25)=0;
    rLocalLHS(26,26)=0;
    rLocalLHS(26,27)=0;
    rLocalLHS(26,28)=0;
    rLocalLHS(26,29)=0;
    rLocalLHS(26,30)=0;
    rLocalLHS(26,31)=0;
    rLocalLHS(26,32)=0;
    rLocalLHS(26,33)=0;
    rLocalLHS(26,34)=0;
    rLocalLHS(26,35)=0;
    rLocalLHS(27,0)=clhs1;
    rLocalLHS(27,1)=0;
    rLocalLHS(27,2)=0;
    rLocalLHS(27,3)=clhs5;
    rLocalLHS(27,4)=0;
    rLocalLHS(27,5)=0;
    rLocalLHS(27,6)=clhs9;
    rLocalLHS(27,7)=0;
    rLocalLHS(27,8)=0;
    rLocalLHS(27,9)=clhs13;
    rLocalLHS(27,10)=0;
    rLocalLHS(27,11)=0;
    rLocalLHS(27,12)=DOperator(1,0);
    rLocalLHS(27,13)=0;
    rLocalLHS(27,14)=0;
    rLocalLHS(27,15)=DOperator(1,1);
    rLocalLHS(27,16)=0;
    rLocalLHS(27,17)=0;
    rLocalLHS(27,18)=DOperator(1,2);
    rLocalLHS(27,19)=0;
    rLocalLHS(27,20)=0;
    rLocalLHS(27,21)=DOperator(1,3);
    rLocalLHS(27,22)=0;
    rLocalLHS(27,23)=0;
    rLocalLHS(27,24)=0;
    rLocalLHS(27,25)=0;
    rLocalLHS(27,26)=0;
    rLocalLHS(27,27)=0;
    rLocalLHS(27,28)=0;
    rLocalLHS(27,29)=0;
    rLocalLHS(27,30)=0;
    rLocalLHS(27,31)=0;
    rLocalLHS(27,32)=0;
    rLocalLHS(27,33)=0;
    rLocalLHS(27,34)=0;
    rLocalLHS(27,35)=0;
    rLocalLHS(28,0)=0;
    rLocalLHS(28,1)=clhs1;
    rLocalLHS(28,2)=0;
    rLocalLHS(28,3)=0;
    rLocalLHS(28,4)=clhs5;
    rLocalLHS(28,5)=0;
    rLocalLHS(28,6)=0;
    rLocalLHS(28,7)=clhs9;
    rLocalLHS(28,8)=0;
    rLocalLHS(28,9)=0;
    rLocalLHS(28,10)=clhs13;
    rLocalLHS(28,11)=0;
    rLocalLHS(28,12)=0;
    rLocalLHS(28,13)=DOperator(1,0);
    rLocalLHS(28,14)=0;
    rLocalLHS(28,15)=0;
    rLocalLHS(28,16)=DOperator(1,1);
    rLocalLHS(28,17)=0;
    rLocalLHS(28,18)=0;
    rLocalLHS(28,19)=DOperator(1,2);
    rLocalLHS(28,20)=0;
    rLocalLHS(28,21)=0;
    rLocalLHS(28,22)=DOperator(1,3);
    rLocalLHS(28,23)=0;
    rLocalLHS(28,24)=0;
    rLocalLHS(28,25)=0;
    rLocalLHS(28,26)=0;
    rLocalLHS(28,27)=0;
    rLocalLHS(28,28)=0;
    rLocalLHS(28,29)=0;
    rLocalLHS(28,30)=0;
    rLocalLHS(28,31)=0;
    rLocalLHS(28,32)=0;
    rLocalLHS(28,33)=0;
    rLocalLHS(28,34)=0;
    rLocalLHS(28,35)=0;
    rLocalLHS(29,0)=0;
    rLocalLHS(29,1)=0;
    rLocalLHS(29,2)=clhs1;
    rLocalLHS(29,3)=0;
    rLocalLHS(29,4)=0;
    rLocalLHS(29,5)=clhs5;
    rLocalLHS(29,6)=0;
    rLocalLHS(29,7)=0;
    rLocalLHS(29,8)=clhs9;
    rLocalLHS(29,9)=0;
    rLocalLHS(29,10)=0;
    rLocalLHS(29,11)=clhs13;
    rLocalLHS(29,12)=0;
    rLocalLHS(29,13)=0;
    rLocalLHS(29,14)=DOperator(1,0);
    rLocalLHS(29,15)=0;
    rLocalLHS(29,16)=0;
    rLocalLHS(29,17)=DOperator(1,1);
    rLocalLHS(29,18)=0;
    rLocalLHS(29,19)=0;
    rLocalLHS(29,20)=DOperator(1,2);
    rLocalLHS(29,21)=0;
    rLocalLHS(29,22)=0;
    rLocalLHS(29,23)=DOperator(1,3);
    rLocalLHS(29,24)=0;
    rLocalLHS(29,25)=0;
    rLocalLHS(29,26)=0;
    rLocalLHS(29,27)=0;
    rLocalLHS(29,28)=0;
    rLocalLHS(29,29)=0;
    rLocalLHS(29,30)=0;
    rLocalLHS(29,31)=0;
    rLocalLHS(29,32)=0;
    rLocalLHS(29,33)=0;
    rLocalLHS(29,34)=0;
    rLocalLHS(29,35)=0;
    rLocalLHS(30,0)=clhs2;
    rLocalLHS(30,1)=0;
    rLocalLHS(30,2)=0;
    rLocalLHS(30,3)=clhs6;
    rLocalLHS(30,4)=0;
    rLocalLHS(30,5)=0;
    rLocalLHS(30,6)=clhs10;
    rLocalLHS(30,7)=0;
    rLocalLHS(30,8)=0;
    rLocalLHS(30,9)=clhs14;
    rLocalLHS(30,10)=0;
    rLocalLHS(30,11)=0;
    rLocalLHS(30,12)=DOperator(2,0);
    rLocalLHS(30,13)=0;
    rLocalLHS(30,14)=0;
    rLocalLHS(30,15)=DOperator(2,1);
    rLocalLHS(30,16)=0;
    rLocalLHS(30,17)=0;
    rLocalLHS(30,18)=DOperator(2,2);
    rLocalLHS(30,19)=0;
    rLocalLHS(30,20)=0;
    rLocalLHS(30,21)=DOperator(2,3);
    rLocalLHS(30,22)=0;
    rLocalLHS(30,23)=0;
    rLocalLHS(30,24)=0;
    rLocalLHS(30,25)=0;
    rLocalLHS(30,26)=0;
    rLocalLHS(30,27)=0;
    rLocalLHS(30,28)=0;
    rLocalLHS(30,29)=0;
    rLocalLHS(30,30)=0;
    rLocalLHS(30,31)=0;
    rLocalLHS(30,32)=0;
    rLocalLHS(30,33)=0;
    rLocalLHS(30,34)=0;
    rLocalLHS(30,35)=0;
    rLocalLHS(31,0)=0;
    rLocalLHS(31,1)=clhs2;
    rLocalLHS(31,2)=0;
    rLocalLHS(31,3)=0;
    rLocalLHS(31,4)=clhs6;
    rLocalLHS(31,5)=0;
    rLocalLHS(31,6)=0;
    rLocalLHS(31,7)=clhs10;
    rLocalLHS(31,8)=0;
    rLocalLHS(31,9)=0;
    rLocalLHS(31,10)=clhs14;
    rLocalLHS(31,11)=0;
    rLocalLHS(31,12)=0;
    rLocalLHS(31,13)=DOperator(2,0);
    rLocalLHS(31,14)=0;
    rLocalLHS(31,15)=0;
    rLocalLHS(31,16)=DOperator(2,1);
    rLocalLHS(31,17)=0;
    rLocalLHS(31,18)=0;
    rLocalLHS(31,19)=DOperator(2,2);
    rLocalLHS(31,20)=0;
    rLocalLHS(31,21)=0;
    rLocalLHS(31,22)=DOperator(2,3);
    rLocalLHS(31,23)=0;
    rLocalLHS(31,24)=0;
    rLocalLHS(31,25)=0;
    rLocalLHS(31,26)=0;
    rLocalLHS(31,27)=0;
    rLocalLHS(31,28)=0;
    rLocalLHS(31,29)=0;
    rLocalLHS(31,30)=0;
    rLocalLHS(31,31)=0;
    rLocalLHS(31,32)=0;
    rLocalLHS(31,33)=0;
    rLocalLHS(31,34)=0;
    rLocalLHS(31,35)=0;
    rLocalLHS(32,0)=0;
    rLocalLHS(32,1)=0;
    rLocalLHS(32,2)=clhs2;
    rLocalLHS(32,3)=0;
    rLocalLHS(32,4)=0;
    rLocalLHS(32,5)=clhs6;
    rLocalLHS(32,6)=0;
    rLocalLHS(32,7)=0;
    rLocalLHS(32,8)=clhs10;
    rLocalLHS(32,9)=0;
    rLocalLHS(32,10)=0;
    rLocalLHS(32,11)=clhs14;
    rLocalLHS(32,12)=0;
    rLocalLHS(32,13)=0;
    rLocalLHS(32,14)=DOperator(2,0);
    rLocalLHS(32,15)=0;
    rLocalLHS(32,16)=0;
    rLocalLHS(32,17)=DOperator(2,1);
    rLocalLHS(32,18)=0;
    rLocalLHS(32,19)=0;
    rLocalLHS(32,20)=DOperator(2,2);
    rLocalLHS(32,21)=0;
    rLocalLHS(32,22)=0;
    rLocalLHS(32,23)=DOperator(2,3);
    rLocalLHS(32,24)=0;
    rLocalLHS(32,25)=0;
    rLocalLHS(32,26)=0;
    rLocalLHS(32,27)=0;
    rLocalLHS(32,28)=0;
    rLocalLHS(32,29)=0;
    rLocalLHS(32,30)=0;
    rLocalLHS(32,31)=0;
    rLocalLHS(32,32)=0;
    rLocalLHS(32,33)=0;
    rLocalLHS(32,34)=0;
    rLocalLHS(32,35)=0;
    rLocalLHS(33,0)=clhs3;
    rLocalLHS(33,1)=0;
    rLocalLHS(33,2)=0;
    rLocalLHS(33,3)=clhs7;
    rLocalLHS(33,4)=0;
    rLocalLHS(33,5)=0;
    rLocalLHS(33,6)=clhs11;
    rLocalLHS(33,7)=0;
    rLocalLHS(33,8)=0;
    rLocalLHS(33,9)=clhs15;
    rLocalLHS(33,10)=0;
    rLocalLHS(33,11)=0;
    rLocalLHS(33,12)=DOperator(3,0);
    rLocalLHS(33,13)=0;
    rLocalLHS(33,14)=0;
    rLocalLHS(33,15)=DOperator(3,1);
    rLocalLHS(33,16)=0;
    rLocalLHS(33,17)=0;
    rLocalLHS(33,18)=DOperator(3,2);
    rLocalLHS(33,19)=0;
    rLocalLHS(33,20)=0;
    rLocalLHS(33,21)=DOperator(3,3);
    rLocalLHS(33,22)=0;
    rLocalLHS(33,23)=0;
    rLocalLHS(33,24)=0;
    rLocalLHS(33,25)=0;
    rLocalLHS(33,26)=0;
    rLocalLHS(33,27)=0;
    rLocalLHS(33,28)=0;
    rLocalLHS(33,29)=0;
    rLocalLHS(33,30)=0;
    rLocalLHS(33,31)=0;
    rLocalLHS(33,32)=0;
    rLocalLHS(33,33)=0;
    rLocalLHS(33,34)=0;
    rLocalLHS(33,35)=0;
    rLocalLHS(34,0)=0;
    rLocalLHS(34,1)=clhs3;
    rLocalLHS(34,2)=0;
    rLocalLHS(34,3)=0;
    rLocalLHS(34,4)=clhs7;
    rLocalLHS(34,5)=0;
    rLocalLHS(34,6)=0;
    rLocalLHS(34,7)=clhs11;
    rLocalLHS(34,8)=0;
    rLocalLHS(34,9)=0;
    rLocalLHS(34,10)=clhs15;
    rLocalLHS(34,11)=0;
    rLocalLHS(34,12)=0;
    rLocalLHS(34,13)=DOperator(3,0);
    rLocalLHS(34,14)=0;
    rLocalLHS(34,15)=0;
    rLocalLHS(34,16)=DOperator(3,1);
    rLocalLHS(34,17)=0;
    rLocalLHS(34,18)=0;
    rLocalLHS(34,19)=DOperator(3,2);
    rLocalLHS(34,20)=0;
    rLocalLHS(34,21)=0;
    rLocalLHS(34,22)=DOperator(3,3);
    rLocalLHS(34,23)=0;
    rLocalLHS(34,24)=0;
    rLocalLHS(34,25)=0;
    rLocalLHS(34,26)=0;
    rLocalLHS(34,27)=0;
    rLocalLHS(34,28)=0;
    rLocalLHS(34,29)=0;
    rLocalLHS(34,30)=0;
    rLocalLHS(34,31)=0;
    rLocalLHS(34,32)=0;
    rLocalLHS(34,33)=0;
    rLocalLHS(34,34)=0;
    rLocalLHS(34,35)=0;
    rLocalLHS(35,0)=0;
    rLocalLHS(35,1)=0;
    rLocalLHS(35,2)=clhs3;
    rLocalLHS(35,3)=0;
    rLocalLHS(35,4)=0;
    rLocalLHS(35,5)=clhs7;
    rLocalLHS(35,6)=0;
    rLocalLHS(35,7)=0;
    rLocalLHS(35,8)=clhs11;
    rLocalLHS(35,9)=0;
    rLocalLHS(35,10)=0;
    rLocalLHS(35,11)=clhs15;
    rLocalLHS(35,12)=0;
    rLocalLHS(35,13)=0;
    rLocalLHS(35,14)=DOperator(3,0);
    rLocalLHS(35,15)=0;
    rLocalLHS(35,16)=0;
    rLocalLHS(35,17)=DOperator(3,1);
    rLocalLHS(35,18)=0;
    rLocalLHS(35,19)=0;
    rLocalLHS(35,20)=DOperator(3,2);
    rLocalLHS(35,21)=0;
    rLocalLHS(35,22)=0;
    rLocalLHS(35,23)=DOperator(3,3);
    rLocalLHS(35,24)=0;
    rLocalLHS(35,25)=0;
    rLocalLHS(35,26)=0;
    rLocalLHS(35,27)=0;
    rLocalLHS(35,28)=0;
    rLocalLHS(35,29)=0;
    rLocalLHS(35,30)=0;
    rLocalLHS(35,31)=0;
    rLocalLHS(35,32)=0;
    rLocalLHS(35,33)=0;
    rLocalLHS(35,34)=0;
    rLocalLHS(35,35)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,4,8>::CalculateLocalLHS<MeshTyingMortarCondition<3,4,8>::ScalarValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 3, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(0,1);
    const double clhs4 =     -MOperator(1,1);
    const double clhs5 =     -MOperator(2,1);
    const double clhs6 =     -MOperator(0,2);
    const double clhs7 =     -MOperator(1,2);
    const double clhs8 =     -MOperator(2,2);
    const double clhs9 =     -MOperator(0,3);
    const double clhs10 =     -MOperator(1,3);
    const double clhs11 =     -MOperator(2,3);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=clhs0;
    rLocalLHS(0,8)=clhs1;
    rLocalLHS(0,9)=clhs2;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=clhs3;
    rLocalLHS(1,8)=clhs4;
    rLocalLHS(1,9)=clhs5;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=clhs6;
    rLocalLHS(2,8)=clhs7;
    rLocalLHS(2,9)=clhs8;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=clhs9;
    rLocalLHS(3,8)=clhs10;
    rLocalLHS(3,9)=clhs11;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=DOperator(0,0);
    rLocalLHS(4,8)=DOperator(1,0);
    rLocalLHS(4,9)=DOperator(2,0);
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=DOperator(0,1);
    rLocalLHS(5,8)=DOperator(1,1);
    rLocalLHS(5,9)=DOperator(2,1);
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=DOperator(0,2);
    rLocalLHS(6,8)=DOperator(1,2);
    rLocalLHS(6,9)=DOperator(2,2);
    rLocalLHS(7,0)=clhs0;
    rLocalLHS(7,1)=clhs3;
    rLocalLHS(7,2)=clhs6;
    rLocalLHS(7,3)=clhs9;
    rLocalLHS(7,4)=DOperator(0,0);
    rLocalLHS(7,5)=DOperator(0,1);
    rLocalLHS(7,6)=DOperator(0,2);
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=0;
    rLocalLHS(8,0)=clhs1;
    rLocalLHS(8,1)=clhs4;
    rLocalLHS(8,2)=clhs7;
    rLocalLHS(8,3)=clhs10;
    rLocalLHS(8,4)=DOperator(1,0);
    rLocalLHS(8,5)=DOperator(1,1);
    rLocalLHS(8,6)=DOperator(1,2);
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(9,0)=clhs2;
    rLocalLHS(9,1)=clhs5;
    rLocalLHS(9,2)=clhs8;
    rLocalLHS(9,3)=clhs11;
    rLocalLHS(9,4)=DOperator(2,0);
    rLocalLHS(9,5)=DOperator(2,1);
    rLocalLHS(9,6)=DOperator(2,2);
    rLocalLHS(9,7)=0;
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,4,8>::CalculateLocalLHS<MeshTyingMortarCondition<3,4,8>::Vector3DValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 3, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(0,1);
    const double clhs4 =     -MOperator(1,1);
    const double clhs5 =     -MOperator(2,1);
    const double clhs6 =     -MOperator(0,2);
    const double clhs7 =     -MOperator(1,2);
    const double clhs8 =     -MOperator(2,2);
    const double clhs9 =     -MOperator(0,3);
    const double clhs10 =     -MOperator(1,3);
    const double clhs11 =     -MOperator(2,3);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=0;
    rLocalLHS(0,8)=0;
    rLocalLHS(0,9)=0;
    rLocalLHS(0,10)=0;
    rLocalLHS(0,11)=0;
    rLocalLHS(0,12)=0;
    rLocalLHS(0,13)=0;
    rLocalLHS(0,14)=0;
    rLocalLHS(0,15)=0;
    rLocalLHS(0,16)=0;
    rLocalLHS(0,17)=0;
    rLocalLHS(0,18)=0;
    rLocalLHS(0,19)=0;
    rLocalLHS(0,20)=0;
    rLocalLHS(0,21)=clhs0;
    rLocalLHS(0,22)=0;
    rLocalLHS(0,23)=0;
    rLocalLHS(0,24)=clhs1;
    rLocalLHS(0,25)=0;
    rLocalLHS(0,26)=0;
    rLocalLHS(0,27)=clhs2;
    rLocalLHS(0,28)=0;
    rLocalLHS(0,29)=0;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=0;
    rLocalLHS(1,8)=0;
    rLocalLHS(1,9)=0;
    rLocalLHS(1,10)=0;
    rLocalLHS(1,11)=0;
    rLocalLHS(1,12)=0;
    rLocalLHS(1,13)=0;
    rLocalLHS(1,14)=0;
    rLocalLHS(1,15)=0;
    rLocalLHS(1,16)=0;
    rLocalLHS(1,17)=0;
    rLocalLHS(1,18)=0;
    rLocalLHS(1,19)=0;
    rLocalLHS(1,20)=0;
    rLocalLHS(1,21)=0;
    rLocalLHS(1,22)=clhs0;
    rLocalLHS(1,23)=0;
    rLocalLHS(1,24)=0;
    rLocalLHS(1,25)=clhs1;
    rLocalLHS(1,26)=0;
    rLocalLHS(1,27)=0;
    rLocalLHS(1,28)=clhs2;
    rLocalLHS(1,29)=0;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=0;
    rLocalLHS(2,8)=0;
    rLocalLHS(2,9)=0;
    rLocalLHS(2,10)=0;
    rLocalLHS(2,11)=0;
    rLocalLHS(2,12)=0;
    rLocalLHS(2,13)=0;
    rLocalLHS(2,14)=0;
    rLocalLHS(2,15)=0;
    rLocalLHS(2,16)=0;
    rLocalLHS(2,17)=0;
    rLocalLHS(2,18)=0;
    rLocalLHS(2,19)=0;
    rLocalLHS(2,20)=0;
    rLocalLHS(2,21)=0;
    rLocalLHS(2,22)=0;
    rLocalLHS(2,23)=clhs0;
    rLocalLHS(2,24)=0;
    rLocalLHS(2,25)=0;
    rLocalLHS(2,26)=clhs1;
    rLocalLHS(2,27)=0;
    rLocalLHS(2,28)=0;
    rLocalLHS(2,29)=clhs2;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=0;
    rLocalLHS(3,8)=0;
    rLocalLHS(3,9)=0;
    rLocalLHS(3,10)=0;
    rLocalLHS(3,11)=0;
    rLocalLHS(3,12)=0;
    rLocalLHS(3,13)=0;
    rLocalLHS(3,14)=0;
    rLocalLHS(3,15)=0;
    rLocalLHS(3,16)=0;
    rLocalLHS(3,17)=0;
    rLocalLHS(3,18)=0;
    rLocalLHS(3,19)=0;
    rLocalLHS(3,20)=0;
    rLocalLHS(3,21)=clhs3;
    rLocalLHS(3,22)=0;
    rLocalLHS(3,23)=0;
    rLocalLHS(3,24)=clhs4;
    rLocalLHS(3,25)=0;
    rLocalLHS(3,26)=0;
    rLocalLHS(3,27)=clhs5;
    rLocalLHS(3,28)=0;
    rLocalLHS(3,29)=0;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=0;
    rLocalLHS(4,8)=0;
    rLocalLHS(4,9)=0;
    rLocalLHS(4,10)=0;
    rLocalLHS(4,11)=0;
    rLocalLHS(4,12)=0;
    rLocalLHS(4,13)=0;
    rLocalLHS(4,14)=0;
    rLocalLHS(4,15)=0;
    rLocalLHS(4,16)=0;
    rLocalLHS(4,17)=0;
    rLocalLHS(4,18)=0;
    rLocalLHS(4,19)=0;
    rLocalLHS(4,20)=0;
    rLocalLHS(4,21)=0;
    rLocalLHS(4,22)=clhs3;
    rLocalLHS(4,23)=0;
    rLocalLHS(4,24)=0;
    rLocalLHS(4,25)=clhs4;
    rLocalLHS(4,26)=0;
    rLocalLHS(4,27)=0;
    rLocalLHS(4,28)=clhs5;
    rLocalLHS(4,29)=0;
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=0;
    rLocalLHS(5,8)=0;
    rLocalLHS(5,9)=0;
    rLocalLHS(5,10)=0;
    rLocalLHS(5,11)=0;
    rLocalLHS(5,12)=0;
    rLocalLHS(5,13)=0;
    rLocalLHS(5,14)=0;
    rLocalLHS(5,15)=0;
    rLocalLHS(5,16)=0;
    rLocalLHS(5,17)=0;
    rLocalLHS(5,18)=0;
    rLocalLHS(5,19)=0;
    rLocalLHS(5,20)=0;
    rLocalLHS(5,21)=0;
    rLocalLHS(5,22)=0;
    rLocalLHS(5,23)=clhs3;
    rLocalLHS(5,24)=0;
    rLocalLHS(5,25)=0;
    rLocalLHS(5,26)=clhs4;
    rLocalLHS(5,27)=0;
    rLocalLHS(5,28)=0;
    rLocalLHS(5,29)=clhs5;
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=0;
    rLocalLHS(6,9)=0;
    rLocalLHS(6,10)=0;
    rLocalLHS(6,11)=0;
    rLocalLHS(6,12)=0;
    rLocalLHS(6,13)=0;
    rLocalLHS(6,14)=0;
    rLocalLHS(6,15)=0;
    rLocalLHS(6,16)=0;
    rLocalLHS(6,17)=0;
    rLocalLHS(6,18)=0;
    rLocalLHS(6,19)=0;
    rLocalLHS(6,20)=0;
    rLocalLHS(6,21)=clhs6;
    rLocalLHS(6,22)=0;
    rLocalLHS(6,23)=0;
    rLocalLHS(6,24)=clhs7;
    rLocalLHS(6,25)=0;
    rLocalLHS(6,26)=0;
    rLocalLHS(6,27)=clhs8;
    rLocalLHS(6,28)=0;
    rLocalLHS(6,29)=0;
    rLocalLHS(7,0)=0;
    rLocalLHS(7,1)=0;
    rLocalLHS(7,2)=0;
    rLocalLHS(7,3)=0;
    rLocalLHS(7,4)=0;
    rLocalLHS(7,5)=0;
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=0;
    rLocalLHS(7,10)=0;
    rLocalLHS(7,11)=0;
    rLocalLHS(7,12)=0;
    rLocalLHS(7,13)=0;
    rLocalLHS(7,14)=0;
    rLocalLHS(7,15)=0;
    rLocalLHS(7,16)=0;
    rLocalLHS(7,17)=0;
    rLocalLHS(7,18)=0;
    rLocalLHS(7,19)=0;
    rLocalLHS(7,20)=0;
    rLocalLHS(7,21)=0;
    rLocalLHS(7,22)=clhs6;
    rLocalLHS(7,23)=0;
    rLocalLHS(7,24)=0;
    rLocalLHS(7,25)=clhs7;
    rLocalLHS(7,26)=0;
    rLocalLHS(7,27)=0;
    rLocalLHS(7,28)=clhs8;
    rLocalLHS(7,29)=0;
    rLocalLHS(8,0)=0;
    rLocalLHS(8,1)=0;
    rLocalLHS(8,2)=0;
    rLocalLHS(8,3)=0;
    rLocalLHS(8,4)=0;
    rLocalLHS(8,5)=0;
    rLocalLHS(8,6)=0;
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(8,11)=0;
    rLocalLHS(8,12)=0;
    rLocalLHS(8,13)=0;
    rLocalLHS(8,14)=0;
    rLocalLHS(8,15)=0;
    rLocalLHS(8,16)=0;
    rLocalLHS(8,17)=0;
    rLocalLHS(8,18)=0;
    rLocalLHS(8,19)=0;
    rLocalLHS(8,20)=0;
    rLocalLHS(8,21)=0;
    rLocalLHS(8,22)=0;
    rLocalLHS(8,23)=clhs6;
    rLocalLHS(8,24)=0;
    rLocalLHS(8,25)=0;
    rLocalLHS(8,26)=clhs7;
    rLocalLHS(8,27)=0;
    rLocalLHS(8,28)=0;
    rLocalLHS(8,29)=clhs8;
    rLocalLHS(9,0)=0;
    rLocalLHS(9,1)=0;
    rLocalLHS(9,2)=0;
    rLocalLHS(9,3)=0;
    rLocalLHS(9,4)=0;
    rLocalLHS(9,5)=0;
    rLocalLHS(9,6)=0;
    rLocalLHS(9,7)=0;
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(9,11)=0;
    rLocalLHS(9,12)=0;
    rLocalLHS(9,13)=0;
    rLocalLHS(9,14)=0;
    rLocalLHS(9,15)=0;
    rLocalLHS(9,16)=0;
    rLocalLHS(9,17)=0;
    rLocalLHS(9,18)=0;
    rLocalLHS(9,19)=0;
    rLocalLHS(9,20)=0;
    rLocalLHS(9,21)=clhs9;
    rLocalLHS(9,22)=0;
    rLocalLHS(9,23)=0;
    rLocalLHS(9,24)=clhs10;
    rLocalLHS(9,25)=0;
    rLocalLHS(9,26)=0;
    rLocalLHS(9,27)=clhs11;
    rLocalLHS(9,28)=0;
    rLocalLHS(9,29)=0;
    rLocalLHS(10,0)=0;
    rLocalLHS(10,1)=0;
    rLocalLHS(10,2)=0;
    rLocalLHS(10,3)=0;
    rLocalLHS(10,4)=0;
    rLocalLHS(10,5)=0;
    rLocalLHS(10,6)=0;
    rLocalLHS(10,7)=0;
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;
    rLocalLHS(10,11)=0;
    rLocalLHS(10,12)=0;
    rLocalLHS(10,13)=0;
    rLocalLHS(10,14)=0;
    rLocalLHS(10,15)=0;
    rLocalLHS(10,16)=0;
    rLocalLHS(10,17)=0;
    rLocalLHS(10,18)=0;
    rLocalLHS(10,19)=0;
    rLocalLHS(10,20)=0;
    rLocalLHS(10,21)=0;
    rLocalLHS(10,22)=clhs9;
    rLocalLHS(10,23)=0;
    rLocalLHS(10,24)=0;
    rLocalLHS(10,25)=clhs10;
    rLocalLHS(10,26)=0;
    rLocalLHS(10,27)=0;
    rLocalLHS(10,28)=clhs11;
    rLocalLHS(10,29)=0;
    rLocalLHS(11,0)=0;
    rLocalLHS(11,1)=0;
    rLocalLHS(11,2)=0;
    rLocalLHS(11,3)=0;
    rLocalLHS(11,4)=0;
    rLocalLHS(11,5)=0;
    rLocalLHS(11,6)=0;
    rLocalLHS(11,7)=0;
    rLocalLHS(11,8)=0;
    rLocalLHS(11,9)=0;
    rLocalLHS(11,10)=0;
    rLocalLHS(11,11)=0;
    rLocalLHS(11,12)=0;
    rLocalLHS(11,13)=0;
    rLocalLHS(11,14)=0;
    rLocalLHS(11,15)=0;
    rLocalLHS(11,16)=0;
    rLocalLHS(11,17)=0;
    rLocalLHS(11,18)=0;
    rLocalLHS(11,19)=0;
    rLocalLHS(11,20)=0;
    rLocalLHS(11,21)=0;
    rLocalLHS(11,22)=0;
    rLocalLHS(11,23)=clhs9;
    rLocalLHS(11,24)=0;
    rLocalLHS(11,25)=0;
    rLocalLHS(11,26)=clhs10;
    rLocalLHS(11,27)=0;
    rLocalLHS(11,28)=0;
    rLocalLHS(11,29)=clhs11;
    rLocalLHS(12,0)=0;
    rLocalLHS(12,1)=0;
    rLocalLHS(12,2)=0;
    rLocalLHS(12,3)=0;
    rLocalLHS(12,4)=0;
    rLocalLHS(12,5)=0;
    rLocalLHS(12,6)=0;
    rLocalLHS(12,7)=0;
    rLocalLHS(12,8)=0;
    rLocalLHS(12,9)=0;
    rLocalLHS(12,10)=0;
    rLocalLHS(12,11)=0;
    rLocalLHS(12,12)=0;
    rLocalLHS(12,13)=0;
    rLocalLHS(12,14)=0;
    rLocalLHS(12,15)=0;
    rLocalLHS(12,16)=0;
    rLocalLHS(12,17)=0;
    rLocalLHS(12,18)=0;
    rLocalLHS(12,19)=0;
    rLocalLHS(12,20)=0;
    rLocalLHS(12,21)=DOperator(0,0);
    rLocalLHS(12,22)=0;
    rLocalLHS(12,23)=0;
    rLocalLHS(12,24)=DOperator(1,0);
    rLocalLHS(12,25)=0;
    rLocalLHS(12,26)=0;
    rLocalLHS(12,27)=DOperator(2,0);
    rLocalLHS(12,28)=0;
    rLocalLHS(12,29)=0;
    rLocalLHS(13,0)=0;
    rLocalLHS(13,1)=0;
    rLocalLHS(13,2)=0;
    rLocalLHS(13,3)=0;
    rLocalLHS(13,4)=0;
    rLocalLHS(13,5)=0;
    rLocalLHS(13,6)=0;
    rLocalLHS(13,7)=0;
    rLocalLHS(13,8)=0;
    rLocalLHS(13,9)=0;
    rLocalLHS(13,10)=0;
    rLocalLHS(13,11)=0;
    rLocalLHS(13,12)=0;
    rLocalLHS(13,13)=0;
    rLocalLHS(13,14)=0;
    rLocalLHS(13,15)=0;
    rLocalLHS(13,16)=0;
    rLocalLHS(13,17)=0;
    rLocalLHS(13,18)=0;
    rLocalLHS(13,19)=0;
    rLocalLHS(13,20)=0;
    rLocalLHS(13,21)=0;
    rLocalLHS(13,22)=DOperator(0,0);
    rLocalLHS(13,23)=0;
    rLocalLHS(13,24)=0;
    rLocalLHS(13,25)=DOperator(1,0);
    rLocalLHS(13,26)=0;
    rLocalLHS(13,27)=0;
    rLocalLHS(13,28)=DOperator(2,0);
    rLocalLHS(13,29)=0;
    rLocalLHS(14,0)=0;
    rLocalLHS(14,1)=0;
    rLocalLHS(14,2)=0;
    rLocalLHS(14,3)=0;
    rLocalLHS(14,4)=0;
    rLocalLHS(14,5)=0;
    rLocalLHS(14,6)=0;
    rLocalLHS(14,7)=0;
    rLocalLHS(14,8)=0;
    rLocalLHS(14,9)=0;
    rLocalLHS(14,10)=0;
    rLocalLHS(14,11)=0;
    rLocalLHS(14,12)=0;
    rLocalLHS(14,13)=0;
    rLocalLHS(14,14)=0;
    rLocalLHS(14,15)=0;
    rLocalLHS(14,16)=0;
    rLocalLHS(14,17)=0;
    rLocalLHS(14,18)=0;
    rLocalLHS(14,19)=0;
    rLocalLHS(14,20)=0;
    rLocalLHS(14,21)=0;
    rLocalLHS(14,22)=0;
    rLocalLHS(14,23)=DOperator(0,0);
    rLocalLHS(14,24)=0;
    rLocalLHS(14,25)=0;
    rLocalLHS(14,26)=DOperator(1,0);
    rLocalLHS(14,27)=0;
    rLocalLHS(14,28)=0;
    rLocalLHS(14,29)=DOperator(2,0);
    rLocalLHS(15,0)=0;
    rLocalLHS(15,1)=0;
    rLocalLHS(15,2)=0;
    rLocalLHS(15,3)=0;
    rLocalLHS(15,4)=0;
    rLocalLHS(15,5)=0;
    rLocalLHS(15,6)=0;
    rLocalLHS(15,7)=0;
    rLocalLHS(15,8)=0;
    rLocalLHS(15,9)=0;
    rLocalLHS(15,10)=0;
    rLocalLHS(15,11)=0;
    rLocalLHS(15,12)=0;
    rLocalLHS(15,13)=0;
    rLocalLHS(15,14)=0;
    rLocalLHS(15,15)=0;
    rLocalLHS(15,16)=0;
    rLocalLHS(15,17)=0;
    rLocalLHS(15,18)=0;
    rLocalLHS(15,19)=0;
    rLocalLHS(15,20)=0;
    rLocalLHS(15,21)=DOperator(0,1);
    rLocalLHS(15,22)=0;
    rLocalLHS(15,23)=0;
    rLocalLHS(15,24)=DOperator(1,1);
    rLocalLHS(15,25)=0;
    rLocalLHS(15,26)=0;
    rLocalLHS(15,27)=DOperator(2,1);
    rLocalLHS(15,28)=0;
    rLocalLHS(15,29)=0;
    rLocalLHS(16,0)=0;
    rLocalLHS(16,1)=0;
    rLocalLHS(16,2)=0;
    rLocalLHS(16,3)=0;
    rLocalLHS(16,4)=0;
    rLocalLHS(16,5)=0;
    rLocalLHS(16,6)=0;
    rLocalLHS(16,7)=0;
    rLocalLHS(16,8)=0;
    rLocalLHS(16,9)=0;
    rLocalLHS(16,10)=0;
    rLocalLHS(16,11)=0;
    rLocalLHS(16,12)=0;
    rLocalLHS(16,13)=0;
    rLocalLHS(16,14)=0;
    rLocalLHS(16,15)=0;
    rLocalLHS(16,16)=0;
    rLocalLHS(16,17)=0;
    rLocalLHS(16,18)=0;
    rLocalLHS(16,19)=0;
    rLocalLHS(16,20)=0;
    rLocalLHS(16,21)=0;
    rLocalLHS(16,22)=DOperator(0,1);
    rLocalLHS(16,23)=0;
    rLocalLHS(16,24)=0;
    rLocalLHS(16,25)=DOperator(1,1);
    rLocalLHS(16,26)=0;
    rLocalLHS(16,27)=0;
    rLocalLHS(16,28)=DOperator(2,1);
    rLocalLHS(16,29)=0;
    rLocalLHS(17,0)=0;
    rLocalLHS(17,1)=0;
    rLocalLHS(17,2)=0;
    rLocalLHS(17,3)=0;
    rLocalLHS(17,4)=0;
    rLocalLHS(17,5)=0;
    rLocalLHS(17,6)=0;
    rLocalLHS(17,7)=0;
    rLocalLHS(17,8)=0;
    rLocalLHS(17,9)=0;
    rLocalLHS(17,10)=0;
    rLocalLHS(17,11)=0;
    rLocalLHS(17,12)=0;
    rLocalLHS(17,13)=0;
    rLocalLHS(17,14)=0;
    rLocalLHS(17,15)=0;
    rLocalLHS(17,16)=0;
    rLocalLHS(17,17)=0;
    rLocalLHS(17,18)=0;
    rLocalLHS(17,19)=0;
    rLocalLHS(17,20)=0;
    rLocalLHS(17,21)=0;
    rLocalLHS(17,22)=0;
    rLocalLHS(17,23)=DOperator(0,1);
    rLocalLHS(17,24)=0;
    rLocalLHS(17,25)=0;
    rLocalLHS(17,26)=DOperator(1,1);
    rLocalLHS(17,27)=0;
    rLocalLHS(17,28)=0;
    rLocalLHS(17,29)=DOperator(2,1);
    rLocalLHS(18,0)=0;
    rLocalLHS(18,1)=0;
    rLocalLHS(18,2)=0;
    rLocalLHS(18,3)=0;
    rLocalLHS(18,4)=0;
    rLocalLHS(18,5)=0;
    rLocalLHS(18,6)=0;
    rLocalLHS(18,7)=0;
    rLocalLHS(18,8)=0;
    rLocalLHS(18,9)=0;
    rLocalLHS(18,10)=0;
    rLocalLHS(18,11)=0;
    rLocalLHS(18,12)=0;
    rLocalLHS(18,13)=0;
    rLocalLHS(18,14)=0;
    rLocalLHS(18,15)=0;
    rLocalLHS(18,16)=0;
    rLocalLHS(18,17)=0;
    rLocalLHS(18,18)=0;
    rLocalLHS(18,19)=0;
    rLocalLHS(18,20)=0;
    rLocalLHS(18,21)=DOperator(0,2);
    rLocalLHS(18,22)=0;
    rLocalLHS(18,23)=0;
    rLocalLHS(18,24)=DOperator(1,2);
    rLocalLHS(18,25)=0;
    rLocalLHS(18,26)=0;
    rLocalLHS(18,27)=DOperator(2,2);
    rLocalLHS(18,28)=0;
    rLocalLHS(18,29)=0;
    rLocalLHS(19,0)=0;
    rLocalLHS(19,1)=0;
    rLocalLHS(19,2)=0;
    rLocalLHS(19,3)=0;
    rLocalLHS(19,4)=0;
    rLocalLHS(19,5)=0;
    rLocalLHS(19,6)=0;
    rLocalLHS(19,7)=0;
    rLocalLHS(19,8)=0;
    rLocalLHS(19,9)=0;
    rLocalLHS(19,10)=0;
    rLocalLHS(19,11)=0;
    rLocalLHS(19,12)=0;
    rLocalLHS(19,13)=0;
    rLocalLHS(19,14)=0;
    rLocalLHS(19,15)=0;
    rLocalLHS(19,16)=0;
    rLocalLHS(19,17)=0;
    rLocalLHS(19,18)=0;
    rLocalLHS(19,19)=0;
    rLocalLHS(19,20)=0;
    rLocalLHS(19,21)=0;
    rLocalLHS(19,22)=DOperator(0,2);
    rLocalLHS(19,23)=0;
    rLocalLHS(19,24)=0;
    rLocalLHS(19,25)=DOperator(1,2);
    rLocalLHS(19,26)=0;
    rLocalLHS(19,27)=0;
    rLocalLHS(19,28)=DOperator(2,2);
    rLocalLHS(19,29)=0;
    rLocalLHS(20,0)=0;
    rLocalLHS(20,1)=0;
    rLocalLHS(20,2)=0;
    rLocalLHS(20,3)=0;
    rLocalLHS(20,4)=0;
    rLocalLHS(20,5)=0;
    rLocalLHS(20,6)=0;
    rLocalLHS(20,7)=0;
    rLocalLHS(20,8)=0;
    rLocalLHS(20,9)=0;
    rLocalLHS(20,10)=0;
    rLocalLHS(20,11)=0;
    rLocalLHS(20,12)=0;
    rLocalLHS(20,13)=0;
    rLocalLHS(20,14)=0;
    rLocalLHS(20,15)=0;
    rLocalLHS(20,16)=0;
    rLocalLHS(20,17)=0;
    rLocalLHS(20,18)=0;
    rLocalLHS(20,19)=0;
    rLocalLHS(20,20)=0;
    rLocalLHS(20,21)=0;
    rLocalLHS(20,22)=0;
    rLocalLHS(20,23)=DOperator(0,2);
    rLocalLHS(20,24)=0;
    rLocalLHS(20,25)=0;
    rLocalLHS(20,26)=DOperator(1,2);
    rLocalLHS(20,27)=0;
    rLocalLHS(20,28)=0;
    rLocalLHS(20,29)=DOperator(2,2);
    rLocalLHS(21,0)=clhs0;
    rLocalLHS(21,1)=0;
    rLocalLHS(21,2)=0;
    rLocalLHS(21,3)=clhs3;
    rLocalLHS(21,4)=0;
    rLocalLHS(21,5)=0;
    rLocalLHS(21,6)=clhs6;
    rLocalLHS(21,7)=0;
    rLocalLHS(21,8)=0;
    rLocalLHS(21,9)=clhs9;
    rLocalLHS(21,10)=0;
    rLocalLHS(21,11)=0;
    rLocalLHS(21,12)=DOperator(0,0);
    rLocalLHS(21,13)=0;
    rLocalLHS(21,14)=0;
    rLocalLHS(21,15)=DOperator(0,1);
    rLocalLHS(21,16)=0;
    rLocalLHS(21,17)=0;
    rLocalLHS(21,18)=DOperator(0,2);
    rLocalLHS(21,19)=0;
    rLocalLHS(21,20)=0;
    rLocalLHS(21,21)=0;
    rLocalLHS(21,22)=0;
    rLocalLHS(21,23)=0;
    rLocalLHS(21,24)=0;
    rLocalLHS(21,25)=0;
    rLocalLHS(21,26)=0;
    rLocalLHS(21,27)=0;
    rLocalLHS(21,28)=0;
    rLocalLHS(21,29)=0;
    rLocalLHS(22,0)=0;
    rLocalLHS(22,1)=clhs0;
    rLocalLHS(22,2)=0;
    rLocalLHS(22,3)=0;
    rLocalLHS(22,4)=clhs3;
    rLocalLHS(22,5)=0;
    rLocalLHS(22,6)=0;
    rLocalLHS(22,7)=clhs6;
    rLocalLHS(22,8)=0;
    rLocalLHS(22,9)=0;
    rLocalLHS(22,10)=clhs9;
    rLocalLHS(22,11)=0;
    rLocalLHS(22,12)=0;
    rLocalLHS(22,13)=DOperator(0,0);
    rLocalLHS(22,14)=0;
    rLocalLHS(22,15)=0;
    rLocalLHS(22,16)=DOperator(0,1);
    rLocalLHS(22,17)=0;
    rLocalLHS(22,18)=0;
    rLocalLHS(22,19)=DOperator(0,2);
    rLocalLHS(22,20)=0;
    rLocalLHS(22,21)=0;
    rLocalLHS(22,22)=0;
    rLocalLHS(22,23)=0;
    rLocalLHS(22,24)=0;
    rLocalLHS(22,25)=0;
    rLocalLHS(22,26)=0;
    rLocalLHS(22,27)=0;
    rLocalLHS(22,28)=0;
    rLocalLHS(22,29)=0;
    rLocalLHS(23,0)=0;
    rLocalLHS(23,1)=0;
    rLocalLHS(23,2)=clhs0;
    rLocalLHS(23,3)=0;
    rLocalLHS(23,4)=0;
    rLocalLHS(23,5)=clhs3;
    rLocalLHS(23,6)=0;
    rLocalLHS(23,7)=0;
    rLocalLHS(23,8)=clhs6;
    rLocalLHS(23,9)=0;
    rLocalLHS(23,10)=0;
    rLocalLHS(23,11)=clhs9;
    rLocalLHS(23,12)=0;
    rLocalLHS(23,13)=0;
    rLocalLHS(23,14)=DOperator(0,0);
    rLocalLHS(23,15)=0;
    rLocalLHS(23,16)=0;
    rLocalLHS(23,17)=DOperator(0,1);
    rLocalLHS(23,18)=0;
    rLocalLHS(23,19)=0;
    rLocalLHS(23,20)=DOperator(0,2);
    rLocalLHS(23,21)=0;
    rLocalLHS(23,22)=0;
    rLocalLHS(23,23)=0;
    rLocalLHS(23,24)=0;
    rLocalLHS(23,25)=0;
    rLocalLHS(23,26)=0;
    rLocalLHS(23,27)=0;
    rLocalLHS(23,28)=0;
    rLocalLHS(23,29)=0;
    rLocalLHS(24,0)=clhs1;
    rLocalLHS(24,1)=0;
    rLocalLHS(24,2)=0;
    rLocalLHS(24,3)=clhs4;
    rLocalLHS(24,4)=0;
    rLocalLHS(24,5)=0;
    rLocalLHS(24,6)=clhs7;
    rLocalLHS(24,7)=0;
    rLocalLHS(24,8)=0;
    rLocalLHS(24,9)=clhs10;
    rLocalLHS(24,10)=0;
    rLocalLHS(24,11)=0;
    rLocalLHS(24,12)=DOperator(1,0);
    rLocalLHS(24,13)=0;
    rLocalLHS(24,14)=0;
    rLocalLHS(24,15)=DOperator(1,1);
    rLocalLHS(24,16)=0;
    rLocalLHS(24,17)=0;
    rLocalLHS(24,18)=DOperator(1,2);
    rLocalLHS(24,19)=0;
    rLocalLHS(24,20)=0;
    rLocalLHS(24,21)=0;
    rLocalLHS(24,22)=0;
    rLocalLHS(24,23)=0;
    rLocalLHS(24,24)=0;
    rLocalLHS(24,25)=0;
    rLocalLHS(24,26)=0;
    rLocalLHS(24,27)=0;
    rLocalLHS(24,28)=0;
    rLocalLHS(24,29)=0;
    rLocalLHS(25,0)=0;
    rLocalLHS(25,1)=clhs1;
    rLocalLHS(25,2)=0;
    rLocalLHS(25,3)=0;
    rLocalLHS(25,4)=clhs4;
    rLocalLHS(25,5)=0;
    rLocalLHS(25,6)=0;
    rLocalLHS(25,7)=clhs7;
    rLocalLHS(25,8)=0;
    rLocalLHS(25,9)=0;
    rLocalLHS(25,10)=clhs10;
    rLocalLHS(25,11)=0;
    rLocalLHS(25,12)=0;
    rLocalLHS(25,13)=DOperator(1,0);
    rLocalLHS(25,14)=0;
    rLocalLHS(25,15)=0;
    rLocalLHS(25,16)=DOperator(1,1);
    rLocalLHS(25,17)=0;
    rLocalLHS(25,18)=0;
    rLocalLHS(25,19)=DOperator(1,2);
    rLocalLHS(25,20)=0;
    rLocalLHS(25,21)=0;
    rLocalLHS(25,22)=0;
    rLocalLHS(25,23)=0;
    rLocalLHS(25,24)=0;
    rLocalLHS(25,25)=0;
    rLocalLHS(25,26)=0;
    rLocalLHS(25,27)=0;
    rLocalLHS(25,28)=0;
    rLocalLHS(25,29)=0;
    rLocalLHS(26,0)=0;
    rLocalLHS(26,1)=0;
    rLocalLHS(26,2)=clhs1;
    rLocalLHS(26,3)=0;
    rLocalLHS(26,4)=0;
    rLocalLHS(26,5)=clhs4;
    rLocalLHS(26,6)=0;
    rLocalLHS(26,7)=0;
    rLocalLHS(26,8)=clhs7;
    rLocalLHS(26,9)=0;
    rLocalLHS(26,10)=0;
    rLocalLHS(26,11)=clhs10;
    rLocalLHS(26,12)=0;
    rLocalLHS(26,13)=0;
    rLocalLHS(26,14)=DOperator(1,0);
    rLocalLHS(26,15)=0;
    rLocalLHS(26,16)=0;
    rLocalLHS(26,17)=DOperator(1,1);
    rLocalLHS(26,18)=0;
    rLocalLHS(26,19)=0;
    rLocalLHS(26,20)=DOperator(1,2);
    rLocalLHS(26,21)=0;
    rLocalLHS(26,22)=0;
    rLocalLHS(26,23)=0;
    rLocalLHS(26,24)=0;
    rLocalLHS(26,25)=0;
    rLocalLHS(26,26)=0;
    rLocalLHS(26,27)=0;
    rLocalLHS(26,28)=0;
    rLocalLHS(26,29)=0;
    rLocalLHS(27,0)=clhs2;
    rLocalLHS(27,1)=0;
    rLocalLHS(27,2)=0;
    rLocalLHS(27,3)=clhs5;
    rLocalLHS(27,4)=0;
    rLocalLHS(27,5)=0;
    rLocalLHS(27,6)=clhs8;
    rLocalLHS(27,7)=0;
    rLocalLHS(27,8)=0;
    rLocalLHS(27,9)=clhs11;
    rLocalLHS(27,10)=0;
    rLocalLHS(27,11)=0;
    rLocalLHS(27,12)=DOperator(2,0);
    rLocalLHS(27,13)=0;
    rLocalLHS(27,14)=0;
    rLocalLHS(27,15)=DOperator(2,1);
    rLocalLHS(27,16)=0;
    rLocalLHS(27,17)=0;
    rLocalLHS(27,18)=DOperator(2,2);
    rLocalLHS(27,19)=0;
    rLocalLHS(27,20)=0;
    rLocalLHS(27,21)=0;
    rLocalLHS(27,22)=0;
    rLocalLHS(27,23)=0;
    rLocalLHS(27,24)=0;
    rLocalLHS(27,25)=0;
    rLocalLHS(27,26)=0;
    rLocalLHS(27,27)=0;
    rLocalLHS(27,28)=0;
    rLocalLHS(27,29)=0;
    rLocalLHS(28,0)=0;
    rLocalLHS(28,1)=clhs2;
    rLocalLHS(28,2)=0;
    rLocalLHS(28,3)=0;
    rLocalLHS(28,4)=clhs5;
    rLocalLHS(28,5)=0;
    rLocalLHS(28,6)=0;
    rLocalLHS(28,7)=clhs8;
    rLocalLHS(28,8)=0;
    rLocalLHS(28,9)=0;
    rLocalLHS(28,10)=clhs11;
    rLocalLHS(28,11)=0;
    rLocalLHS(28,12)=0;
    rLocalLHS(28,13)=DOperator(2,0);
    rLocalLHS(28,14)=0;
    rLocalLHS(28,15)=0;
    rLocalLHS(28,16)=DOperator(2,1);
    rLocalLHS(28,17)=0;
    rLocalLHS(28,18)=0;
    rLocalLHS(28,19)=DOperator(2,2);
    rLocalLHS(28,20)=0;
    rLocalLHS(28,21)=0;
    rLocalLHS(28,22)=0;
    rLocalLHS(28,23)=0;
    rLocalLHS(28,24)=0;
    rLocalLHS(28,25)=0;
    rLocalLHS(28,26)=0;
    rLocalLHS(28,27)=0;
    rLocalLHS(28,28)=0;
    rLocalLHS(28,29)=0;
    rLocalLHS(29,0)=0;
    rLocalLHS(29,1)=0;
    rLocalLHS(29,2)=clhs2;
    rLocalLHS(29,3)=0;
    rLocalLHS(29,4)=0;
    rLocalLHS(29,5)=clhs5;
    rLocalLHS(29,6)=0;
    rLocalLHS(29,7)=0;
    rLocalLHS(29,8)=clhs8;
    rLocalLHS(29,9)=0;
    rLocalLHS(29,10)=0;
    rLocalLHS(29,11)=clhs11;
    rLocalLHS(29,12)=0;
    rLocalLHS(29,13)=0;
    rLocalLHS(29,14)=DOperator(2,0);
    rLocalLHS(29,15)=0;
    rLocalLHS(29,16)=0;
    rLocalLHS(29,17)=DOperator(2,1);
    rLocalLHS(29,18)=0;
    rLocalLHS(29,19)=0;
    rLocalLHS(29,20)=DOperator(2,2);
    rLocalLHS(29,21)=0;
    rLocalLHS(29,22)=0;
    rLocalLHS(29,23)=0;
    rLocalLHS(29,24)=0;
    rLocalLHS(29,25)=0;
    rLocalLHS(29,26)=0;
    rLocalLHS(29,27)=0;
    rLocalLHS(29,28)=0;
    rLocalLHS(29,29)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,8,4>::CalculateLocalLHS<MeshTyingMortarCondition<3,8,4>::ScalarValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 4, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(3,0);
    const double clhs4 =     -MOperator(0,1);
    const double clhs5 =     -MOperator(1,1);
    const double clhs6 =     -MOperator(2,1);
    const double clhs7 =     -MOperator(3,1);
    const double clhs8 =     -MOperator(0,2);
    const double clhs9 =     -MOperator(1,2);
    const double clhs10 =     -MOperator(2,2);
    const double clhs11 =     -MOperator(3,2);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=clhs0;
    rLocalLHS(0,8)=clhs1;
    rLocalLHS(0,9)=clhs2;
    rLocalLHS(0,10)=clhs3;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=clhs4;
    rLocalLHS(1,8)=clhs5;
    rLocalLHS(1,9)=clhs6;
    rLocalLHS(1,10)=clhs7;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=clhs8;
    rLocalLHS(2,8)=clhs9;
    rLocalLHS(2,9)=clhs10;
    rLocalLHS(2,10)=clhs11;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=DOperator(0,0);
    rLocalLHS(3,8)=DOperator(1,0);
    rLocalLHS(3,9)=DOperator(2,0);
    rLocalLHS(3,10)=DOperator(3,0);
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=DOperator(0,1);
    rLocalLHS(4,8)=DOperator(1,1);
    rLocalLHS(4,9)=DOperator(2,1);
    rLocalLHS(4,10)=DOperator(3,1);
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=DOperator(0,2);
    rLocalLHS(5,8)=DOperator(1,2);
    rLocalLHS(5,9)=DOperator(2,2);
    rLocalLHS(5,10)=DOperator(3,2);
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=DOperator(0,3);
    rLocalLHS(6,8)=DOperator(1,3);
    rLocalLHS(6,9)=DOperator(2,3);
    rLocalLHS(6,10)=DOperator(3,3);
    rLocalLHS(7,0)=clhs0;
    rLocalLHS(7,1)=clhs4;
    rLocalLHS(7,2)=clhs8;
    rLocalLHS(7,3)=DOperator(0,0);
    rLocalLHS(7,4)=DOperator(0,1);
    rLocalLHS(7,5)=DOperator(0,2);
    rLocalLHS(7,6)=DOperator(0,3);
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=0;
    rLocalLHS(7,10)=0;
    rLocalLHS(8,0)=clhs1;
    rLocalLHS(8,1)=clhs5;
    rLocalLHS(8,2)=clhs9;
    rLocalLHS(8,3)=DOperator(1,0);
    rLocalLHS(8,4)=DOperator(1,1);
    rLocalLHS(8,5)=DOperator(1,2);
    rLocalLHS(8,6)=DOperator(1,3);
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(9,0)=clhs2;
    rLocalLHS(9,1)=clhs6;
    rLocalLHS(9,2)=clhs10;
    rLocalLHS(9,3)=DOperator(2,0);
    rLocalLHS(9,4)=DOperator(2,1);
    rLocalLHS(9,5)=DOperator(2,2);
    rLocalLHS(9,6)=DOperator(2,3);
    rLocalLHS(9,7)=0;
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(10,0)=clhs3;
    rLocalLHS(10,1)=clhs7;
    rLocalLHS(10,2)=clhs11;
    rLocalLHS(10,3)=DOperator(3,0);
    rLocalLHS(10,4)=DOperator(3,1);
    rLocalLHS(10,5)=DOperator(3,2);
    rLocalLHS(10,6)=DOperator(3,3);
    rLocalLHS(10,7)=0;
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;

}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MeshTyingMortarCondition<3,8,4>::CalculateLocalLHS<MeshTyingMortarCondition<3,8,4>::Vector3DValue>(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // We get the mortar operators
    const BoundedMatrix<double, 4, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(3,0);
    const double clhs4 =     -MOperator(0,1);
    const double clhs5 =     -MOperator(1,1);
    const double clhs6 =     -MOperator(2,1);
    const double clhs7 =     -MOperator(3,1);
    const double clhs8 =     -MOperator(0,2);
    const double clhs9 =     -MOperator(1,2);
    const double clhs10 =     -MOperator(2,2);
    const double clhs11 =     -MOperator(3,2);

    rLocalLHS(0,0)=0;
    rLocalLHS(0,1)=0;
    rLocalLHS(0,2)=0;
    rLocalLHS(0,3)=0;
    rLocalLHS(0,4)=0;
    rLocalLHS(0,5)=0;
    rLocalLHS(0,6)=0;
    rLocalLHS(0,7)=0;
    rLocalLHS(0,8)=0;
    rLocalLHS(0,9)=0;
    rLocalLHS(0,10)=0;
    rLocalLHS(0,11)=0;
    rLocalLHS(0,12)=0;
    rLocalLHS(0,13)=0;
    rLocalLHS(0,14)=0;
    rLocalLHS(0,15)=0;
    rLocalLHS(0,16)=0;
    rLocalLHS(0,17)=0;
    rLocalLHS(0,18)=0;
    rLocalLHS(0,19)=0;
    rLocalLHS(0,20)=0;
    rLocalLHS(0,21)=clhs0;
    rLocalLHS(0,22)=0;
    rLocalLHS(0,23)=0;
    rLocalLHS(0,24)=clhs1;
    rLocalLHS(0,25)=0;
    rLocalLHS(0,26)=0;
    rLocalLHS(0,27)=clhs2;
    rLocalLHS(0,28)=0;
    rLocalLHS(0,29)=0;
    rLocalLHS(0,30)=clhs3;
    rLocalLHS(0,31)=0;
    rLocalLHS(0,32)=0;
    rLocalLHS(1,0)=0;
    rLocalLHS(1,1)=0;
    rLocalLHS(1,2)=0;
    rLocalLHS(1,3)=0;
    rLocalLHS(1,4)=0;
    rLocalLHS(1,5)=0;
    rLocalLHS(1,6)=0;
    rLocalLHS(1,7)=0;
    rLocalLHS(1,8)=0;
    rLocalLHS(1,9)=0;
    rLocalLHS(1,10)=0;
    rLocalLHS(1,11)=0;
    rLocalLHS(1,12)=0;
    rLocalLHS(1,13)=0;
    rLocalLHS(1,14)=0;
    rLocalLHS(1,15)=0;
    rLocalLHS(1,16)=0;
    rLocalLHS(1,17)=0;
    rLocalLHS(1,18)=0;
    rLocalLHS(1,19)=0;
    rLocalLHS(1,20)=0;
    rLocalLHS(1,21)=0;
    rLocalLHS(1,22)=clhs0;
    rLocalLHS(1,23)=0;
    rLocalLHS(1,24)=0;
    rLocalLHS(1,25)=clhs1;
    rLocalLHS(1,26)=0;
    rLocalLHS(1,27)=0;
    rLocalLHS(1,28)=clhs2;
    rLocalLHS(1,29)=0;
    rLocalLHS(1,30)=0;
    rLocalLHS(1,31)=clhs3;
    rLocalLHS(1,32)=0;
    rLocalLHS(2,0)=0;
    rLocalLHS(2,1)=0;
    rLocalLHS(2,2)=0;
    rLocalLHS(2,3)=0;
    rLocalLHS(2,4)=0;
    rLocalLHS(2,5)=0;
    rLocalLHS(2,6)=0;
    rLocalLHS(2,7)=0;
    rLocalLHS(2,8)=0;
    rLocalLHS(2,9)=0;
    rLocalLHS(2,10)=0;
    rLocalLHS(2,11)=0;
    rLocalLHS(2,12)=0;
    rLocalLHS(2,13)=0;
    rLocalLHS(2,14)=0;
    rLocalLHS(2,15)=0;
    rLocalLHS(2,16)=0;
    rLocalLHS(2,17)=0;
    rLocalLHS(2,18)=0;
    rLocalLHS(2,19)=0;
    rLocalLHS(2,20)=0;
    rLocalLHS(2,21)=0;
    rLocalLHS(2,22)=0;
    rLocalLHS(2,23)=clhs0;
    rLocalLHS(2,24)=0;
    rLocalLHS(2,25)=0;
    rLocalLHS(2,26)=clhs1;
    rLocalLHS(2,27)=0;
    rLocalLHS(2,28)=0;
    rLocalLHS(2,29)=clhs2;
    rLocalLHS(2,30)=0;
    rLocalLHS(2,31)=0;
    rLocalLHS(2,32)=clhs3;
    rLocalLHS(3,0)=0;
    rLocalLHS(3,1)=0;
    rLocalLHS(3,2)=0;
    rLocalLHS(3,3)=0;
    rLocalLHS(3,4)=0;
    rLocalLHS(3,5)=0;
    rLocalLHS(3,6)=0;
    rLocalLHS(3,7)=0;
    rLocalLHS(3,8)=0;
    rLocalLHS(3,9)=0;
    rLocalLHS(3,10)=0;
    rLocalLHS(3,11)=0;
    rLocalLHS(3,12)=0;
    rLocalLHS(3,13)=0;
    rLocalLHS(3,14)=0;
    rLocalLHS(3,15)=0;
    rLocalLHS(3,16)=0;
    rLocalLHS(3,17)=0;
    rLocalLHS(3,18)=0;
    rLocalLHS(3,19)=0;
    rLocalLHS(3,20)=0;
    rLocalLHS(3,21)=clhs4;
    rLocalLHS(3,22)=0;
    rLocalLHS(3,23)=0;
    rLocalLHS(3,24)=clhs5;
    rLocalLHS(3,25)=0;
    rLocalLHS(3,26)=0;
    rLocalLHS(3,27)=clhs6;
    rLocalLHS(3,28)=0;
    rLocalLHS(3,29)=0;
    rLocalLHS(3,30)=clhs7;
    rLocalLHS(3,31)=0;
    rLocalLHS(3,32)=0;
    rLocalLHS(4,0)=0;
    rLocalLHS(4,1)=0;
    rLocalLHS(4,2)=0;
    rLocalLHS(4,3)=0;
    rLocalLHS(4,4)=0;
    rLocalLHS(4,5)=0;
    rLocalLHS(4,6)=0;
    rLocalLHS(4,7)=0;
    rLocalLHS(4,8)=0;
    rLocalLHS(4,9)=0;
    rLocalLHS(4,10)=0;
    rLocalLHS(4,11)=0;
    rLocalLHS(4,12)=0;
    rLocalLHS(4,13)=0;
    rLocalLHS(4,14)=0;
    rLocalLHS(4,15)=0;
    rLocalLHS(4,16)=0;
    rLocalLHS(4,17)=0;
    rLocalLHS(4,18)=0;
    rLocalLHS(4,19)=0;
    rLocalLHS(4,20)=0;
    rLocalLHS(4,21)=0;
    rLocalLHS(4,22)=clhs4;
    rLocalLHS(4,23)=0;
    rLocalLHS(4,24)=0;
    rLocalLHS(4,25)=clhs5;
    rLocalLHS(4,26)=0;
    rLocalLHS(4,27)=0;
    rLocalLHS(4,28)=clhs6;
    rLocalLHS(4,29)=0;
    rLocalLHS(4,30)=0;
    rLocalLHS(4,31)=clhs7;
    rLocalLHS(4,32)=0;
    rLocalLHS(5,0)=0;
    rLocalLHS(5,1)=0;
    rLocalLHS(5,2)=0;
    rLocalLHS(5,3)=0;
    rLocalLHS(5,4)=0;
    rLocalLHS(5,5)=0;
    rLocalLHS(5,6)=0;
    rLocalLHS(5,7)=0;
    rLocalLHS(5,8)=0;
    rLocalLHS(5,9)=0;
    rLocalLHS(5,10)=0;
    rLocalLHS(5,11)=0;
    rLocalLHS(5,12)=0;
    rLocalLHS(5,13)=0;
    rLocalLHS(5,14)=0;
    rLocalLHS(5,15)=0;
    rLocalLHS(5,16)=0;
    rLocalLHS(5,17)=0;
    rLocalLHS(5,18)=0;
    rLocalLHS(5,19)=0;
    rLocalLHS(5,20)=0;
    rLocalLHS(5,21)=0;
    rLocalLHS(5,22)=0;
    rLocalLHS(5,23)=clhs4;
    rLocalLHS(5,24)=0;
    rLocalLHS(5,25)=0;
    rLocalLHS(5,26)=clhs5;
    rLocalLHS(5,27)=0;
    rLocalLHS(5,28)=0;
    rLocalLHS(5,29)=clhs6;
    rLocalLHS(5,30)=0;
    rLocalLHS(5,31)=0;
    rLocalLHS(5,32)=clhs7;
    rLocalLHS(6,0)=0;
    rLocalLHS(6,1)=0;
    rLocalLHS(6,2)=0;
    rLocalLHS(6,3)=0;
    rLocalLHS(6,4)=0;
    rLocalLHS(6,5)=0;
    rLocalLHS(6,6)=0;
    rLocalLHS(6,7)=0;
    rLocalLHS(6,8)=0;
    rLocalLHS(6,9)=0;
    rLocalLHS(6,10)=0;
    rLocalLHS(6,11)=0;
    rLocalLHS(6,12)=0;
    rLocalLHS(6,13)=0;
    rLocalLHS(6,14)=0;
    rLocalLHS(6,15)=0;
    rLocalLHS(6,16)=0;
    rLocalLHS(6,17)=0;
    rLocalLHS(6,18)=0;
    rLocalLHS(6,19)=0;
    rLocalLHS(6,20)=0;
    rLocalLHS(6,21)=clhs8;
    rLocalLHS(6,22)=0;
    rLocalLHS(6,23)=0;
    rLocalLHS(6,24)=clhs9;
    rLocalLHS(6,25)=0;
    rLocalLHS(6,26)=0;
    rLocalLHS(6,27)=clhs10;
    rLocalLHS(6,28)=0;
    rLocalLHS(6,29)=0;
    rLocalLHS(6,30)=clhs11;
    rLocalLHS(6,31)=0;
    rLocalLHS(6,32)=0;
    rLocalLHS(7,0)=0;
    rLocalLHS(7,1)=0;
    rLocalLHS(7,2)=0;
    rLocalLHS(7,3)=0;
    rLocalLHS(7,4)=0;
    rLocalLHS(7,5)=0;
    rLocalLHS(7,6)=0;
    rLocalLHS(7,7)=0;
    rLocalLHS(7,8)=0;
    rLocalLHS(7,9)=0;
    rLocalLHS(7,10)=0;
    rLocalLHS(7,11)=0;
    rLocalLHS(7,12)=0;
    rLocalLHS(7,13)=0;
    rLocalLHS(7,14)=0;
    rLocalLHS(7,15)=0;
    rLocalLHS(7,16)=0;
    rLocalLHS(7,17)=0;
    rLocalLHS(7,18)=0;
    rLocalLHS(7,19)=0;
    rLocalLHS(7,20)=0;
    rLocalLHS(7,21)=0;
    rLocalLHS(7,22)=clhs8;
    rLocalLHS(7,23)=0;
    rLocalLHS(7,24)=0;
    rLocalLHS(7,25)=clhs9;
    rLocalLHS(7,26)=0;
    rLocalLHS(7,27)=0;
    rLocalLHS(7,28)=clhs10;
    rLocalLHS(7,29)=0;
    rLocalLHS(7,30)=0;
    rLocalLHS(7,31)=clhs11;
    rLocalLHS(7,32)=0;
    rLocalLHS(8,0)=0;
    rLocalLHS(8,1)=0;
    rLocalLHS(8,2)=0;
    rLocalLHS(8,3)=0;
    rLocalLHS(8,4)=0;
    rLocalLHS(8,5)=0;
    rLocalLHS(8,6)=0;
    rLocalLHS(8,7)=0;
    rLocalLHS(8,8)=0;
    rLocalLHS(8,9)=0;
    rLocalLHS(8,10)=0;
    rLocalLHS(8,11)=0;
    rLocalLHS(8,12)=0;
    rLocalLHS(8,13)=0;
    rLocalLHS(8,14)=0;
    rLocalLHS(8,15)=0;
    rLocalLHS(8,16)=0;
    rLocalLHS(8,17)=0;
    rLocalLHS(8,18)=0;
    rLocalLHS(8,19)=0;
    rLocalLHS(8,20)=0;
    rLocalLHS(8,21)=0;
    rLocalLHS(8,22)=0;
    rLocalLHS(8,23)=clhs8;
    rLocalLHS(8,24)=0;
    rLocalLHS(8,25)=0;
    rLocalLHS(8,26)=clhs9;
    rLocalLHS(8,27)=0;
    rLocalLHS(8,28)=0;
    rLocalLHS(8,29)=clhs10;
    rLocalLHS(8,30)=0;
    rLocalLHS(8,31)=0;
    rLocalLHS(8,32)=clhs11;
    rLocalLHS(9,0)=0;
    rLocalLHS(9,1)=0;
    rLocalLHS(9,2)=0;
    rLocalLHS(9,3)=0;
    rLocalLHS(9,4)=0;
    rLocalLHS(9,5)=0;
    rLocalLHS(9,6)=0;
    rLocalLHS(9,7)=0;
    rLocalLHS(9,8)=0;
    rLocalLHS(9,9)=0;
    rLocalLHS(9,10)=0;
    rLocalLHS(9,11)=0;
    rLocalLHS(9,12)=0;
    rLocalLHS(9,13)=0;
    rLocalLHS(9,14)=0;
    rLocalLHS(9,15)=0;
    rLocalLHS(9,16)=0;
    rLocalLHS(9,17)=0;
    rLocalLHS(9,18)=0;
    rLocalLHS(9,19)=0;
    rLocalLHS(9,20)=0;
    rLocalLHS(9,21)=DOperator(0,0);
    rLocalLHS(9,22)=0;
    rLocalLHS(9,23)=0;
    rLocalLHS(9,24)=DOperator(1,0);
    rLocalLHS(9,25)=0;
    rLocalLHS(9,26)=0;
    rLocalLHS(9,27)=DOperator(2,0);
    rLocalLHS(9,28)=0;
    rLocalLHS(9,29)=0;
    rLocalLHS(9,30)=DOperator(3,0);
    rLocalLHS(9,31)=0;
    rLocalLHS(9,32)=0;
    rLocalLHS(10,0)=0;
    rLocalLHS(10,1)=0;
    rLocalLHS(10,2)=0;
    rLocalLHS(10,3)=0;
    rLocalLHS(10,4)=0;
    rLocalLHS(10,5)=0;
    rLocalLHS(10,6)=0;
    rLocalLHS(10,7)=0;
    rLocalLHS(10,8)=0;
    rLocalLHS(10,9)=0;
    rLocalLHS(10,10)=0;
    rLocalLHS(10,11)=0;
    rLocalLHS(10,12)=0;
    rLocalLHS(10,13)=0;
    rLocalLHS(10,14)=0;
    rLocalLHS(10,15)=0;
    rLocalLHS(10,16)=0;
    rLocalLHS(10,17)=0;
    rLocalLHS(10,18)=0;
    rLocalLHS(10,19)=0;
    rLocalLHS(10,20)=0;
    rLocalLHS(10,21)=0;
    rLocalLHS(10,22)=DOperator(0,0);
    rLocalLHS(10,23)=0;
    rLocalLHS(10,24)=0;
    rLocalLHS(10,25)=DOperator(1,0);
    rLocalLHS(10,26)=0;
    rLocalLHS(10,27)=0;
    rLocalLHS(10,28)=DOperator(2,0);
    rLocalLHS(10,29)=0;
    rLocalLHS(10,30)=0;
    rLocalLHS(10,31)=DOperator(3,0);
    rLocalLHS(10,32)=0;
    rLocalLHS(11,0)=0;
    rLocalLHS(11,1)=0;
    rLocalLHS(11,2)=0;
    rLocalLHS(11,3)=0;
    rLocalLHS(11,4)=0;
    rLocalLHS(11,5)=0;
    rLocalLHS(11,6)=0;
    rLocalLHS(11,7)=0;
    rLocalLHS(11,8)=0;
    rLocalLHS(11,9)=0;
    rLocalLHS(11,10)=0;
    rLocalLHS(11,11)=0;
    rLocalLHS(11,12)=0;
    rLocalLHS(11,13)=0;
    rLocalLHS(11,14)=0;
    rLocalLHS(11,15)=0;
    rLocalLHS(11,16)=0;
    rLocalLHS(11,17)=0;
    rLocalLHS(11,18)=0;
    rLocalLHS(11,19)=0;
    rLocalLHS(11,20)=0;
    rLocalLHS(11,21)=0;
    rLocalLHS(11,22)=0;
    rLocalLHS(11,23)=DOperator(0,0);
    rLocalLHS(11,24)=0;
    rLocalLHS(11,25)=0;
    rLocalLHS(11,26)=DOperator(1,0);
    rLocalLHS(11,27)=0;
    rLocalLHS(11,28)=0;
    rLocalLHS(11,29)=DOperator(2,0);
    rLocalLHS(11,30)=0;
    rLocalLHS(11,31)=0;
    rLocalLHS(11,32)=DOperator(3,0);
    rLocalLHS(12,0)=0;
    rLocalLHS(12,1)=0;
    rLocalLHS(12,2)=0;
    rLocalLHS(12,3)=0;
    rLocalLHS(12,4)=0;
    rLocalLHS(12,5)=0;
    rLocalLHS(12,6)=0;
    rLocalLHS(12,7)=0;
    rLocalLHS(12,8)=0;
    rLocalLHS(12,9)=0;
    rLocalLHS(12,10)=0;
    rLocalLHS(12,11)=0;
    rLocalLHS(12,12)=0;
    rLocalLHS(12,13)=0;
    rLocalLHS(12,14)=0;
    rLocalLHS(12,15)=0;
    rLocalLHS(12,16)=0;
    rLocalLHS(12,17)=0;
    rLocalLHS(12,18)=0;
    rLocalLHS(12,19)=0;
    rLocalLHS(12,20)=0;
    rLocalLHS(12,21)=DOperator(0,1);
    rLocalLHS(12,22)=0;
    rLocalLHS(12,23)=0;
    rLocalLHS(12,24)=DOperator(1,1);
    rLocalLHS(12,25)=0;
    rLocalLHS(12,26)=0;
    rLocalLHS(12,27)=DOperator(2,1);
    rLocalLHS(12,28)=0;
    rLocalLHS(12,29)=0;
    rLocalLHS(12,30)=DOperator(3,1);
    rLocalLHS(12,31)=0;
    rLocalLHS(12,32)=0;
    rLocalLHS(13,0)=0;
    rLocalLHS(13,1)=0;
    rLocalLHS(13,2)=0;
    rLocalLHS(13,3)=0;
    rLocalLHS(13,4)=0;
    rLocalLHS(13,5)=0;
    rLocalLHS(13,6)=0;
    rLocalLHS(13,7)=0;
    rLocalLHS(13,8)=0;
    rLocalLHS(13,9)=0;
    rLocalLHS(13,10)=0;
    rLocalLHS(13,11)=0;
    rLocalLHS(13,12)=0;
    rLocalLHS(13,13)=0;
    rLocalLHS(13,14)=0;
    rLocalLHS(13,15)=0;
    rLocalLHS(13,16)=0;
    rLocalLHS(13,17)=0;
    rLocalLHS(13,18)=0;
    rLocalLHS(13,19)=0;
    rLocalLHS(13,20)=0;
    rLocalLHS(13,21)=0;
    rLocalLHS(13,22)=DOperator(0,1);
    rLocalLHS(13,23)=0;
    rLocalLHS(13,24)=0;
    rLocalLHS(13,25)=DOperator(1,1);
    rLocalLHS(13,26)=0;
    rLocalLHS(13,27)=0;
    rLocalLHS(13,28)=DOperator(2,1);
    rLocalLHS(13,29)=0;
    rLocalLHS(13,30)=0;
    rLocalLHS(13,31)=DOperator(3,1);
    rLocalLHS(13,32)=0;
    rLocalLHS(14,0)=0;
    rLocalLHS(14,1)=0;
    rLocalLHS(14,2)=0;
    rLocalLHS(14,3)=0;
    rLocalLHS(14,4)=0;
    rLocalLHS(14,5)=0;
    rLocalLHS(14,6)=0;
    rLocalLHS(14,7)=0;
    rLocalLHS(14,8)=0;
    rLocalLHS(14,9)=0;
    rLocalLHS(14,10)=0;
    rLocalLHS(14,11)=0;
    rLocalLHS(14,12)=0;
    rLocalLHS(14,13)=0;
    rLocalLHS(14,14)=0;
    rLocalLHS(14,15)=0;
    rLocalLHS(14,16)=0;
    rLocalLHS(14,17)=0;
    rLocalLHS(14,18)=0;
    rLocalLHS(14,19)=0;
    rLocalLHS(14,20)=0;
    rLocalLHS(14,21)=0;
    rLocalLHS(14,22)=0;
    rLocalLHS(14,23)=DOperator(0,1);
    rLocalLHS(14,24)=0;
    rLocalLHS(14,25)=0;
    rLocalLHS(14,26)=DOperator(1,1);
    rLocalLHS(14,27)=0;
    rLocalLHS(14,28)=0;
    rLocalLHS(14,29)=DOperator(2,1);
    rLocalLHS(14,30)=0;
    rLocalLHS(14,31)=0;
    rLocalLHS(14,32)=DOperator(3,1);
    rLocalLHS(15,0)=0;
    rLocalLHS(15,1)=0;
    rLocalLHS(15,2)=0;
    rLocalLHS(15,3)=0;
    rLocalLHS(15,4)=0;
    rLocalLHS(15,5)=0;
    rLocalLHS(15,6)=0;
    rLocalLHS(15,7)=0;
    rLocalLHS(15,8)=0;
    rLocalLHS(15,9)=0;
    rLocalLHS(15,10)=0;
    rLocalLHS(15,11)=0;
    rLocalLHS(15,12)=0;
    rLocalLHS(15,13)=0;
    rLocalLHS(15,14)=0;
    rLocalLHS(15,15)=0;
    rLocalLHS(15,16)=0;
    rLocalLHS(15,17)=0;
    rLocalLHS(15,18)=0;
    rLocalLHS(15,19)=0;
    rLocalLHS(15,20)=0;
    rLocalLHS(15,21)=DOperator(0,2);
    rLocalLHS(15,22)=0;
    rLocalLHS(15,23)=0;
    rLocalLHS(15,24)=DOperator(1,2);
    rLocalLHS(15,25)=0;
    rLocalLHS(15,26)=0;
    rLocalLHS(15,27)=DOperator(2,2);
    rLocalLHS(15,28)=0;
    rLocalLHS(15,29)=0;
    rLocalLHS(15,30)=DOperator(3,2);
    rLocalLHS(15,31)=0;
    rLocalLHS(15,32)=0;
    rLocalLHS(16,0)=0;
    rLocalLHS(16,1)=0;
    rLocalLHS(16,2)=0;
    rLocalLHS(16,3)=0;
    rLocalLHS(16,4)=0;
    rLocalLHS(16,5)=0;
    rLocalLHS(16,6)=0;
    rLocalLHS(16,7)=0;
    rLocalLHS(16,8)=0;
    rLocalLHS(16,9)=0;
    rLocalLHS(16,10)=0;
    rLocalLHS(16,11)=0;
    rLocalLHS(16,12)=0;
    rLocalLHS(16,13)=0;
    rLocalLHS(16,14)=0;
    rLocalLHS(16,15)=0;
    rLocalLHS(16,16)=0;
    rLocalLHS(16,17)=0;
    rLocalLHS(16,18)=0;
    rLocalLHS(16,19)=0;
    rLocalLHS(16,20)=0;
    rLocalLHS(16,21)=0;
    rLocalLHS(16,22)=DOperator(0,2);
    rLocalLHS(16,23)=0;
    rLocalLHS(16,24)=0;
    rLocalLHS(16,25)=DOperator(1,2);
    rLocalLHS(16,26)=0;
    rLocalLHS(16,27)=0;
    rLocalLHS(16,28)=DOperator(2,2);
    rLocalLHS(16,29)=0;
    rLocalLHS(16,30)=0;
    rLocalLHS(16,31)=DOperator(3,2);
    rLocalLHS(16,32)=0;
    rLocalLHS(17,0)=0;
    rLocalLHS(17,1)=0;
    rLocalLHS(17,2)=0;
    rLocalLHS(17,3)=0;
    rLocalLHS(17,4)=0;
    rLocalLHS(17,5)=0;
    rLocalLHS(17,6)=0;
    rLocalLHS(17,7)=0;
    rLocalLHS(17,8)=0;
    rLocalLHS(17,9)=0;
    rLocalLHS(17,10)=0;
    rLocalLHS(17,11)=0;
    rLocalLHS(17,12)=0;
    rLocalLHS(17,13)=0;
    rLocalLHS(17,14)=0;
    rLocalLHS(17,15)=0;
    rLocalLHS(17,16)=0;
    rLocalLHS(17,17)=0;
    rLocalLHS(17,18)=0;
    rLocalLHS(17,19)=0;
    rLocalLHS(17,20)=0;
    rLocalLHS(17,21)=0;
    rLocalLHS(17,22)=0;
    rLocalLHS(17,23)=DOperator(0,2);
    rLocalLHS(17,24)=0;
    rLocalLHS(17,25)=0;
    rLocalLHS(17,26)=DOperator(1,2);
    rLocalLHS(17,27)=0;
    rLocalLHS(17,28)=0;
    rLocalLHS(17,29)=DOperator(2,2);
    rLocalLHS(17,30)=0;
    rLocalLHS(17,31)=0;
    rLocalLHS(17,32)=DOperator(3,2);
    rLocalLHS(18,0)=0;
    rLocalLHS(18,1)=0;
    rLocalLHS(18,2)=0;
    rLocalLHS(18,3)=0;
    rLocalLHS(18,4)=0;
    rLocalLHS(18,5)=0;
    rLocalLHS(18,6)=0;
    rLocalLHS(18,7)=0;
    rLocalLHS(18,8)=0;
    rLocalLHS(18,9)=0;
    rLocalLHS(18,10)=0;
    rLocalLHS(18,11)=0;
    rLocalLHS(18,12)=0;
    rLocalLHS(18,13)=0;
    rLocalLHS(18,14)=0;
    rLocalLHS(18,15)=0;
    rLocalLHS(18,16)=0;
    rLocalLHS(18,17)=0;
    rLocalLHS(18,18)=0;
    rLocalLHS(18,19)=0;
    rLocalLHS(18,20)=0;
    rLocalLHS(18,21)=DOperator(0,3);
    rLocalLHS(18,22)=0;
    rLocalLHS(18,23)=0;
    rLocalLHS(18,24)=DOperator(1,3);
    rLocalLHS(18,25)=0;
    rLocalLHS(18,26)=0;
    rLocalLHS(18,27)=DOperator(2,3);
    rLocalLHS(18,28)=0;
    rLocalLHS(18,29)=0;
    rLocalLHS(18,30)=DOperator(3,3);
    rLocalLHS(18,31)=0;
    rLocalLHS(18,32)=0;
    rLocalLHS(19,0)=0;
    rLocalLHS(19,1)=0;
    rLocalLHS(19,2)=0;
    rLocalLHS(19,3)=0;
    rLocalLHS(19,4)=0;
    rLocalLHS(19,5)=0;
    rLocalLHS(19,6)=0;
    rLocalLHS(19,7)=0;
    rLocalLHS(19,8)=0;
    rLocalLHS(19,9)=0;
    rLocalLHS(19,10)=0;
    rLocalLHS(19,11)=0;
    rLocalLHS(19,12)=0;
    rLocalLHS(19,13)=0;
    rLocalLHS(19,14)=0;
    rLocalLHS(19,15)=0;
    rLocalLHS(19,16)=0;
    rLocalLHS(19,17)=0;
    rLocalLHS(19,18)=0;
    rLocalLHS(19,19)=0;
    rLocalLHS(19,20)=0;
    rLocalLHS(19,21)=0;
    rLocalLHS(19,22)=DOperator(0,3);
    rLocalLHS(19,23)=0;
    rLocalLHS(19,24)=0;
    rLocalLHS(19,25)=DOperator(1,3);
    rLocalLHS(19,26)=0;
    rLocalLHS(19,27)=0;
    rLocalLHS(19,28)=DOperator(2,3);
    rLocalLHS(19,29)=0;
    rLocalLHS(19,30)=0;
    rLocalLHS(19,31)=DOperator(3,3);
    rLocalLHS(19,32)=0;
    rLocalLHS(20,0)=0;
    rLocalLHS(20,1)=0;
    rLocalLHS(20,2)=0;
    rLocalLHS(20,3)=0;
    rLocalLHS(20,4)=0;
    rLocalLHS(20,5)=0;
    rLocalLHS(20,6)=0;
    rLocalLHS(20,7)=0;
    rLocalLHS(20,8)=0;
    rLocalLHS(20,9)=0;
    rLocalLHS(20,10)=0;
    rLocalLHS(20,11)=0;
    rLocalLHS(20,12)=0;
    rLocalLHS(20,13)=0;
    rLocalLHS(20,14)=0;
    rLocalLHS(20,15)=0;
    rLocalLHS(20,16)=0;
    rLocalLHS(20,17)=0;
    rLocalLHS(20,18)=0;
    rLocalLHS(20,19)=0;
    rLocalLHS(20,20)=0;
    rLocalLHS(20,21)=0;
    rLocalLHS(20,22)=0;
    rLocalLHS(20,23)=DOperator(0,3);
    rLocalLHS(20,24)=0;
    rLocalLHS(20,25)=0;
    rLocalLHS(20,26)=DOperator(1,3);
    rLocalLHS(20,27)=0;
    rLocalLHS(20,28)=0;
    rLocalLHS(20,29)=DOperator(2,3);
    rLocalLHS(20,30)=0;
    rLocalLHS(20,31)=0;
    rLocalLHS(20,32)=DOperator(3,3);
    rLocalLHS(21,0)=clhs0;
    rLocalLHS(21,1)=0;
    rLocalLHS(21,2)=0;
    rLocalLHS(21,3)=clhs4;
    rLocalLHS(21,4)=0;
    rLocalLHS(21,5)=0;
    rLocalLHS(21,6)=clhs8;
    rLocalLHS(21,7)=0;
    rLocalLHS(21,8)=0;
    rLocalLHS(21,9)=DOperator(0,0);
    rLocalLHS(21,10)=0;
    rLocalLHS(21,11)=0;
    rLocalLHS(21,12)=DOperator(0,1);
    rLocalLHS(21,13)=0;
    rLocalLHS(21,14)=0;
    rLocalLHS(21,15)=DOperator(0,2);
    rLocalLHS(21,16)=0;
    rLocalLHS(21,17)=0;
    rLocalLHS(21,18)=DOperator(0,3);
    rLocalLHS(21,19)=0;
    rLocalLHS(21,20)=0;
    rLocalLHS(21,21)=0;
    rLocalLHS(21,22)=0;
    rLocalLHS(21,23)=0;
    rLocalLHS(21,24)=0;
    rLocalLHS(21,25)=0;
    rLocalLHS(21,26)=0;
    rLocalLHS(21,27)=0;
    rLocalLHS(21,28)=0;
    rLocalLHS(21,29)=0;
    rLocalLHS(21,30)=0;
    rLocalLHS(21,31)=0;
    rLocalLHS(21,32)=0;
    rLocalLHS(22,0)=0;
    rLocalLHS(22,1)=clhs0;
    rLocalLHS(22,2)=0;
    rLocalLHS(22,3)=0;
    rLocalLHS(22,4)=clhs4;
    rLocalLHS(22,5)=0;
    rLocalLHS(22,6)=0;
    rLocalLHS(22,7)=clhs8;
    rLocalLHS(22,8)=0;
    rLocalLHS(22,9)=0;
    rLocalLHS(22,10)=DOperator(0,0);
    rLocalLHS(22,11)=0;
    rLocalLHS(22,12)=0;
    rLocalLHS(22,13)=DOperator(0,1);
    rLocalLHS(22,14)=0;
    rLocalLHS(22,15)=0;
    rLocalLHS(22,16)=DOperator(0,2);
    rLocalLHS(22,17)=0;
    rLocalLHS(22,18)=0;
    rLocalLHS(22,19)=DOperator(0,3);
    rLocalLHS(22,20)=0;
    rLocalLHS(22,21)=0;
    rLocalLHS(22,22)=0;
    rLocalLHS(22,23)=0;
    rLocalLHS(22,24)=0;
    rLocalLHS(22,25)=0;
    rLocalLHS(22,26)=0;
    rLocalLHS(22,27)=0;
    rLocalLHS(22,28)=0;
    rLocalLHS(22,29)=0;
    rLocalLHS(22,30)=0;
    rLocalLHS(22,31)=0;
    rLocalLHS(22,32)=0;
    rLocalLHS(23,0)=0;
    rLocalLHS(23,1)=0;
    rLocalLHS(23,2)=clhs0;
    rLocalLHS(23,3)=0;
    rLocalLHS(23,4)=0;
    rLocalLHS(23,5)=clhs4;
    rLocalLHS(23,6)=0;
    rLocalLHS(23,7)=0;
    rLocalLHS(23,8)=clhs8;
    rLocalLHS(23,9)=0;
    rLocalLHS(23,10)=0;
    rLocalLHS(23,11)=DOperator(0,0);
    rLocalLHS(23,12)=0;
    rLocalLHS(23,13)=0;
    rLocalLHS(23,14)=DOperator(0,1);
    rLocalLHS(23,15)=0;
    rLocalLHS(23,16)=0;
    rLocalLHS(23,17)=DOperator(0,2);
    rLocalLHS(23,18)=0;
    rLocalLHS(23,19)=0;
    rLocalLHS(23,20)=DOperator(0,3);
    rLocalLHS(23,21)=0;
    rLocalLHS(23,22)=0;
    rLocalLHS(23,23)=0;
    rLocalLHS(23,24)=0;
    rLocalLHS(23,25)=0;
    rLocalLHS(23,26)=0;
    rLocalLHS(23,27)=0;
    rLocalLHS(23,28)=0;
    rLocalLHS(23,29)=0;
    rLocalLHS(23,30)=0;
    rLocalLHS(23,31)=0;
    rLocalLHS(23,32)=0;
    rLocalLHS(24,0)=clhs1;
    rLocalLHS(24,1)=0;
    rLocalLHS(24,2)=0;
    rLocalLHS(24,3)=clhs5;
    rLocalLHS(24,4)=0;
    rLocalLHS(24,5)=0;
    rLocalLHS(24,6)=clhs9;
    rLocalLHS(24,7)=0;
    rLocalLHS(24,8)=0;
    rLocalLHS(24,9)=DOperator(1,0);
    rLocalLHS(24,10)=0;
    rLocalLHS(24,11)=0;
    rLocalLHS(24,12)=DOperator(1,1);
    rLocalLHS(24,13)=0;
    rLocalLHS(24,14)=0;
    rLocalLHS(24,15)=DOperator(1,2);
    rLocalLHS(24,16)=0;
    rLocalLHS(24,17)=0;
    rLocalLHS(24,18)=DOperator(1,3);
    rLocalLHS(24,19)=0;
    rLocalLHS(24,20)=0;
    rLocalLHS(24,21)=0;
    rLocalLHS(24,22)=0;
    rLocalLHS(24,23)=0;
    rLocalLHS(24,24)=0;
    rLocalLHS(24,25)=0;
    rLocalLHS(24,26)=0;
    rLocalLHS(24,27)=0;
    rLocalLHS(24,28)=0;
    rLocalLHS(24,29)=0;
    rLocalLHS(24,30)=0;
    rLocalLHS(24,31)=0;
    rLocalLHS(24,32)=0;
    rLocalLHS(25,0)=0;
    rLocalLHS(25,1)=clhs1;
    rLocalLHS(25,2)=0;
    rLocalLHS(25,3)=0;
    rLocalLHS(25,4)=clhs5;
    rLocalLHS(25,5)=0;
    rLocalLHS(25,6)=0;
    rLocalLHS(25,7)=clhs9;
    rLocalLHS(25,8)=0;
    rLocalLHS(25,9)=0;
    rLocalLHS(25,10)=DOperator(1,0);
    rLocalLHS(25,11)=0;
    rLocalLHS(25,12)=0;
    rLocalLHS(25,13)=DOperator(1,1);
    rLocalLHS(25,14)=0;
    rLocalLHS(25,15)=0;
    rLocalLHS(25,16)=DOperator(1,2);
    rLocalLHS(25,17)=0;
    rLocalLHS(25,18)=0;
    rLocalLHS(25,19)=DOperator(1,3);
    rLocalLHS(25,20)=0;
    rLocalLHS(25,21)=0;
    rLocalLHS(25,22)=0;
    rLocalLHS(25,23)=0;
    rLocalLHS(25,24)=0;
    rLocalLHS(25,25)=0;
    rLocalLHS(25,26)=0;
    rLocalLHS(25,27)=0;
    rLocalLHS(25,28)=0;
    rLocalLHS(25,29)=0;
    rLocalLHS(25,30)=0;
    rLocalLHS(25,31)=0;
    rLocalLHS(25,32)=0;
    rLocalLHS(26,0)=0;
    rLocalLHS(26,1)=0;
    rLocalLHS(26,2)=clhs1;
    rLocalLHS(26,3)=0;
    rLocalLHS(26,4)=0;
    rLocalLHS(26,5)=clhs5;
    rLocalLHS(26,6)=0;
    rLocalLHS(26,7)=0;
    rLocalLHS(26,8)=clhs9;
    rLocalLHS(26,9)=0;
    rLocalLHS(26,10)=0;
    rLocalLHS(26,11)=DOperator(1,0);
    rLocalLHS(26,12)=0;
    rLocalLHS(26,13)=0;
    rLocalLHS(26,14)=DOperator(1,1);
    rLocalLHS(26,15)=0;
    rLocalLHS(26,16)=0;
    rLocalLHS(26,17)=DOperator(1,2);
    rLocalLHS(26,18)=0;
    rLocalLHS(26,19)=0;
    rLocalLHS(26,20)=DOperator(1,3);
    rLocalLHS(26,21)=0;
    rLocalLHS(26,22)=0;
    rLocalLHS(26,23)=0;
    rLocalLHS(26,24)=0;
    rLocalLHS(26,25)=0;
    rLocalLHS(26,26)=0;
    rLocalLHS(26,27)=0;
    rLocalLHS(26,28)=0;
    rLocalLHS(26,29)=0;
    rLocalLHS(26,30)=0;
    rLocalLHS(26,31)=0;
    rLocalLHS(26,32)=0;
    rLocalLHS(27,0)=clhs2;
    rLocalLHS(27,1)=0;
    rLocalLHS(27,2)=0;
    rLocalLHS(27,3)=clhs6;
    rLocalLHS(27,4)=0;
    rLocalLHS(27,5)=0;
    rLocalLHS(27,6)=clhs10;
    rLocalLHS(27,7)=0;
    rLocalLHS(27,8)=0;
    rLocalLHS(27,9)=DOperator(2,0);
    rLocalLHS(27,10)=0;
    rLocalLHS(27,11)=0;
    rLocalLHS(27,12)=DOperator(2,1);
    rLocalLHS(27,13)=0;
    rLocalLHS(27,14)=0;
    rLocalLHS(27,15)=DOperator(2,2);
    rLocalLHS(27,16)=0;
    rLocalLHS(27,17)=0;
    rLocalLHS(27,18)=DOperator(2,3);
    rLocalLHS(27,19)=0;
    rLocalLHS(27,20)=0;
    rLocalLHS(27,21)=0;
    rLocalLHS(27,22)=0;
    rLocalLHS(27,23)=0;
    rLocalLHS(27,24)=0;
    rLocalLHS(27,25)=0;
    rLocalLHS(27,26)=0;
    rLocalLHS(27,27)=0;
    rLocalLHS(27,28)=0;
    rLocalLHS(27,29)=0;
    rLocalLHS(27,30)=0;
    rLocalLHS(27,31)=0;
    rLocalLHS(27,32)=0;
    rLocalLHS(28,0)=0;
    rLocalLHS(28,1)=clhs2;
    rLocalLHS(28,2)=0;
    rLocalLHS(28,3)=0;
    rLocalLHS(28,4)=clhs6;
    rLocalLHS(28,5)=0;
    rLocalLHS(28,6)=0;
    rLocalLHS(28,7)=clhs10;
    rLocalLHS(28,8)=0;
    rLocalLHS(28,9)=0;
    rLocalLHS(28,10)=DOperator(2,0);
    rLocalLHS(28,11)=0;
    rLocalLHS(28,12)=0;
    rLocalLHS(28,13)=DOperator(2,1);
    rLocalLHS(28,14)=0;
    rLocalLHS(28,15)=0;
    rLocalLHS(28,16)=DOperator(2,2);
    rLocalLHS(28,17)=0;
    rLocalLHS(28,18)=0;
    rLocalLHS(28,19)=DOperator(2,3);
    rLocalLHS(28,20)=0;
    rLocalLHS(28,21)=0;
    rLocalLHS(28,22)=0;
    rLocalLHS(28,23)=0;
    rLocalLHS(28,24)=0;
    rLocalLHS(28,25)=0;
    rLocalLHS(28,26)=0;
    rLocalLHS(28,27)=0;
    rLocalLHS(28,28)=0;
    rLocalLHS(28,29)=0;
    rLocalLHS(28,30)=0;
    rLocalLHS(28,31)=0;
    rLocalLHS(28,32)=0;
    rLocalLHS(29,0)=0;
    rLocalLHS(29,1)=0;
    rLocalLHS(29,2)=clhs2;
    rLocalLHS(29,3)=0;
    rLocalLHS(29,4)=0;
    rLocalLHS(29,5)=clhs6;
    rLocalLHS(29,6)=0;
    rLocalLHS(29,7)=0;
    rLocalLHS(29,8)=clhs10;
    rLocalLHS(29,9)=0;
    rLocalLHS(29,10)=0;
    rLocalLHS(29,11)=DOperator(2,0);
    rLocalLHS(29,12)=0;
    rLocalLHS(29,13)=0;
    rLocalLHS(29,14)=DOperator(2,1);
    rLocalLHS(29,15)=0;
    rLocalLHS(29,16)=0;
    rLocalLHS(29,17)=DOperator(2,2);
    rLocalLHS(29,18)=0;
    rLocalLHS(29,19)=0;
    rLocalLHS(29,20)=DOperator(2,3);
    rLocalLHS(29,21)=0;
    rLocalLHS(29,22)=0;
    rLocalLHS(29,23)=0;
    rLocalLHS(29,24)=0;
    rLocalLHS(29,25)=0;
    rLocalLHS(29,26)=0;
    rLocalLHS(29,27)=0;
    rLocalLHS(29,28)=0;
    rLocalLHS(29,29)=0;
    rLocalLHS(29,30)=0;
    rLocalLHS(29,31)=0;
    rLocalLHS(29,32)=0;
    rLocalLHS(30,0)=clhs3;
    rLocalLHS(30,1)=0;
    rLocalLHS(30,2)=0;
    rLocalLHS(30,3)=clhs7;
    rLocalLHS(30,4)=0;
    rLocalLHS(30,5)=0;
    rLocalLHS(30,6)=clhs11;
    rLocalLHS(30,7)=0;
    rLocalLHS(30,8)=0;
    rLocalLHS(30,9)=DOperator(3,0);
    rLocalLHS(30,10)=0;
    rLocalLHS(30,11)=0;
    rLocalLHS(30,12)=DOperator(3,1);
    rLocalLHS(30,13)=0;
    rLocalLHS(30,14)=0;
    rLocalLHS(30,15)=DOperator(3,2);
    rLocalLHS(30,16)=0;
    rLocalLHS(30,17)=0;
    rLocalLHS(30,18)=DOperator(3,3);
    rLocalLHS(30,19)=0;
    rLocalLHS(30,20)=0;
    rLocalLHS(30,21)=0;
    rLocalLHS(30,22)=0;
    rLocalLHS(30,23)=0;
    rLocalLHS(30,24)=0;
    rLocalLHS(30,25)=0;
    rLocalLHS(30,26)=0;
    rLocalLHS(30,27)=0;
    rLocalLHS(30,28)=0;
    rLocalLHS(30,29)=0;
    rLocalLHS(30,30)=0;
    rLocalLHS(30,31)=0;
    rLocalLHS(30,32)=0;
    rLocalLHS(31,0)=0;
    rLocalLHS(31,1)=clhs3;
    rLocalLHS(31,2)=0;
    rLocalLHS(31,3)=0;
    rLocalLHS(31,4)=clhs7;
    rLocalLHS(31,5)=0;
    rLocalLHS(31,6)=0;
    rLocalLHS(31,7)=clhs11;
    rLocalLHS(31,8)=0;
    rLocalLHS(31,9)=0;
    rLocalLHS(31,10)=DOperator(3,0);
    rLocalLHS(31,11)=0;
    rLocalLHS(31,12)=0;
    rLocalLHS(31,13)=DOperator(3,1);
    rLocalLHS(31,14)=0;
    rLocalLHS(31,15)=0;
    rLocalLHS(31,16)=DOperator(3,2);
    rLocalLHS(31,17)=0;
    rLocalLHS(31,18)=0;
    rLocalLHS(31,19)=DOperator(3,3);
    rLocalLHS(31,20)=0;
    rLocalLHS(31,21)=0;
    rLocalLHS(31,22)=0;
    rLocalLHS(31,23)=0;
    rLocalLHS(31,24)=0;
    rLocalLHS(31,25)=0;
    rLocalLHS(31,26)=0;
    rLocalLHS(31,27)=0;
    rLocalLHS(31,28)=0;
    rLocalLHS(31,29)=0;
    rLocalLHS(31,30)=0;
    rLocalLHS(31,31)=0;
    rLocalLHS(31,32)=0;
    rLocalLHS(32,0)=0;
    rLocalLHS(32,1)=0;
    rLocalLHS(32,2)=clhs3;
    rLocalLHS(32,3)=0;
    rLocalLHS(32,4)=0;
    rLocalLHS(32,5)=clhs7;
    rLocalLHS(32,6)=0;
    rLocalLHS(32,7)=0;
    rLocalLHS(32,8)=clhs11;
    rLocalLHS(32,9)=0;
    rLocalLHS(32,10)=0;
    rLocalLHS(32,11)=DOperator(3,0);
    rLocalLHS(32,12)=0;
    rLocalLHS(32,13)=0;
    rLocalLHS(32,14)=DOperator(3,1);
    rLocalLHS(32,15)=0;
    rLocalLHS(32,16)=0;
    rLocalLHS(32,17)=DOperator(3,2);
    rLocalLHS(32,18)=0;
    rLocalLHS(32,19)=0;
    rLocalLHS(32,20)=DOperator(3,3);
    rLocalLHS(32,21)=0;
    rLocalLHS(32,22)=0;
    rLocalLHS(32,23)=0;
    rLocalLHS(32,24)=0;
    rLocalLHS(32,25)=0;
    rLocalLHS(32,26)=0;
    rLocalLHS(32,27)=0;
    rLocalLHS(32,28)=0;
    rLocalLHS(32,29)=0;
    rLocalLHS(32,30)=0;
    rLocalLHS(32,31)=0;
    rLocalLHS(32,32)=0;

}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<2,3, 3>::CalculateLocalRHS<MeshTyingMortarCondition<2,3,3>::ScalarValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 2, ScalarValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 2, ScalarValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 2, ScalarValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rLocalRHS[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rLocalRHS[2]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rLocalRHS[3]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rLocalRHS[4]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rLocalRHS[5]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<2,3, 3>::CalculateLocalRHS<MeshTyingMortarCondition<2,3,3>::Vector2DValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector2DValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 2, Vector2DValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 2, Vector2DValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 2, Vector2DValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rLocalRHS[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1);
    rLocalRHS[2]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rLocalRHS[3]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1);
    rLocalRHS[4]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rLocalRHS[5]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1));
    rLocalRHS[6]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rLocalRHS[7]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1));
    rLocalRHS[8]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rLocalRHS[9]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1);
    rLocalRHS[10]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);
    rLocalRHS[11]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<2,4, 4>::CalculateLocalRHS<MeshTyingMortarCondition<2,4,4>::ScalarValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 2, ScalarValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 2, ScalarValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 2, ScalarValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rLocalRHS[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rLocalRHS[2]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rLocalRHS[3]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rLocalRHS[4]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rLocalRHS[5]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<2,4, 4>::CalculateLocalRHS<MeshTyingMortarCondition<2,4,4>::Vector2DValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector2DValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 2, Vector2DValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 2, Vector2DValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 2, Vector2DValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rLocalRHS[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1);
    rLocalRHS[2]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rLocalRHS[3]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1);
    rLocalRHS[4]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rLocalRHS[5]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1));
    rLocalRHS[6]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rLocalRHS[7]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1));
    rLocalRHS[8]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rLocalRHS[9]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1);
    rLocalRHS[10]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);
    rLocalRHS[11]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,4, 4>::CalculateLocalRHS<MeshTyingMortarCondition<3,4,4>::ScalarValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 3, ScalarValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 3, ScalarValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 3, ScalarValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0);
    rLocalRHS[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0);
    rLocalRHS[2]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0);
    rLocalRHS[3]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0));
    rLocalRHS[4]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0));
    rLocalRHS[5]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0));
    rLocalRHS[6]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0);
    rLocalRHS[7]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0);
    rLocalRHS[8]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,4, 4>::CalculateLocalRHS<MeshTyingMortarCondition<3,4,4>::Vector3DValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 3, Vector3DValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 3, Vector3DValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 3, Vector3DValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0);
    rLocalRHS[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1) + MOperator(2,0)*lm(2,1);
    rLocalRHS[2]=MOperator(0,0)*lm(0,2) + MOperator(1,0)*lm(1,2) + MOperator(2,0)*lm(2,2);
    rLocalRHS[3]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0);
    rLocalRHS[4]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1) + MOperator(2,1)*lm(2,1);
    rLocalRHS[5]=MOperator(0,1)*lm(0,2) + MOperator(1,1)*lm(1,2) + MOperator(2,1)*lm(2,2);
    rLocalRHS[6]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0);
    rLocalRHS[7]=MOperator(0,2)*lm(0,1) + MOperator(1,2)*lm(1,1) + MOperator(2,2)*lm(2,1);
    rLocalRHS[8]=MOperator(0,2)*lm(0,2) + MOperator(1,2)*lm(1,2) + MOperator(2,2)*lm(2,2);
    rLocalRHS[9]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0));
    rLocalRHS[10]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1) + DOperator(2,0)*lm(2,1));
    rLocalRHS[11]=-(DOperator(0,0)*lm(0,2) + DOperator(1,0)*lm(1,2) + DOperator(2,0)*lm(2,2));
    rLocalRHS[12]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0));
    rLocalRHS[13]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1) + DOperator(2,1)*lm(2,1));
    rLocalRHS[14]=-(DOperator(0,1)*lm(0,2) + DOperator(1,1)*lm(1,2) + DOperator(2,1)*lm(2,2));
    rLocalRHS[15]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0));
    rLocalRHS[16]=-(DOperator(0,2)*lm(0,1) + DOperator(1,2)*lm(1,1) + DOperator(2,2)*lm(2,1));
    rLocalRHS[17]=-(DOperator(0,2)*lm(0,2) + DOperator(1,2)*lm(1,2) + DOperator(2,2)*lm(2,2));
    rLocalRHS[18]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0);
    rLocalRHS[19]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) - DOperator(0,2)*u1(2,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1) + MOperator(0,2)*u2(2,1);
    rLocalRHS[20]=-DOperator(0,0)*u1(0,2) - DOperator(0,1)*u1(1,2) - DOperator(0,2)*u1(2,2) + MOperator(0,0)*u2(0,2) + MOperator(0,1)*u2(1,2) + MOperator(0,2)*u2(2,2);
    rLocalRHS[21]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0);
    rLocalRHS[22]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) - DOperator(1,2)*u1(2,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1) + MOperator(1,2)*u2(2,1);
    rLocalRHS[23]=-DOperator(1,0)*u1(0,2) - DOperator(1,1)*u1(1,2) - DOperator(1,2)*u1(2,2) + MOperator(1,0)*u2(0,2) + MOperator(1,1)*u2(1,2) + MOperator(1,2)*u2(2,2);
    rLocalRHS[24]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0);
    rLocalRHS[25]=-DOperator(2,0)*u1(0,1) - DOperator(2,1)*u1(1,1) - DOperator(2,2)*u1(2,1) + MOperator(2,0)*u2(0,1) + MOperator(2,1)*u2(1,1) + MOperator(2,2)*u2(2,1);
    rLocalRHS[26]=-DOperator(2,0)*u1(0,2) - DOperator(2,1)*u1(1,2) - DOperator(2,2)*u1(2,2) + MOperator(2,0)*u2(0,2) + MOperator(2,1)*u2(1,2) + MOperator(2,2)*u2(2,2);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,8, 8>::CalculateLocalRHS<MeshTyingMortarCondition<3,8,8>::ScalarValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 4, ScalarValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 4, ScalarValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 4, ScalarValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0) + MOperator(3,0)*lm(3,0);
    rLocalRHS[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0) + MOperator(3,1)*lm(3,0);
    rLocalRHS[2]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0) + MOperator(3,2)*lm(3,0);
    rLocalRHS[3]=MOperator(0,3)*lm(0,0) + MOperator(1,3)*lm(1,0) + MOperator(2,3)*lm(2,0) + MOperator(3,3)*lm(3,0);
    rLocalRHS[4]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0) + DOperator(3,0)*lm(3,0));
    rLocalRHS[5]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0) + DOperator(3,1)*lm(3,0));
    rLocalRHS[6]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0) + DOperator(3,2)*lm(3,0));
    rLocalRHS[7]=-(DOperator(0,3)*lm(0,0) + DOperator(1,3)*lm(1,0) + DOperator(2,3)*lm(2,0) + DOperator(3,3)*lm(3,0));
    rLocalRHS[8]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) - DOperator(0,3)*u1(3,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0) + MOperator(0,3)*u2(3,0);
    rLocalRHS[9]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) - DOperator(1,3)*u1(3,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0) + MOperator(1,3)*u2(3,0);
    rLocalRHS[10]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) - DOperator(2,3)*u1(3,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0) + MOperator(2,3)*u2(3,0);
    rLocalRHS[11]=-DOperator(3,0)*u1(0,0) - DOperator(3,1)*u1(1,0) - DOperator(3,2)*u1(2,0) - DOperator(3,3)*u1(3,0) + MOperator(3,0)*u2(0,0) + MOperator(3,1)*u2(1,0) + MOperator(3,2)*u2(2,0) + MOperator(3,3)*u2(3,0);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,8, 8>::CalculateLocalRHS<MeshTyingMortarCondition<3,8,8>::Vector3DValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 4, Vector3DValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 4, Vector3DValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 4, Vector3DValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0) + MOperator(3,0)*lm(3,0);
    rLocalRHS[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1) + MOperator(2,0)*lm(2,1) + MOperator(3,0)*lm(3,1);
    rLocalRHS[2]=MOperator(0,0)*lm(0,2) + MOperator(1,0)*lm(1,2) + MOperator(2,0)*lm(2,2) + MOperator(3,0)*lm(3,2);
    rLocalRHS[3]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0) + MOperator(3,1)*lm(3,0);
    rLocalRHS[4]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1) + MOperator(2,1)*lm(2,1) + MOperator(3,1)*lm(3,1);
    rLocalRHS[5]=MOperator(0,1)*lm(0,2) + MOperator(1,1)*lm(1,2) + MOperator(2,1)*lm(2,2) + MOperator(3,1)*lm(3,2);
    rLocalRHS[6]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0) + MOperator(3,2)*lm(3,0);
    rLocalRHS[7]=MOperator(0,2)*lm(0,1) + MOperator(1,2)*lm(1,1) + MOperator(2,2)*lm(2,1) + MOperator(3,2)*lm(3,1);
    rLocalRHS[8]=MOperator(0,2)*lm(0,2) + MOperator(1,2)*lm(1,2) + MOperator(2,2)*lm(2,2) + MOperator(3,2)*lm(3,2);
    rLocalRHS[9]=MOperator(0,3)*lm(0,0) + MOperator(1,3)*lm(1,0) + MOperator(2,3)*lm(2,0) + MOperator(3,3)*lm(3,0);
    rLocalRHS[10]=MOperator(0,3)*lm(0,1) + MOperator(1,3)*lm(1,1) + MOperator(2,3)*lm(2,1) + MOperator(3,3)*lm(3,1);
    rLocalRHS[11]=MOperator(0,3)*lm(0,2) + MOperator(1,3)*lm(1,2) + MOperator(2,3)*lm(2,2) + MOperator(3,3)*lm(3,2);
    rLocalRHS[12]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0) + DOperator(3,0)*lm(3,0));
    rLocalRHS[13]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1) + DOperator(2,0)*lm(2,1) + DOperator(3,0)*lm(3,1));
    rLocalRHS[14]=-(DOperator(0,0)*lm(0,2) + DOperator(1,0)*lm(1,2) + DOperator(2,0)*lm(2,2) + DOperator(3,0)*lm(3,2));
    rLocalRHS[15]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0) + DOperator(3,1)*lm(3,0));
    rLocalRHS[16]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1) + DOperator(2,1)*lm(2,1) + DOperator(3,1)*lm(3,1));
    rLocalRHS[17]=-(DOperator(0,1)*lm(0,2) + DOperator(1,1)*lm(1,2) + DOperator(2,1)*lm(2,2) + DOperator(3,1)*lm(3,2));
    rLocalRHS[18]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0) + DOperator(3,2)*lm(3,0));
    rLocalRHS[19]=-(DOperator(0,2)*lm(0,1) + DOperator(1,2)*lm(1,1) + DOperator(2,2)*lm(2,1) + DOperator(3,2)*lm(3,1));
    rLocalRHS[20]=-(DOperator(0,2)*lm(0,2) + DOperator(1,2)*lm(1,2) + DOperator(2,2)*lm(2,2) + DOperator(3,2)*lm(3,2));
    rLocalRHS[21]=-(DOperator(0,3)*lm(0,0) + DOperator(1,3)*lm(1,0) + DOperator(2,3)*lm(2,0) + DOperator(3,3)*lm(3,0));
    rLocalRHS[22]=-(DOperator(0,3)*lm(0,1) + DOperator(1,3)*lm(1,1) + DOperator(2,3)*lm(2,1) + DOperator(3,3)*lm(3,1));
    rLocalRHS[23]=-(DOperator(0,3)*lm(0,2) + DOperator(1,3)*lm(1,2) + DOperator(2,3)*lm(2,2) + DOperator(3,3)*lm(3,2));
    rLocalRHS[24]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) - DOperator(0,3)*u1(3,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0) + MOperator(0,3)*u2(3,0);
    rLocalRHS[25]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) - DOperator(0,2)*u1(2,1) - DOperator(0,3)*u1(3,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1) + MOperator(0,2)*u2(2,1) + MOperator(0,3)*u2(3,1);
    rLocalRHS[26]=-DOperator(0,0)*u1(0,2) - DOperator(0,1)*u1(1,2) - DOperator(0,2)*u1(2,2) - DOperator(0,3)*u1(3,2) + MOperator(0,0)*u2(0,2) + MOperator(0,1)*u2(1,2) + MOperator(0,2)*u2(2,2) + MOperator(0,3)*u2(3,2);
    rLocalRHS[27]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) - DOperator(1,3)*u1(3,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0) + MOperator(1,3)*u2(3,0);
    rLocalRHS[28]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) - DOperator(1,2)*u1(2,1) - DOperator(1,3)*u1(3,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1) + MOperator(1,2)*u2(2,1) + MOperator(1,3)*u2(3,1);
    rLocalRHS[29]=-DOperator(1,0)*u1(0,2) - DOperator(1,1)*u1(1,2) - DOperator(1,2)*u1(2,2) - DOperator(1,3)*u1(3,2) + MOperator(1,0)*u2(0,2) + MOperator(1,1)*u2(1,2) + MOperator(1,2)*u2(2,2) + MOperator(1,3)*u2(3,2);
    rLocalRHS[30]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) - DOperator(2,3)*u1(3,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0) + MOperator(2,3)*u2(3,0);
    rLocalRHS[31]=-DOperator(2,0)*u1(0,1) - DOperator(2,1)*u1(1,1) - DOperator(2,2)*u1(2,1) - DOperator(2,3)*u1(3,1) + MOperator(2,0)*u2(0,1) + MOperator(2,1)*u2(1,1) + MOperator(2,2)*u2(2,1) + MOperator(2,3)*u2(3,1);
    rLocalRHS[32]=-DOperator(2,0)*u1(0,2) - DOperator(2,1)*u1(1,2) - DOperator(2,2)*u1(2,2) - DOperator(2,3)*u1(3,2) + MOperator(2,0)*u2(0,2) + MOperator(2,1)*u2(1,2) + MOperator(2,2)*u2(2,2) + MOperator(2,3)*u2(3,2);
    rLocalRHS[33]=-DOperator(3,0)*u1(0,0) - DOperator(3,1)*u1(1,0) - DOperator(3,2)*u1(2,0) - DOperator(3,3)*u1(3,0) + MOperator(3,0)*u2(0,0) + MOperator(3,1)*u2(1,0) + MOperator(3,2)*u2(2,0) + MOperator(3,3)*u2(3,0);
    rLocalRHS[34]=-DOperator(3,0)*u1(0,1) - DOperator(3,1)*u1(1,1) - DOperator(3,2)*u1(2,1) - DOperator(3,3)*u1(3,1) + MOperator(3,0)*u2(0,1) + MOperator(3,1)*u2(1,1) + MOperator(3,2)*u2(2,1) + MOperator(3,3)*u2(3,1);
    rLocalRHS[35]=-DOperator(3,0)*u1(0,2) - DOperator(3,1)*u1(1,2) - DOperator(3,2)*u1(2,2) - DOperator(3,3)*u1(3,2) + MOperator(3,0)*u2(0,2) + MOperator(3,1)*u2(1,2) + MOperator(3,2)*u2(2,2) + MOperator(3,3)*u2(3,2);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,4, 8>::CalculateLocalRHS<MeshTyingMortarCondition<3,4,8>::ScalarValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 3, ScalarValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 4, ScalarValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 3, ScalarValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 3, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0);
    rLocalRHS[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0);
    rLocalRHS[2]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0);
    rLocalRHS[3]=MOperator(0,3)*lm(0,0) + MOperator(1,3)*lm(1,0) + MOperator(2,3)*lm(2,0);
    rLocalRHS[4]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0));
    rLocalRHS[5]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0));
    rLocalRHS[6]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0));
    rLocalRHS[7]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0) + MOperator(0,3)*u2(3,0);
    rLocalRHS[8]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0) + MOperator(1,3)*u2(3,0);
    rLocalRHS[9]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0) + MOperator(2,3)*u2(3,0);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,4, 8>::CalculateLocalRHS<MeshTyingMortarCondition<3,4,8>::Vector3DValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 3, Vector3DValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 4, Vector3DValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 3, Vector3DValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 3, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0);
    rLocalRHS[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1) + MOperator(2,0)*lm(2,1);
    rLocalRHS[2]=MOperator(0,0)*lm(0,2) + MOperator(1,0)*lm(1,2) + MOperator(2,0)*lm(2,2);
    rLocalRHS[3]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0);
    rLocalRHS[4]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1) + MOperator(2,1)*lm(2,1);
    rLocalRHS[5]=MOperator(0,1)*lm(0,2) + MOperator(1,1)*lm(1,2) + MOperator(2,1)*lm(2,2);
    rLocalRHS[6]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0);
    rLocalRHS[7]=MOperator(0,2)*lm(0,1) + MOperator(1,2)*lm(1,1) + MOperator(2,2)*lm(2,1);
    rLocalRHS[8]=MOperator(0,2)*lm(0,2) + MOperator(1,2)*lm(1,2) + MOperator(2,2)*lm(2,2);
    rLocalRHS[9]=MOperator(0,3)*lm(0,0) + MOperator(1,3)*lm(1,0) + MOperator(2,3)*lm(2,0);
    rLocalRHS[10]=MOperator(0,3)*lm(0,1) + MOperator(1,3)*lm(1,1) + MOperator(2,3)*lm(2,1);
    rLocalRHS[11]=MOperator(0,3)*lm(0,2) + MOperator(1,3)*lm(1,2) + MOperator(2,3)*lm(2,2);
    rLocalRHS[12]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0));
    rLocalRHS[13]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1) + DOperator(2,0)*lm(2,1));
    rLocalRHS[14]=-(DOperator(0,0)*lm(0,2) + DOperator(1,0)*lm(1,2) + DOperator(2,0)*lm(2,2));
    rLocalRHS[15]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0));
    rLocalRHS[16]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1) + DOperator(2,1)*lm(2,1));
    rLocalRHS[17]=-(DOperator(0,1)*lm(0,2) + DOperator(1,1)*lm(1,2) + DOperator(2,1)*lm(2,2));
    rLocalRHS[18]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0));
    rLocalRHS[19]=-(DOperator(0,2)*lm(0,1) + DOperator(1,2)*lm(1,1) + DOperator(2,2)*lm(2,1));
    rLocalRHS[20]=-(DOperator(0,2)*lm(0,2) + DOperator(1,2)*lm(1,2) + DOperator(2,2)*lm(2,2));
    rLocalRHS[21]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0) + MOperator(0,3)*u2(3,0);
    rLocalRHS[22]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) - DOperator(0,2)*u1(2,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1) + MOperator(0,2)*u2(2,1) + MOperator(0,3)*u2(3,1);
    rLocalRHS[23]=-DOperator(0,0)*u1(0,2) - DOperator(0,1)*u1(1,2) - DOperator(0,2)*u1(2,2) + MOperator(0,0)*u2(0,2) + MOperator(0,1)*u2(1,2) + MOperator(0,2)*u2(2,2) + MOperator(0,3)*u2(3,2);
    rLocalRHS[24]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0) + MOperator(1,3)*u2(3,0);
    rLocalRHS[25]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) - DOperator(1,2)*u1(2,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1) + MOperator(1,2)*u2(2,1) + MOperator(1,3)*u2(3,1);
    rLocalRHS[26]=-DOperator(1,0)*u1(0,2) - DOperator(1,1)*u1(1,2) - DOperator(1,2)*u1(2,2) + MOperator(1,0)*u2(0,2) + MOperator(1,1)*u2(1,2) + MOperator(1,2)*u2(2,2) + MOperator(1,3)*u2(3,2);
    rLocalRHS[27]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0) + MOperator(2,3)*u2(3,0);
    rLocalRHS[28]=-DOperator(2,0)*u1(0,1) - DOperator(2,1)*u1(1,1) - DOperator(2,2)*u1(2,1) + MOperator(2,0)*u2(0,1) + MOperator(2,1)*u2(1,1) + MOperator(2,2)*u2(2,1) + MOperator(2,3)*u2(3,1);
    rLocalRHS[29]=-DOperator(2,0)*u1(0,2) - DOperator(2,1)*u1(1,2) - DOperator(2,2)*u1(2,2) + MOperator(2,0)*u2(0,2) + MOperator(2,1)*u2(1,2) + MOperator(2,2)*u2(2,2) + MOperator(2,3)*u2(3,2);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,8, 4>::CalculateLocalRHS<MeshTyingMortarCondition<3,8,4>::ScalarValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<ScalarValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 4, ScalarValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 3, ScalarValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 4, ScalarValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 4, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0) + MOperator(3,0)*lm(3,0);
    rLocalRHS[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0) + MOperator(3,1)*lm(3,0);
    rLocalRHS[2]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0) + MOperator(3,2)*lm(3,0);
    rLocalRHS[3]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0) + DOperator(3,0)*lm(3,0));
    rLocalRHS[4]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0) + DOperator(3,1)*lm(3,0));
    rLocalRHS[5]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0) + DOperator(3,2)*lm(3,0));
    rLocalRHS[6]=-(DOperator(0,3)*lm(0,0) + DOperator(1,3)*lm(1,0) + DOperator(2,3)*lm(2,0) + DOperator(3,3)*lm(3,0));
    rLocalRHS[7]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) - DOperator(0,3)*u1(3,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0);
    rLocalRHS[8]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) - DOperator(1,3)*u1(3,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0);
    rLocalRHS[9]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) - DOperator(2,3)*u1(3,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0);
    rLocalRHS[10]=-DOperator(3,0)*u1(0,0) - DOperator(3,1)*u1(1,0) - DOperator(3,2)*u1(2,0) - DOperator(3,3)*u1(3,0) + MOperator(3,0)*u2(0,0) + MOperator(3,1)*u2(1,0) + MOperator(3,2)*u2(2,0);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
void MeshTyingMortarCondition<3,8, 4>::CalculateLocalRHS<MeshTyingMortarCondition<3,8,4>::Vector3DValue>(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData<Vector3DValue>& rDofData
    )
{
    // Initialize values
    const BoundedMatrix<double, 4, Vector3DValue> u1 = rDofData.u1;
    const BoundedMatrix<double, 3, Vector3DValue> u2 = rDofData.u2;

    const BoundedMatrix<double, 4, Vector3DValue> lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const BoundedMatrix<double, 4, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;


    rLocalRHS[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0) + MOperator(3,0)*lm(3,0);
    rLocalRHS[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1) + MOperator(2,0)*lm(2,1) + MOperator(3,0)*lm(3,1);
    rLocalRHS[2]=MOperator(0,0)*lm(0,2) + MOperator(1,0)*lm(1,2) + MOperator(2,0)*lm(2,2) + MOperator(3,0)*lm(3,2);
    rLocalRHS[3]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0) + MOperator(3,1)*lm(3,0);
    rLocalRHS[4]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1) + MOperator(2,1)*lm(2,1) + MOperator(3,1)*lm(3,1);
    rLocalRHS[5]=MOperator(0,1)*lm(0,2) + MOperator(1,1)*lm(1,2) + MOperator(2,1)*lm(2,2) + MOperator(3,1)*lm(3,2);
    rLocalRHS[6]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0) + MOperator(3,2)*lm(3,0);
    rLocalRHS[7]=MOperator(0,2)*lm(0,1) + MOperator(1,2)*lm(1,1) + MOperator(2,2)*lm(2,1) + MOperator(3,2)*lm(3,1);
    rLocalRHS[8]=MOperator(0,2)*lm(0,2) + MOperator(1,2)*lm(1,2) + MOperator(2,2)*lm(2,2) + MOperator(3,2)*lm(3,2);
    rLocalRHS[9]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0) + DOperator(3,0)*lm(3,0));
    rLocalRHS[10]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1) + DOperator(2,0)*lm(2,1) + DOperator(3,0)*lm(3,1));
    rLocalRHS[11]=-(DOperator(0,0)*lm(0,2) + DOperator(1,0)*lm(1,2) + DOperator(2,0)*lm(2,2) + DOperator(3,0)*lm(3,2));
    rLocalRHS[12]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0) + DOperator(3,1)*lm(3,0));
    rLocalRHS[13]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1) + DOperator(2,1)*lm(2,1) + DOperator(3,1)*lm(3,1));
    rLocalRHS[14]=-(DOperator(0,1)*lm(0,2) + DOperator(1,1)*lm(1,2) + DOperator(2,1)*lm(2,2) + DOperator(3,1)*lm(3,2));
    rLocalRHS[15]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0) + DOperator(3,2)*lm(3,0));
    rLocalRHS[16]=-(DOperator(0,2)*lm(0,1) + DOperator(1,2)*lm(1,1) + DOperator(2,2)*lm(2,1) + DOperator(3,2)*lm(3,1));
    rLocalRHS[17]=-(DOperator(0,2)*lm(0,2) + DOperator(1,2)*lm(1,2) + DOperator(2,2)*lm(2,2) + DOperator(3,2)*lm(3,2));
    rLocalRHS[18]=-(DOperator(0,3)*lm(0,0) + DOperator(1,3)*lm(1,0) + DOperator(2,3)*lm(2,0) + DOperator(3,3)*lm(3,0));
    rLocalRHS[19]=-(DOperator(0,3)*lm(0,1) + DOperator(1,3)*lm(1,1) + DOperator(2,3)*lm(2,1) + DOperator(3,3)*lm(3,1));
    rLocalRHS[20]=-(DOperator(0,3)*lm(0,2) + DOperator(1,3)*lm(1,2) + DOperator(2,3)*lm(2,2) + DOperator(3,3)*lm(3,2));
    rLocalRHS[21]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) - DOperator(0,3)*u1(3,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0);
    rLocalRHS[22]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) - DOperator(0,2)*u1(2,1) - DOperator(0,3)*u1(3,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1) + MOperator(0,2)*u2(2,1);
    rLocalRHS[23]=-DOperator(0,0)*u1(0,2) - DOperator(0,1)*u1(1,2) - DOperator(0,2)*u1(2,2) - DOperator(0,3)*u1(3,2) + MOperator(0,0)*u2(0,2) + MOperator(0,1)*u2(1,2) + MOperator(0,2)*u2(2,2);
    rLocalRHS[24]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) - DOperator(1,3)*u1(3,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0);
    rLocalRHS[25]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) - DOperator(1,2)*u1(2,1) - DOperator(1,3)*u1(3,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1) + MOperator(1,2)*u2(2,1);
    rLocalRHS[26]=-DOperator(1,0)*u1(0,2) - DOperator(1,1)*u1(1,2) - DOperator(1,2)*u1(2,2) - DOperator(1,3)*u1(3,2) + MOperator(1,0)*u2(0,2) + MOperator(1,1)*u2(1,2) + MOperator(1,2)*u2(2,2);
    rLocalRHS[27]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) - DOperator(2,3)*u1(3,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0);
    rLocalRHS[28]=-DOperator(2,0)*u1(0,1) - DOperator(2,1)*u1(1,1) - DOperator(2,2)*u1(2,1) - DOperator(2,3)*u1(3,1) + MOperator(2,0)*u2(0,1) + MOperator(2,1)*u2(1,1) + MOperator(2,2)*u2(2,1);
    rLocalRHS[29]=-DOperator(2,0)*u1(0,2) - DOperator(2,1)*u1(1,2) - DOperator(2,2)*u1(2,2) - DOperator(2,3)*u1(3,2) + MOperator(2,0)*u2(0,2) + MOperator(2,1)*u2(1,2) + MOperator(2,2)*u2(2,2);
    rLocalRHS[30]=-DOperator(3,0)*u1(0,0) - DOperator(3,1)*u1(1,0) - DOperator(3,2)*u1(2,0) - DOperator(3,3)*u1(3,0) + MOperator(3,0)*u2(0,0) + MOperator(3,1)*u2(1,0) + MOperator(3,2)*u2(2,0);
    rLocalRHS[31]=-DOperator(3,0)*u1(0,1) - DOperator(3,1)*u1(1,1) - DOperator(3,2)*u1(2,1) - DOperator(3,3)*u1(3,1) + MOperator(3,0)*u2(0,1) + MOperator(3,1)*u2(1,1) + MOperator(3,2)*u2(2,1);
    rLocalRHS[32]=-DOperator(3,0)*u1(0,2) - DOperator(3,1)*u1(1,2) - DOperator(3,2)*u1(2,2) - DOperator(3,3)*u1(3,2) + MOperator(3,0)*u2(0,2) + MOperator(3,1)*u2(1,2) + MOperator(3,2)*u2(2,2);

}


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
        const std::string& variable_name = (mArray1DVariables[0]).Name();
        const Array1DComponentsType& var_x = KratosComponents<Array1DComponentsType>::Get(variable_name+"_X");
        const Array1DComponentsType& var_y = KratosComponents<Array1DComponentsType>::Get(variable_name+"_Y");
        const Array1DComponentsType& var_z = KratosComponents<Array1DComponentsType>::Get(variable_name+"_Z");
        for ( IndexType i_master = 0; i_master < NumNodesMaster; ++i_master ) {
            NodeType& r_master_node = r_current_master[i_master];
            rResult[index++] = r_master_node.GetDof( var_x ).EquationId( );
            rResult[index++] = r_master_node.GetDof( var_y ).EquationId( );
            if (TDim == 3) rResult[index++] = r_master_node.GetDof( var_z ).EquationId( );
        }
    }

    // Slave Nodes DoF Equation IDs
    if (tensor_value == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            rResult[index++] = this->GetGeometry()[i_slave].GetDof( mDoubleVariables[0] ).EquationId( );
        }
    } else {
        const std::string& variable_name = (mArray1DVariables[0]).Name();
        const Array1DComponentsType& var_x = KratosComponents<Array1DComponentsType>::Get(variable_name+"_X");
        const Array1DComponentsType& var_y = KratosComponents<Array1DComponentsType>::Get(variable_name+"_Y");
        const Array1DComponentsType& var_z = KratosComponents<Array1DComponentsType>::Get(variable_name+"_Z");
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rResult[index++] = slave_node.GetDof( var_x ).EquationId( );
            rResult[index++] = slave_node.GetDof( var_y ).EquationId( );
            if (TDim == 3) rResult[index++] = slave_node.GetDof( var_z ).EquationId( );
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
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
            if (TDim == 3) rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
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
        const std::string& variable_name = (mArray1DVariables[0]).Name();
        const Array1DComponentsType& var_x = KratosComponents< Array1DComponentsType>::Get(variable_name+"_X");
        const Array1DComponentsType& var_y = KratosComponents< Array1DComponentsType>::Get(variable_name+"_Y");
        const Array1DComponentsType& var_z = KratosComponents< Array1DComponentsType>::Get(variable_name+"_Z");
        for ( IndexType i_master = 0; i_master < NumNodesMaster; ++i_master ) {
            NodeType& master_node = current_master[i_master];
            rConditionalDofList[index++] = master_node.pGetDof( var_x );
            rConditionalDofList[index++] = master_node.pGetDof( var_y );
            if (TDim == 3) rConditionalDofList[index++] = master_node.pGetDof( var_z );
        }
    }

    // Slave Nodes DoF Equation IDs
    if (tensor_value == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            rConditionalDofList[index++] = this->GetGeometry()[i_slave].pGetDof( mDoubleVariables[0] );
        }
    } else {
        const std::string& variable_name = (mArray1DVariables[0]).Name();
        const Array1DComponentsType& var_x = KratosComponents< Array1DComponentsType>::Get(variable_name+"_X");
        const Array1DComponentsType& var_y = KratosComponents< Array1DComponentsType>::Get(variable_name+"_Y");
        const Array1DComponentsType& var_z = KratosComponents< Array1DComponentsType>::Get(variable_name+"_Z");
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rConditionalDofList[index++] = slave_node.pGetDof( var_x );
            rConditionalDofList[index++] = slave_node.pGetDof( var_y );
            if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( var_z );
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
            rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
            if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
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
