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
#include <algorithm>

// External includes

// Project includes
#include "custom_conditions/mesh_tying_mortar_condition.h"
/* Mortar includes */
#include "custom_conditions/mesh_tying_mortar_condition.h"

namespace Kratos 
{
/**
 * Flags related to the condition computation 
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor> >( NewId, pGeom, pProperties, pMasterGeom);
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::~MeshTyingMortarCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Initialize( ) 
{
    KRATOS_TRY;
    
    BaseType::Initialize();
    
    mIntegrationOrder = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    
    // The slave geometry
    GeometryType& slave_geometry = this->GetGeometry();
    const array_1d<double, 3>& normal_slave = this->GetValue(NORMAL);
    
    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Create the current DoF data
    DofData rDofData;
    
    // The master geometry
    GeometryType& master_geometry = this->GetPairedGeometry();
    const array_1d<double, 3>& normal_master = this->GetValue(PAIRED_NORMAL);
    // Initialize general variables for the current master element
    rVariables.Initialize();
    
    // Initialize the mortar operators
    mrThisMortarConditionMatrices.Initialize();
    
    // We call the exact integration utility
    IntegrationUtility integration_utility = IntegrationUtility (mIntegrationOrder);
    
    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, normal_master, conditions_points_slave);
    
    double integration_area;
    integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);
    
    if ((is_inside == true) && ((integration_area/slave_geometry.Area()) > 1.0e-3 * slave_geometry.Area())) {
        IntegrationMethod this_integration_method = GetIntegrationMethod();
        
        // Initialize general variables for the current master element
        rVariables.Initialize();
        
        // Initialize the mortar operators
        mrThisMortarConditionMatrices.Initialize();
        
        const bool dual_LM = CalculateAe(normal_master, rDofData, rVariables, conditions_points_slave, this_integration_method);
            
        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array[i_node] = PointType::Pointer( new PointType(global_point) );
            }
            
            DecompositionType decomp_geom( points_array );
            
            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
            
            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                
                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates 
                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                    
                    // Calculate the kinematic variables
                    this->CalculateKinematics( rVariables, rDofData, normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);
                    
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

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // TODO: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // Calculation flags
    mCalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );
    mCalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true );
    
    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
            rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
    
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );
    
    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculation flags
    mCalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );
    mCalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, false);

    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
        rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
    
    // Creating an auxiliar vector
    VectorType aux_right_hand_side_vector = Vector();

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculation flags
    mCalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, false);
    mCalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true);

    // Creating an auxiliar matrix
    MatrixType aux_left_hand_side_matrix = Matrix();
    
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );

    // Calculate condition system
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateMassMatrix( 
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

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateDampingMatrix( 
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

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim, TNumNodesElem, TTensor>::CalculateConditionSystem( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
        
    // Create the current DoF data
    DofData rDofData;
      
    // Initialize the DoF data
    this->InitializeDofData(rDofData);
    
    // Update slave element info
    rDofData.UpdateMasterPair(this->GetPairedGeometry());
    
    // Assemble of the matrix is required
    if ( mCalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX ) ) {
        // Calculate the local contribution
        this->CalculateLocalLHS(rLeftHandSideMatrix, mrThisMortarConditionMatrices, rDofData);
    }
    
    // Assemble of the vector is required
    if ( mCalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR )) {
        // Calculate the local contribution
        this->CalculateLocalRHS( rRightHandSideVector, mrThisMortarConditionMatrices, rDofData);
    }
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
bool MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateAe(
    const array_1d<double, 3>& NormalMaster,
    DofData& rDofData,
    GeneralVariables& rVariables,
    ConditionArrayListType& ConditionsPointsSlave,
    IntegrationMethod ThisIntegrationMethod
    )
{
    // We initilize the Ae components
    AeData rAeData;
    rAeData.Initialize();
    
    rDofData.InitializeAeComponents();
    
    // The slave geometry
    GeometryType& slave_geometry = GetGeometry();
    
    // Initialize general variables for the current master element
    rVariables.Initialize();
        
    // Calculating the proportion between the integrated area and segment area
    for (IndexType i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            slave_geometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = PointType::Pointer( new PointType(global_point) );
        }
        
        DecompositionType decomp_geom( points_array );
        
        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
        
        if (bad_shape == false) {
            const GeometryType::IntegrationPointsArrayType integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );
            
            // Integrating the mortar operators
            for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                // We compute the local coordinates 
                const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                PointType local_point_parent;
                PointType gp_global;
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                
                // Calculate the kinematic variables
                // We compute the current configuration
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                
                this->CalculateKinematics( rVariables, rDofData, NormalMaster, local_point_decomp, local_point_parent, decomp_geom, false);
                
                // Integrate
                const double integration_weight = integration_points_slave[point_number].Weight();
        
                rAeData.CalculateAeComponents(rVariables, integration_weight);
            }
        }
    }
    
    return rAeData.CalculateAe(rDofData.Ae);
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeDofData(DofData& rDofData)
{
    // Slave element info
    rDofData.Initialize(GetGeometry());
    
    if (TTensor == 1) {
        for (IndexType i_node = 0; i_node < NumNodes; i_node++) {
            const double value = GetGeometry()[i_node].FastGetSolutionStepValue(TEMPERATURE);
            const double lm = GetGeometry()[i_node].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
            rDofData.u1(i_node, 0) = value;
            rDofData.LagrangeMultipliers(i_node, 0) = lm;
        }
    } else {
        for (IndexType i_node = 0; i_node < NumNodes; i_node++) {
            const array_1d<double, 3>& value = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& lm = GetGeometry()[i_node].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
            for (IndexType iDof = 0; iDof < TTensor; iDof++) {
                rDofData.u1(i_node, iDof) = value[iDof];
                rDofData.LagrangeMultipliers(i_node, iDof) = lm[iDof];
            }
        }
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DofData& rDofData,
    const array_1d<double, 3>& NormalMaster,
    const PointType& LocalPointDecomp,
    const PointType& LocalPointParent,
    GeometryPointType& GeometryDecomp,
    const bool DualLM
    )
{       
    /// SLAVE CONDITION ///
    /* SHAPE FUNCTIONS */
    GetGeometry().ShapeFunctionsValues( rVariables.NSlave, LocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = (DualLM == true) ? prod(rDofData.Ae, rVariables.NSlave) : rVariables.NSlave;
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.DetjSlave = GeometryDecomp.DeterminantOfJacobian( LocalPointDecomp );
    
    KRATOS_ERROR_IF(rVariables.DetjSlave < 0.0) << "WARNING:: CONDITION ID: " << this->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;
    
    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, NormalMaster, LocalPointParent);
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const array_1d<double, 3>& NormalMaster,
    const PointType& LocalPoint
    )
{    
    GeometryType& master_geometry = this->GetPairedGeometry();

    PointType projected_gp_global;
    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, LocalPoint );
    GeometricalProjectionUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, NormalMaster, -gp_normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    
    master_geometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;
    
    // SHAPE FUNCTIONS 
    master_geometry.ShapeFunctionsValues( rVariables.NMaster, projected_gp_local );         
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

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;  
    
    if (rResult.size() != MatrixSize) {
        rResult.resize( MatrixSize, false );
    }
    
    IndexType index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LM ] */
    // Master Nodes DoF Equation IDs
    GeometryType& current_master = this->GetPairedGeometry();
    
    if (TTensor == ScalarValue) {
        for ( IndexType i_master = 0; i_master < NumNodes; ++i_master ) {
            NodeType& master_node = current_master[i_master];
            rResult[index++] = master_node.GetDof( TEMPERATURE ).EquationId( );
        }
    } else {
        for ( IndexType i_master = 0; i_master < NumNodes; ++i_master ) {
            NodeType& master_node = current_master[i_master];
            rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3) {
                rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }
    }
    
    // Slave Nodes DoF Equation IDs
    if (TTensor == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rResult[index++] = slave_node.GetDof( TEMPERATURE ).EquationId( );
        }
    } else {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3) rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
        }
    }
    
    // Slave Nodes LM Equation IDs
    if (TTensor == ScalarValue) {
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

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim, TNumNodesElem, TTensor>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    if (rConditionalDofList.size() != MatrixSize) {
        rConditionalDofList.resize( MatrixSize );
    }
    
    IndexType index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LM ] */
    // Master Nodes DoF Equation IDs
    GeometryType& current_master = this->GetPairedGeometry();
    
    if (TTensor == ScalarValue) {
        for ( IndexType i_master = 0; i_master < NumNodes; ++i_master )  {
            NodeType& master_node = current_master[i_master];
            rConditionalDofList[index++] = master_node.pGetDof( TEMPERATURE );
        }
    } else {
        for ( IndexType i_master = 0; i_master < NumNodes; ++i_master ) {
            NodeType& master_node = current_master[i_master];
            rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3) rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
        }
    }
    
    // Slave Nodes DoF Equation IDs
    if (TTensor == ScalarValue) {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rConditionalDofList[index++] = slave_node.pGetDof( TEMPERATURE );
        }
    } else {
        for ( IndexType i_slave = 0; i_slave < NumNodes; ++i_slave ) {
            NodeType& slave_node = this->GetGeometry()[i_slave];
            rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
        }
    }
    
    // Slave Nodes LM Equation IDs
    if (TTensor == ScalarValue) {
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

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
        
    if ( rOutput.size() != integration_points.size() ) {
        rOutput.resize( integration_points.size() );
    }
    
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        rOutput[point_number] = 0.0;
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
                                                                                            
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
        
    if ( rOutput.size() != integration_points.size() ) {
        rOutput.resize( integration_points.size() );
    }
    
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        rOutput[point_number] = ZeroVector(3);
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<Vector>& rVariable, 
    std::vector<Vector>& rOutput, 
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
        
    if ( rOutput.size() != integration_points.size() ) {
        rOutput.resize( integration_points.size() );
    }
    
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        rOutput[point_number] = ZeroVector(3);
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodesElem, TensorValue TTensor>
int MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Check( const ProcessInfo& rCurrentProcessInfo )
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

template class MeshTyingMortarCondition<2, 3, ScalarValue>;   // 2DLine/Triangle for scalar variables
template class MeshTyingMortarCondition<2, 3, Vector2DValue>; // 2DLine/Triangle for components variables
template class MeshTyingMortarCondition<2, 4, ScalarValue>;   // 2DLine/Quadrilateral for scalar variables
template class MeshTyingMortarCondition<2, 4, Vector2DValue>; // 2DLine/Quadrilateral for scalar variables
template class MeshTyingMortarCondition<3, 4, ScalarValue>;   // 3D Triangle/Tetrahedron for scalar variables
template class MeshTyingMortarCondition<3, 4, Vector3DValue>; // 3D Triangle/Tetrahedron for components variables
template class MeshTyingMortarCondition<3, 8, ScalarValue>;   // 3D Quadrilateral/Hexahedra for scalar variables
template class MeshTyingMortarCondition<3, 8, Vector3DValue>; // 3D Quadrilateral/Hexahedra for components variables

} // Namespace Kratos
