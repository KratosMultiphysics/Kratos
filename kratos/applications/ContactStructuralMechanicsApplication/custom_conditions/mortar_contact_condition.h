// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Mohamed Khalil
//                   Vicente Mataix Ferrándiz
// 

#if !defined(KRATOS_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include <vector>

// Project includes
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"

namespace Kratos 
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    typedef Point<3>                                  PointType;
    typedef Node<3>                                    NodeType;
    typedef Geometry<NodeType>                     GeometryType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod   IntegrationMethod;
    
///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** \brief ContactData 
 * This data will be used to compute de external condition files
 */ 
struct ContactData
{
public:
    
    // Master and element geometries
    GeometryType SlaveGeometry;
    GeometryType MasterGeometry;
    
    // Gap function and its derivative variables
    Vector Gaps;
    
    // The current Lagrange Multipliers
    Matrix LagrangeMultipliers;
    
    // The normals of the nodes
    Matrix NormalsMaster;
    Matrix NormalsSlave;
    
    Matrix Tangent1Slave;
    Matrix Tangent2Slave;
    
    // Delta time
    double Dt;
    
    // Augmentation parameter
    double epsilon_normal;
    double epsilon_tangent;
    
    // Initializer method 
    void Initialize(      
            const GeometryType& GeometryInput,          // The geometry of the slave 
            const unsigned int& rNumberOfSlaveNodes,    // Number of nodes of the slave
            const unsigned int& rDimension              // 3D/2D physical space
            )
    {
        SlaveGeometry  = GeometryInput;
            
        // Gap function and its derivative variables
        Gaps = ZeroVector(rNumberOfSlaveNodes);
        
        // The current Lagrange Multipliers
        LagrangeMultipliers = ZeroMatrix(rNumberOfSlaveNodes, rDimension);
        
        // The normals of the nodes
        NormalsSlave  = ZeroMatrix(rNumberOfSlaveNodes, rDimension);
        
        Tangent1Slave = ZeroMatrix(rNumberOfSlaveNodes, rDimension);
        if (rDimension == 3)
        {
            Tangent2Slave = ZeroMatrix(rNumberOfSlaveNodes, rDimension);
        }
        
        // Delta time 
        Dt = 0.0;
        
        // Augmentation parameter
        epsilon_normal  = 0.0;
        epsilon_tangent = 0.0;
    }
    
    // Updating the Master pair
    void UpdateMasterPair(
//         const GeometryType& GeometryInput,          // The geometry of the current master
        const Condition::Pointer& pCond,          // The pointer of the current master
        const unsigned int& rNumberOfMasterNodes,   // Number of nodes of the master
        const unsigned int& rDimension             // 3D/2D physical space
    )
    {
        const GeometryType GeometryInput =  pCond->GetGeometry();
        MasterGeometry = GeometryInput; // Updating the geometry
        
        NormalsMaster = ZeroMatrix(rNumberOfMasterNodes, rDimension);
        
        for (unsigned int iNode = 0; iNode < rNumberOfMasterNodes; iNode++)
        {
            array_1d<double,3> normal = pCond->GetValue(NORMAL); // TODO: To consider an interpolation it is necessary to smooth the surface
//             array_1d<double,3> normal = MasterGeometry[iNode].GetValue(NORMAL);

            for (unsigned int iDof = 0; iDof < rDimension; iDof++)
            {
                NormalsMaster(iNode, iDof) = normal[iDof]; 
            }
        }
    }
};
    
    
/** \brief MortarContactCondition
 * This is a contact condition which employes the mortar method with double lagrange multiplier 
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 */
class MortarContactCondition:
                public Condition 
{
public:
        ///@name Type Definitions
        ///@{
        
        /// Counted pointer of MortarContactCondition
        KRATOS_CLASS_POINTER_DEFINITION( MortarContactCondition );

protected:

    /**
     * Parameters to be used in the Condition as they are.
     */
    struct GeneralVariables
    {
    private:
        // Contact pair information
        GeometryType* pMasterElement;   // Pointer to the master contact segment only
        unsigned int mMasterElementIndex;

    public:
        // Shape functions for contact pair
        Vector N_Master;
        Vector N_Slave;
        Vector Phi_LagrangeMultipliers;

        // Shape functions local derivatives for contact pair
        Matrix DN_De_Master;
        Matrix DN_De_Slave;
        Matrix DPhi_De_LagrangeMultipliers;

        // Determinant of slave cell's jacobian
        double DetJSlave;
        
        /*
         * Jacobians in current configuration on all integration points of slave segment
         * Only those two variables contain info on all GP
         * other variables contain info only on the currently-calculated GP
         */
        Matrix j_Slave;
        
        // Friction coefficient
        double mu;
        
        /********************************************************/
        /******************** STRUCT METHODS ********************/
        /********************************************************/

        // Initializer method 
        void Initialize(
            const unsigned int& rLocalDimensionMaster,  // Xi local coordinate
            const unsigned int& rLocalDimensionSlave,   // Xi local coordinate
            const unsigned int& rNumberOfMasterNodes,   // 2-node line segment
            const unsigned int& rNumberOfSlaveNodes,    // 2-node line segment
            const unsigned int& rDimension,             // 3D/2D physical space
            const unsigned int& rIntegrationPointNumber // Number of integration points to be considered
            )
        {
            pMasterElement = NULL;
            mMasterElementIndex = 0;

            // Shape functions
            N_Master.resize( rNumberOfMasterNodes, false );
            N_Slave.resize( rNumberOfSlaveNodes, false );
            Phi_LagrangeMultipliers.resize( rNumberOfSlaveNodes, false );
            N_Master                     = ZeroVector( rNumberOfMasterNodes );
            N_Slave                      = ZeroVector( rNumberOfSlaveNodes  );
            Phi_LagrangeMultipliers      = ZeroVector( rNumberOfSlaveNodes  ); // Each slave node carries a lambda DOF

            // Shape functions local derivatives
            DN_De_Master.resize( rNumberOfMasterNodes, rLocalDimensionMaster, false);
            DN_De_Slave.resize( rNumberOfSlaveNodes,  rLocalDimensionSlave, false);
            DPhi_De_LagrangeMultipliers.resize( rNumberOfSlaveNodes ,  rLocalDimensionSlave, false);
            DN_De_Master                 = ZeroMatrix( rNumberOfMasterNodes, rLocalDimensionMaster );
            DN_De_Slave                  = ZeroMatrix( rNumberOfSlaveNodes ,  rLocalDimensionSlave );
            DPhi_De_LagrangeMultipliers  = ZeroMatrix( rNumberOfSlaveNodes ,  rLocalDimensionSlave );

            // Jacobian of slave
            DetJSlave = 0.0;
           
            // Jacobians on all integration points
            j_Slave.resize( rNumberOfSlaveNodes, 1, false );
            
            // Friction coefficient
            mu = 0.0;
        }

        /* Setters and getters for the master element */
        void SetMasterElement( GeometryType& rMasterElement ) 
        { 
            pMasterElement = &rMasterElement;
        }
        
        GeometryType& GetMasterElement( ) 
        { 
            return *pMasterElement;
        }

        /* Setters and getters for the master element index */
        void SetMasterElementIndex( const unsigned int& index ) 
        { 
            mMasterElementIndex = index; 
        }
        
        const unsigned int& GetMasterElementIndex( ) 
        { 
            return mMasterElementIndex; 
        }
        
        void print( )
        {
            KRATOS_WATCH( N_Slave );
            KRATOS_WATCH( N_Master );
            KRATOS_WATCH( Phi_LagrangeMultipliers );
//             KRATOS_WATCH( DN_De_Slave );
//             KRATOS_WATCH( DN_De_Master );
//             KRATOS_WATCH( DPhi_De_LagrangeMultipliers );
            
//             KRATOS_WATCH( ColocationWeightCoeff );
//             KRATOS_WATCH( pMasterElement->Center( ).Coordinates( ) );
//             KRATOS_WATCH( mMasterElementIndex );
            
//             KRATOS_WATCH( IntegrationPointNormalGap );
//             KRATOS_WATCH( IntegrationPointNormalVector );
//             KRATOS_WATCH( j_Master );
            KRATOS_WATCH( j_Slave );
            KRATOS_WATCH( DetJSlave );
//             KRATOS_WATCH( SegmentProportion );
        }
    };

   /**
    * This struct is used in the component wise calculation only
    * is defined here and is used to declare a member variable in the component wise condition
    * private pointers can only be accessed by means of set and get functions
    * this allows to set and not copy the local system variables
    */
    struct LocalSystemComponents
    {
    private:
            //for calculation local system with compacted LHS and RHS
            MatrixType *mpLeftHandSideMatrix;
            VectorType *mpRightHandSideVector;

            //for calculation local system with LHS and RHS components
            std::vector<MatrixType> *mpLeftHandSideMatrices;
            std::vector<VectorType> *mpRightHandSideVectors;
            
            //LHS variable components
            const std::vector< Variable< MatrixType > > *mpLeftHandSideVariables;

            //RHS variable components
            const std::vector< Variable< VectorType > > *mpRightHandSideVariables;

    public:
            // Calculation flags
            Flags  CalculationFlags;

           /**
            * Sets the value of a specified pointer variable
            */
            void SetLeftHandSideMatrix( MatrixType& rLeftHandSideMatrix ) { mpLeftHandSideMatrix = &rLeftHandSideMatrix; };
            void SetLeftHandSideMatrices( std::vector<MatrixType>& rLeftHandSideMatrices ) { mpLeftHandSideMatrices = &rLeftHandSideMatrices; };
            void SetLeftHandSideVariables(const std::vector< Variable< MatrixType > >& rLeftHandSideVariables ) { mpLeftHandSideVariables = &rLeftHandSideVariables; };

            void SetRightHandSideVector( VectorType& rRightHandSideVector ) { mpRightHandSideVector = &rRightHandSideVector; };
            void SetRightHandSideVectors( std::vector<VectorType>& rRightHandSideVectors ) { mpRightHandSideVectors = &rRightHandSideVectors; };
            void SetRightHandSideVariables(const std::vector< Variable< VectorType > >& rRightHandSideVariables ) { mpRightHandSideVariables = &rRightHandSideVariables; };

           /**
            * Returns the value of a specified pointer variable
            */
            MatrixType& GetLeftHandSideMatrix() { return *mpLeftHandSideMatrix; };
            std::vector<MatrixType>& GetLeftHandSideMatrices() { return *mpLeftHandSideMatrices; };
            const std::vector< Variable< MatrixType > >& GetLeftHandSideVariables() { return *mpLeftHandSideVariables; };

            VectorType& GetRightHandSideVector() { return *mpRightHandSideVector; };
            std::vector<VectorType>& GetRightHandSideVectors() { return *mpRightHandSideVectors; };
            const std::vector< Variable< VectorType > >& GetRightHandSideVariables() { return *mpRightHandSideVariables; };
    };

    /*
     * Colocation Integration
     */
    struct ColocationIntegration
    {
    private:
        GeometryType::IntegrationPointsArrayType mIntegrationPoints;
        
    public:
        void Initialize( const unsigned int num_integration_points_per_local_dim, const unsigned int local_dimension_slave , const unsigned int number_of_slave_nodes  )
        {   
            if (local_dimension_slave == 1)
            {
                const unsigned int num_integration_points = num_integration_points_per_local_dim;
                mIntegrationPoints.resize( num_integration_points, false );
            
                const double elem_local_length = 2.0;
                const double colocation_weight = elem_local_length / num_integration_points;
                
                for ( unsigned int i_col_point = 0; i_col_point < num_integration_points; ++i_col_point )
                {
                    const double xi = 0.5 * colocation_weight * ( 2 * i_col_point + 1 ) - 1;
                    mIntegrationPoints[i_col_point] = IntegrationPoint<2>( xi, colocation_weight );
                }

            }
            else
            {
                if (number_of_slave_nodes == 3)
                {
                    unsigned int num_integration_points = 0;
                    for ( unsigned int i_xi = 0; i_xi < num_integration_points_per_local_dim; ++i_xi )
                    {
                        for ( unsigned int j_eta = i_xi; j_eta < num_integration_points_per_local_dim; ++j_eta )
                        {
                            num_integration_points += 1;
                        }
                    }
                     
                    mIntegrationPoints.clear( );
                    mIntegrationPoints.resize( num_integration_points, false );
                    const double elem_local_length = 0.5;   // both xi and eta local coords span between -1 and 1
                    const double colocation_weight = elem_local_length / num_integration_points;
                    
                    double xi = 0, eta = 0;
                    unsigned int i_col_pt = 0;
                    for ( unsigned int i_xi = 0; i_xi < num_integration_points_per_local_dim; ++i_xi )
                    {
                        xi = 1.0/(3.0 * num_integration_points_per_local_dim) + (i_xi - 1) * (1.0/(1.0 * num_integration_points_per_local_dim));
                        for ( unsigned int j_eta = i_xi; j_eta < num_integration_points_per_local_dim; ++j_eta )
                        {
                            eta =  1.0/(3.0 * num_integration_points_per_local_dim) + (j_eta - 1) * (1.0/(1.0 * num_integration_points_per_local_dim));
                            mIntegrationPoints[i_col_pt] = IntegrationPoint<2>( xi, eta, colocation_weight ); 
                            i_col_pt += 1.0;
                        }
                    }
                }
                else if (number_of_slave_nodes == 4)
                {
                    const unsigned int num_integration_points = num_integration_points_per_local_dim * num_integration_points_per_local_dim;
                    mIntegrationPoints.clear( );
                    mIntegrationPoints.resize( num_integration_points, false );
                    const double elem_local_length = 2.0;   // both xi and eta local coords span between -1 and 1
                    const double colocation_weight = elem_local_length / num_integration_points_per_local_dim;
                    
                    double xi = 0, eta = 0;
                    unsigned int i_col_pt = 0;
                    for ( unsigned int i_xi = 0; i_xi < num_integration_points_per_local_dim; ++i_xi )
                    {
                        xi = 0.5 * colocation_weight * ( 2 * i_xi + 1 ) - 1;
                        for ( unsigned int j_eta = 0; j_eta < num_integration_points_per_local_dim; ++j_eta )
                        {
                            i_col_pt = i_xi * num_integration_points_per_local_dim + j_eta;
                            eta = 0.5 * colocation_weight * ( 2 * j_eta + 1 ) - 1;
                            mIntegrationPoints[i_col_pt] = IntegrationPoint<2>( xi, eta, colocation_weight*colocation_weight ); // w_xi * w_eta
                        }
                    }
                }
                else
                {
                     KRATOS_THROW_ERROR( std::logic_error,  " Collocation not defined. number_of_slave_nodes: ", number_of_slave_nodes );
                }
            }
        }
        
        const GeometryType::IntegrationPointsArrayType IntegrationPoints( )
        {
            return mIntegrationPoints;
        }
        
        void print( )
        {
            std::cout << "mIntegrationPoints : " << mIntegrationPoints.size( ) << std::endl;
            for ( unsigned int i_vec = 0; i_vec < mIntegrationPoints.size( ); ++i_vec )
            {
                KRATOS_WATCH( mIntegrationPoints[i_vec] );
            }
        }
        
    };
    
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MortarContactCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MortarContactCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    MortarContactCondition( MortarContactCondition const& rOther);

    /// Destructor.
    virtual ~MortarContactCondition();

    /**
     * Flags related to the element computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

   /**
    * Called at the beginning of each solution step
    */
    void Initialize();

   /**
    * Called at the beginning of each solution step
    */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

   /**
    * Called at the beginning of each iteration
    */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

   /**
    * Called at the end of each iteration
    */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod();
    
    /**
    * Initialize System Matrices
    */
    virtual void InitializeSystemMatrices( 
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags
        );

    /**
    * Initialize Mass Matrix
    */
    void CalculateMassMatrix( 
    MatrixType& rMassMatrix, 
    ProcessInfo& rCurrentProcessInfo
    );
    
    /**
    * Initialize Damping Matrix
    */
    void CalculateDampingMatrix( 
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    );
    
    /**
     * Evaluation methods for Lagrange multipliers shape functions
     * and its local derivatives
     */
    const Vector LagrangeMultiplierShapeFunctionValue(
        const double xi_local,
        const double eta_local
        );
    
    /**
     * Evaluation methods for Lagrange multipliers shape functions
     * and its local derivatives
     */
    const Matrix LagrangeMultiplierShapeFunctionLocalGradient( 
        const double xi_local,
        const double eta_local
        );

    /**
     * Creates a new element pointer from an arry of nodes
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create( 
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties 
        ) const;
    
    /**
     * Creates a new element pointer from an existing geometry
     * @param NewId: the ID of the new element
     * @param pGeom: the  geometry taken to create the condition
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const;

    /******************************************************************/
    /*********************** COMPUTING  METHODS ***********************/
    /******************************************************************/

    /**
     * This is called during the assembling process in order
     * to calculate all condition contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rRightHandSideVector: the condition right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem( 
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * This function provides a more general interface to the condition.
     * it is designed so that rLHSvariables and rRHSvariables are passed TO the condition
     * thus telling what is the desired output
     * @param rLeftHandSideMatrices: container with the output left hand side matrices
     * @param rLHSVariables: paramter describing the expected LHSs
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateLocalSystem( 
        std::vector< MatrixType >& rLeftHandSideMatrices,
        const std::vector< Variable< MatrixType > >& rLHSVariables,
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * This is called during the assembling process in order
     * to calculate the condition right hand side vector only
     * @param rRightHandSideVector: the condition right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * This function provides a more general interface to the condition.
     * it is designed so that rRHSvariables are passed TO the condition
     * thus telling what is the desired output
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateRightHandSide(
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * This is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide( 
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * This function provides a more general interface to the condition.
     * it is designed so that rRHSvariables are passed TO the condition
     * thus telling what is the desired output
     * @param rRightHandSideVectors: container for the desired LHS output
     * @param rRHSVariables: parameter describing the expected LHSs
     */
    void CalculateLeftHandSide( 
        std::vector< MatrixType >& rLeftHandSideMatrices,
        const std::vector< Variable< MatrixType > >& rLHSVariables,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * Calculates the condition contribution
     */
    virtual void CalculateConditionSystem( 
        LocalSystemComponents& rLocalSystem,
        const ProcessInfo& CurrentProcessInfo 
        );

    /**
     * Initialize General Variables
     */
    void InitializeGeneralVariables( 
        GeneralVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const unsigned int& rMasterElementIndex
        );
    
    /**
     * Initialize Contact data
     */
    void InitializeContactData( 
        ContactData& rContactData,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Initialize General Variables
     */
    void UpdateContactData( 
        ContactData& rContactData,
        const unsigned int& rMasterElementIndex
        );
    
    /**
     * This function loops over all conditions and calculates the overall number of DOFs
     * total_dofs = SUM( master_u_dofs + 2 * slave_u_dofs) 
     */
    const unsigned int CalculateConditionSize( );
    
    /**
     * Calculate condition kinematics
     */
    bool CalculateKinematics( 
        GeneralVariables& rVariables,
        const double& rPointNumber,
        const unsigned int& rPairIndex,
        const GeometryType::IntegrationPointsArrayType& integration_points
        );

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /*
     * Calculation and addition of the matrices of the LHS of a contact pair
     */
    virtual void CalculateAndAddLHS( 
        LocalSystemComponents& rLocalSystem,
        MatrixType& LHS_contact_pair, 
        const unsigned int rPairIndex,
        const GeometryType& current_master_element
        );

    /*
     * Assembles the contact pair LHS block into the condition's LHS
     */
    void AssembleContactPairLHSToConditionSystem( 
        MatrixType& rPairLHS,
        MatrixType& rConditionLHS,
        const unsigned int rPairIndex
        );

    /*
     * 
     */
    void CalculateLocalLHS(
        MatrixType& rPairLHS,
        GeneralVariables& rVariables,
        const ContactData& rContactData,
        const double& rIntegrationWeight,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        );
    
    /*
     * Calculation and addition fo the vectors of the RHS of a contact pair
     */
    virtual void CalculateAndAddRHS( 
        LocalSystemComponents& rLocalSystem,
        VectorType& RHS_contact_pair, 
        const unsigned int rPairIndex,
        const GeometryType& current_master_element
        );
    
    /*
     * Assembles the contact pair RHS block into the condition's RHS
     */
    void AssembleContactPairRHSToConditionSystem( 
        VectorType& rPairRHS,
        VectorType& rConditionRHS,
        const unsigned int rPairIndex
        );
    
    /*
     * 
     */
    void CalculateLocalRHS(
        VectorType& rPairRHS,
        GeneralVariables& rVariables,
        const ContactData& rContactData,
        const double& rIntegrationWeight,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        );
    
    /***********************************************************************************/
    /**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
    /***********************************************************************************/

    /*
     * Calculates the values of the shape functions for the master element
     */
    bool MasterShapeFunctionValue(
        GeneralVariables& rVariables,
        const PointType& local_point 
    );
        
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    double AugmentedNormalLM(
        const GeneralVariables& rVariables,
        const ContactData& rContactData,
        const double& integration_point_gap,
        const unsigned int& dimension
    );
    
    double AugmentedTangentLM(
        const GeneralVariables& rVariables,
        const ContactData& rContactData,
        const GeometryType& current_master_element, 
        double& integration_point_slip,
        const unsigned int& dimension
    );
    
    void InitializeIntegrationMethod( 
        const unsigned int dimension,
        const unsigned int local_dimension_slave,
        const unsigned int number_of_slave_nodes
    );
    
    void GetNodalDeltaMovements( 
        Vector& rValues,
        const int& rNode 
        );

    Vector& GetCurrentValue( 
        const Variable<array_1d<double,3> >& rVariable,
        Vector& rValue,
        const unsigned int& rNode 
        );

    Vector& GetPreviousValue( 
        const Variable<array_1d<double,3> >& rVariable,
        Vector& rValue,
        const unsigned int& rNode 
        );

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @return rResult: The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector( 
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        );

    /**
     * Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @return rConditionalDofList
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList( 
        DofsVectorType& rConditionalDofList,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * Get on rVariable a double Value
     */
    void GetValueOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo 
        );
    
    /**
     * Get on rVariable a array_1d Value
     */
    void GetValueOnIntegrationPoints( 
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Get on rVariable a Vector Value
     */
    void GetValueOnIntegrationPoints( 
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * Calculate a double Variable
     */
    void CalculateOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo 
        );
    
    /**
     * Calculate a array_1d Variable
     */
    void CalculateOnIntegrationPoints( 
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector< array_1d<double, 3 > >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Calculate a Vector Variable
     */
    void CalculateOnIntegrationPoints( 
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        );

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    IntegrationMethod mThisIntegrationMethod;              // Integration order of the element
    std::vector<Condition::Pointer> mThisMasterElements;   // Vector which contains the pointers to the master elements

    bool mUseManualColocationIntegration;                  // Use the manual collocation integration
    ColocationIntegration mColocationIntegration;          // The manual collocation integration
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    MortarContactCondition( ) { };

    virtual void save( Serializer& rSerializer ) const;

    virtual void load( Serializer& rSerializer );
    ///@}

}; // Class MortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined 
