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

#if !defined(KRATOS_MORTAR_CONTACT_2D_CONDITION_H_INCLUDED )
#define  KRATOS_MORTAR_CONTACT_2D_CONDITION_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include <vector>

// Project includes
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"


namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** \brief MortarContact2DCondition
 * This is a 2D contact condition which employes the mortar method with double lagrange multiplier 
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 */
class MortarContact2DCondition:
                public Condition {
public:
        ///@name Type Definitions
        ///@{
    
        typedef Point<3>                                  PointType;
        typedef Node<3>                                    NodeType;
        typedef Geometry<NodeType>                     GeometryType;
        
        /// Counted pointer of MortarContact2DCondition
        KRATOS_CLASS_POINTER_DEFINITION( MortarContact2DCondition );

protected:

    /**
    * Parameters to be used in the Condition as they are.
    */
    struct GeneralVariables
    {
    private:
        // Contact pair information
        GeometryType* pMasterElement;   // Points to the master contact segment only
        unsigned int mMasterElementIndex;
        std::vector< unsigned int > mThisActiveSlaveNodes;
        std::vector< unsigned int > mThisInactiveSlaveNodes;

    public:
        // Shape functions for contact pair
        Vector N_Master;
        Vector N_Slave;
        Vector Phi_LagrangeMultipliers;

        // Shape functions local derivatives for contact pair
        Matrix DN_De_Master;
        Matrix DN_De_Slave;
        Matrix DPhi_De_LagrangeMultipliers;

        // Gap function and its derivative variables
        Vector IntegrationPointNormalGap;
        Vector IntegrationPointNormalVector;

        // determinant of slave cell's jacobian
        double DetJSlave;

        /*
        * Jacobians in current configuration on all integration points of slave segment
        * Only those two variables contain info on all GP
        * other variables contain info only on the currently-calculated GP
        */
        GeometryType::JacobiansType j_Slave;

        /********************************************************/
        /******************** STRUCT METHODS ********************/
        /********************************************************/

        // Initializer method
        void Initialize( 
            const unsigned int& rLocalDimensionMaster = 1,  // Xi local coordinate
            const unsigned int& rLocalDimensionSlave  = 1,  // Xi local coordinate
            const unsigned int& rNumberOfMasterNodes  = 2,  // 2-node line segment
            const unsigned int& rNumberOfSlaveNodes   = 2,  // 2-node line segment
            const unsigned int& rDimension            = 2   // 2D physical space
            )
        {
            pMasterElement = NULL;
            mMasterElementIndex = 0;

            // Shape functions
            N_Master                     = ZeroVector( rNumberOfMasterNodes );
            N_Slave                      = ZeroVector( rNumberOfSlaveNodes  );
            Phi_LagrangeMultipliers      = ZeroVector( rNumberOfSlaveNodes  ); // Each slave node carries a lambda DOF

            // Shape functions local derivatives
            DN_De_Master                 = ZeroMatrix( rNumberOfMasterNodes, rLocalDimensionMaster );
            DN_De_Slave                  = ZeroMatrix( rNumberOfSlaveNodes ,  rLocalDimensionSlave );
            DPhi_De_LagrangeMultipliers  = ZeroMatrix( rNumberOfSlaveNodes ,  rLocalDimensionSlave );

            // Gap function
            IntegrationPointNormalGap    = ZeroVector( rDimension );
            IntegrationPointNormalVector = ZeroVector( rDimension );

            // Jacobian of slave
            DetJSlave = 0;

            // Jacobians on all integration points
            j_Slave.resize( 1, false );
            j_Slave[0] = ZeroMatrix( rDimension );
        }

        /* Setters and getters for the master element */
        void SetMasterElement( GeometryType& rSlaveElement ) 
        { 
            pMasterElement = &rSlaveElement;
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
      
        /* Setters and getters for the active set */
        void SetActiveSet( const std::vector< unsigned int >& rActiveNodes ) 
        {
            mThisActiveSlaveNodes = rActiveNodes; 
        }
        
        const std::vector< unsigned int >& GetActiveSet( ) 
        {
            return mThisActiveSlaveNodes;
        }
      
        /* Setters and getters for the active set */
        void SetInactiveSet( const std::vector< unsigned int >& rInactiveNodes ) 
        {
            mThisInactiveSlaveNodes = rInactiveNodes; 
        }
        
        const std::vector< unsigned int >& GetInactiveSet( )
        { 
            return mThisInactiveSlaveNodes;
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
    * Mortar condition matrices
    */
    struct MortarConditionMatrices
    {
    private:
        // Vector of projection matrices for all master elements in the contact condition
        std::vector<MatrixType> MortarProjectionMatricesTranspose;

    public:
       /*
        * Struct Member Variables
        */
        // Mortar condition matrices - D and M and their directional derivatives
        Matrix D;
        Matrix M;

       /*
        * Struct Methods
        */
        // Initializer method
        void Initialize( 
            const bool& rCalculateSlaveContributions  = false,
            const unsigned int& rLocalDimensionMaster = 1, // Xi local coordinate
            const unsigned int& rLocalDimensionSlave  = 1, // Xi local coordinate
            const unsigned int& rNumberOfMasterNodes  = 2, // 2-node line segment
            const unsigned int& rNumberOfSlaveNodes   = 2, // 2-node line segment
            const unsigned int& rDimension            = 2  // 2D physical space
            )
        {
            if( rCalculateSlaveContributions )
            {
                D = ZeroMatrix( rDimension * rNumberOfSlaveNodes, rDimension * rNumberOfSlaveNodes );
            }

            M = ZeroMatrix( rDimension * rNumberOfSlaveNodes, rDimension * rNumberOfMasterNodes );
        }
    };

public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MortarContact2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MortarContact2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    MortarContact2DCondition( MortarContact2DCondition const& rOther);

    /// Destructor.
    virtual ~MortarContact2DCondition();

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
    * Initialize System Matrices
    */
    virtual void InitializeSystemMatrices( 
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags
        );

    /**
    * Evaluation methods for Lagrange multipliers shape functions
    * and its local derivatives
    */
    double LagrangeMultiplierShapeFunctionValue( 
        const IndexType& rPointNumber,
        const IndexType& rShapeFunctionIndex 
        );

    const Matrix LagrangeMultiplierShapeFunctionLocalGradient( const IndexType& rPointNumber );

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
    virtual void InitializeGeneralVariables( 
        GeneralVariables& rVariables,
        const unsigned int& rMasterElementIndex 
        );

    /**
    * This function assigns the active and the inactive sets of the current condition
    */
    void InitializeActiveSet( GeneralVariables& rVariables );
                  
    /*
    * This fucntion loops over the slave nodes to determine the
    * active set of nodes and the inactive set of nodes
    * pointers to those nodes are stored in private arrays
    */                  
        void DetermineActiveAndInactiveSets(
        std::vector<unsigned int>& rActiveNodes,
        std::vector<unsigned int>& rInactiveNodes,
        const unsigned int& rCondIndex
        );

    /**
    * 
    */
    void CalculateNormalGapAtIntegrationPoint(
    GeneralVariables& rVariables,
    const unsigned int& rPointNumber
    );
    
    /**
    * Initialize the mortar condition matrices
    * D and M matrices and the directional derivatives matrices
    */
    virtual void InitializeConditionMatrices( 
        GeneralVariables& rVariables,
        double& rIntegrationWeight,
        const unsigned int& rPointNumber,
        const bool& rCalculateSlaveContributions
        );
                  
    /**
    * This function loops over all conditions and calculates the overall number of DOFs
    * total_dofs = SUM( master_u_dofs + slave_u_dofs + active_slave_lambda_dofs )
    */
    const unsigned int CalculateConditionSize( );

    /**
    * Loops over all master elements in the condition and
    * sum up their nodes - used to size the system matrices
    */
    const unsigned int CalculateTotalNumberOfMasterNodes( );
    
    /**
    * This method loops over the vector of active nodes stored in the container
    * correspoding to rPairIndex and counts the number of active nodes there
    */
    const unsigned int CalculateNumberOfActiveNodesInContactPair( const unsigned int& rPairIndex );

    /**
    * Calculate condition kinematics
    */
    void CalculateKinematics( 
        GeneralVariables& rVariables,
        const double& rPointNumber 
        );


    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /*
    * Calculate mortar matrices
    * @param 
    * 
    */
    void CalculateDAndM( 
        GeneralVariables& rVariables,
        double& rIntegrationWeight,
        const bool& rCalculateSlaveContributions 
        );

    /*******************************************************************************/
    /*******************************************************************************/

    /*
    * Calculation and addition of the matrices of the LHS of a contact pair
    */
    virtual void CalculateAndAddLHS( 
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        double& rIntegrationWeight 
        );

    /*
    * Assmbles the contact pair LHS block into the condition's LHS
    */
    void AssembleContactPairLHSToConditionSystem( 
        const unsigned int rMasterElementIndex,
        MatrixType& rPairLHS,
        MatrixType& rConditionLHS 
        );

    /*
    * Calculation and addition fo the vectors of the RHS of a contact pair
    */
    virtual void CalculateAndAddRHS( 
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        double& rIntegrationWeight 
        );
                  
    /*
    * 
    */
    void CalculateAndAddInternalForces(
        VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight
        );
    
    /*
    * Assmbles the contact pair RHS block into the condition's RHS
    */
    void AssembleContactPairRHSToConditionSystem( 
        const unsigned int rPairIndex,
        VectorType& rPairRHS,
        VectorType& rConditionRHS 
        );

    /***********************************************************************************/
    /**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
    /***********************************************************************************/

    /*
    * Calculates B_co = [ -M, D ], where D and M are the mortar condition matrices
    */
    void CalculateAndAddMortarContactOperator( 
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight 
    );
        
    
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

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
    * Calculate a double Variable
    */
    void CalculateOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
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

    IntegrationMethod mThisIntegrationMethod;

    // might be a good idea to move those three vectors to the contact container
    std::vector<Condition*> mThisMasterElements;

    MortarConditionMatrices mThisMortarConditionMatrices;

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
    MortarContact2DCondition( ) { };

    virtual void save( Serializer& rSerializer ) const;

    virtual void load( Serializer& rSerializer );
    ///@}

}; // Class MortarContact2DCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}
  // namespace Kratos.

#endif // KRATOS_MORTAR_CONTACT_2D_CONDITION_H_INCLUDED  defined 