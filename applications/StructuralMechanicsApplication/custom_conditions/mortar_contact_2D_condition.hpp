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

namespace Kratos
{

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

    /// Counted pointer of MortarContact2DCondition
    KRATOS_CLASS_POINTER_DEFINITION( MortarContact2DCondition );

protected:

     /**
     * Flags related to the element computation
     */

     KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
     KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
     KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
     KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    /**
    * Parameters to be used in the Condition as they are.
    */

    struct GeneralVariables
    {
    private:
        // Contact pair information
        GeometryType* pMasterElement;   // Points to the slave contact segment only
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

        // Gap function and its derivative variables
        double IntegrationPointNormalGap;
        Vector IntegrationPointNormalVector;
        Matrix DeltaDiscreteGap;

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
            // fully populated BLOCKS
            DeltaDiscreteGap             = ZeroMatrix( rNumberOfSlaveNodes, rDimension * ( rNumberOfSlaveNodes + rNumberOfMasterNodes ) );
            IntegrationPointNormalGap    = 0.0;
            IntegrationPointNormalVector = ZeroVector( rDimension );

            // Jacobian of slave
            DetJSlave = 0;

            // Jacobians on all integration points
            j_Slave.resize( 1, false );
            j_Slave[0] = ZeroMatrix( rDimension );
        }

        /*
        * Setters and getters for the master element info
        */
        void SetMasterElement( GeometryType& rSlaveElement )
        {
            pMasterElement = &rSlaveElement;
        }

        GeometryType& GetMasterElement( )
        {
            return *pMasterElement;
        }

        void SetMasterElementIndex( const unsigned int& index )
        {
            mMasterElementIndex = index;
        }

        const unsigned int& GetMasterElementIndex( )
        {
            return mMasterElementIndex;
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

        // Directional derivatives
        Matrix DeltaD;
        Matrix DeltaM;
        Matrix DeltaJSlave;
        Matrix DeltaPhiLagrangeMultipliers;
        Matrix DeltaDiscreteGap;

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
                D                       = ZeroMatrix( rDimension * rNumberOfSlaveNodes, rDimension * rNumberOfSlaveNodes );
                DeltaD                  = ZeroMatrix( rDimension * rNumberOfSlaveNodes, rDimension * rNumberOfSlaveNodes );
                DeltaJSlave             = ZeroMatrix( rNumberOfSlaveNodes, rDimension * rNumberOfSlaveNodes );
            }

            M                           = ZeroMatrix( rDimension * rNumberOfSlaveNodes, rDimension * rNumberOfMasterNodes );
            DeltaM                      = ZeroMatrix( rDimension * rNumberOfSlaveNodes, rDimension * rNumberOfMasterNodes );

            DeltaPhiLagrangeMultipliers = ZeroMatrix( rNumberOfSlaveNodes,  rDimension * rNumberOfSlaveNodes );
            DeltaDiscreteGap            = ZeroMatrix( rNumberOfSlaveNodes, rDimension * ( rNumberOfSlaveNodes + rNumberOfMasterNodes ) );
        }

       /*
        * Accessors for MortarProjectionMatrices
        */
        void AppendMortarProjectionMatrixTranspose( const MatrixType& P_transpose )
        {
                MortarProjectionMatricesTranspose.push_back( P_transpose );
        }

        const MatrixType& GetMortarProjectionMatrixTranspose( const unsigned int& rMatrixIndex )
        {
                return MortarProjectionMatricesTranspose[rMatrixIndex];
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

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    //************* STARTING - ENDING  METHODS

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

    /*
    * This fucntion loops over the slave nodes to determine the
    * active set of nodes and the inactive set of nodes
    * pointers to those nodes are stored in private arrays
    */
    void DetermineActiveNodes( );

    /*
    * Evaluation methods for Lagrange multipliers shape functions
    * and its local derivatives
    */
    double LagrangeMultiplierShapeFunctionValue( 
        const IndexType& rPointNumber,
        const IndexType& rShapeFunctionIndex 
        );

    const Matrix LagrangeMultiplierShapeFunctionLocalGradient( const IndexType& rPointNumber );

    /*
    * Loops over all master elements in the condition and
    * sum up their nodes - used to size the system matrices
    */
    const unsigned int CalculateTotalNumberOfMasterNodes( );

    /**
    * Creates a new element pointer
    * @param NewId: the ID of the new element
    * @param ThisNodes: the nodes of the new element
    * @param pProperties: the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Condition::Pointer Create( 
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties 
        ) const;


   /******************************************************************/
    /*********************** COMPUTING  METHODS ***********************/
    /******************************************************************/

    /**
    * this is called during the assembling process in order
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
    * this function provides a more general interface to the condition.
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
    * this is called during the assembling process in order
    * to calculate the condition right hand side vector only
    * @param rRightHandSideVector: the condition right hand side vector
    * @param rCurrentProcessInfo: the current process info instance
    */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
    * this function provides a more general interface to the condition.
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
    * this is called during the assembling process in order
    * to calculate the condition left hand side matrix only
    * @param rLeftHandSideMatrix: the condition left hand side matrix
    * @param rCurrentProcessInfo: the current process info instance
    */
    void CalculateLeftHandSide( 
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
    * this function provides a more general interface to the condition.
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

    /*
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
    * Initialize General Variables
    */
    virtual void InitializeConditionMatrices( 
        GeneralVariables& rVariables,
        double& rIntegrationWeight,
        const unsigned int& rPointNumber,
        const bool& rCalculateSlaveContributions
        );

    /*
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

    /*
    * Calculate mortar matrices directional derivatives
    */
    void CalculateDeltaDAndDeltaM( 
        GeneralVariables& rVariables,
        double& rIntegrationWeight,
        const bool& rCalculateSlaveContributions 
        );

    /*
    * Calculate slave element's Jacobian directional derivative
    */
    void CalculateDeltaJSlave( 
        GeneralVariables& rVariables,
        double& rIntegrationWeight,
        const bool& rCalculateSlaveContributions,
        const IndexType& rPointNumber
        );

    /*
    * Calculate Lagrange multiplier shape functions directional derivative
    */
    void CalculateDeltaPhi( 
        GeneralVariables& rVariables,
        double& rIntegrationWeight,
        const bool& rCalculateSlaveContributions
        );

    /*
     * Calculate the normal gap between a given integration point on the
     * slave element and its projected position on the master element
     */
    void CalculateNormalGapAtIntegrationPoint(
        GeneralVariables& rVariables,
        const IndexType& rPointNumber );


    /*
     * Calculates the linearization of the discrete gap function
     * g_h = normal_vector * ( x_gauss_point - x_gauss_point_projection )
     */
    void CalculateDeltaDiscreteGap(
        GeneralVariables& rVariables,
        const IndexType& rPointNumber );


    void CalculateDeltaIntegrationPointNormal(
        GeneralVariables& rVariables,
        const IndexType& rPointNumber,
        const unsigned int& rIndexSlave,
        MatrixType& rDeltaIntPtNormal );


    void CalculateDeltaNodalNormal(
        GeneralVariables& rVariables,
        const IndexType& rPointNumber,
        MatrixType& rDeltaNodalNormal );


    void CalculateDeltaElementalNormal(
        GeneralVariables& rVariables,
        const IndexType& rPointNumber,
        MatrixType& rDeltaElementalNormal );

    /*
     * Calcualte Mortar Projection Operator tranpose - P' = ( D_{AA}^{-1} * M_{A} )'
     * and appends it to the vector of projection operators
     * This is the LHS contribution related to all Lagrange multipliers DOFs in the contact surfaces
     * This is only applicable for DUAL LAGRANGE MULTIPLIERS
     */
    void CalculateAndAppendMortarProjectionMatrixTranspose( GeneralVariables& rVariables );

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
    * Assmbles the contact pair RHS block into the condition's RHS
    */
    void AssembleContactPairRHSToConditionSystem( 
        const unsigned int rMasterElementIndex,
        VectorType& rPairRHS,
        VectorType& rConditionRHS 
        );

    /***********************************************************************************/
    /**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
    /***********************************************************************************/

    /*
    * Calculate contact stiffness - contribution of the mortar condition to the system stiffness matrix
    * This is the LHS contribution related to all displacement DOFs in the contact surfaces
    */
    void CalculateAndAddKco( 
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight 
    );

    /*
    * Calculate weighted gap derivatives
    * This is the contributions of the derivatives of the weighted gap functions to the system stiffness matrix
    * These contributions are related to the displacemenet DOFs of the nodes of the contact surfaces
    */
    void CalculateAndAddGapDerivatives( 
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight 
        );


    /**************************************************************************/
    /********** AUXILLIARY METHODS FOR CONDITION CONTRIBUTION TO RHS **********/
    /**************************************************************************/

    /*
    * Add the contact forces from the previous solution step to the RHS vector
    */
    void CalculateAndAddContactForces( 
        VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        double& rIntegrationWeight 
        );
    
    /*
    * Add the nodal gaps to the RHS vector
    */
    void CalculateAndAddGapsVector( 
        VectorType& rRightHandSideVector,
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

    /*
    * This is an interface function to access the P matrices stored in the matrices struct
    * the builder will require these matrices to assemble the stiffness matrix of the non-contact elements
    * and the residuals vector in the RHS
    */
    void GetMortarProjectionMatricesTranspose( std::vector< MatrixType >& rMortarProjectionMatricesTranspose );

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

    std::vector<GeometryType::Pointer> mThisMasterElements;
    std::vector< unsigned int > mThisInactiveSlaveNodes;
    std::vector< unsigned int > mThisActiveSlaveNodes;

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