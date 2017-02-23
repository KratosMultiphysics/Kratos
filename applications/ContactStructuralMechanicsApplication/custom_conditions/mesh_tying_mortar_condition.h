// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_MESH_TYING_MORTAR_CONDITION_H_INCLUDED )
#define  KRATOS_MESH_TYING_MORTAR_CONDITION_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include <vector>

// Project includes
#include <contact_structural_mechanics_application_variables.h>
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"

/* Utilities */
#include "custom_utilities/logging_settings.hpp"

namespace Kratos 
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    typedef Point<3>                                             PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                GeometryType;
    
    ///Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsType;
    
///@}
///@name  Enum's
///@{
    
    
#if !defined(TENSOR_VALUE)
#define TENSOR_VALUE
    enum TensorValue {ScalarValue = 1, VectorValue = 3};
#endif
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{
    
/** \brief MeshTyingMortarCondition
 * This is a mesh tying condition which employes the mortar method with dual lagrange multiplier 
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 */

template< unsigned int TDim, unsigned int TNumNodes, TensorValue TTensor>
class MeshTyingMortarCondition: public Condition 
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Counted pointer of MeshTyingMortarCondition
    KRATOS_CLASS_POINTER_DEFINITION( MeshTyingMortarCondition );

    typedef Condition                                                  BaseType;
    
    typedef typename BaseType::VectorType                            VectorType;

    typedef typename BaseType::MatrixType                            MatrixType;

    typedef typename BaseType::IndexType                              IndexType;

    typedef typename BaseType::GeometryType::Pointer        GeometryPointerType;

    typedef typename BaseType::NodesArrayType                    NodesArrayType;

    typedef typename BaseType::PropertiesType::Pointer    PropertiesPointerType;

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
        VectorType N_Master;
        VectorType N_Slave;
        VectorType Phi_LagrangeMultipliers;

        // Determinant of slave cell's jacobian
        double DetJSlave;
        
        /*
         * Jacobians in current configuration on all integration points of slave segment
         * Only those two variables contain info on all GP
         * other variables contain info only on the currently-calculated GP
         */
        MatrixType j_Slave;
        
        /********************************************************/
        /******************** STRUCT METHODS ********************/
        /********************************************************/

        // Initializer method 
        void Initialize()
        {
            pMasterElement = NULL;
            mMasterElementIndex = 0;

            // Shape functions
            N_Master                = ZeroVector(TNumNodes);
            N_Slave                 = ZeroVector(TNumNodes);
            Phi_LagrangeMultipliers = ZeroVector(TNumNodes);
            
            // Jacobian of slave
            DetJSlave = 0.0;
           
            // Jacobians on all integration points
            j_Slave = ZeroMatrix(TDim, TDim - 1);
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
//             LOG_VECTOR_PRETTY( N_Slave );
            KRATOS_WATCH( N_Master );
//             LOG_VECTOR_PRETTY( N_Master );
            KRATOS_WATCH( Phi_LagrangeMultipliers );
//             LOG_VECTOR_PRETTY( Phi_LagrangeMultipliers );
            KRATOS_WATCH( j_Slave );
//             LOG_MATRIX_PRETTY( j_Slave );
            KRATOS_WATCH( DetJSlave );
        }
    };
         
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MeshTyingMortarCondition(): Condition() {}
    
    // Constructor 1
    MeshTyingMortarCondition(IndexType NewId, GeometryType::Pointer pGeometry):Condition(NewId, pGeometry){}
    
    // Constructor 2
    MeshTyingMortarCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):Condition( NewId, pGeometry, pProperties )
    {
        InitializeIntegrationMethod(); 
    }

    ///Copy constructor
    MeshTyingMortarCondition( MeshTyingMortarCondition const& rOther){}

    /// Destructor.
    virtual ~MeshTyingMortarCondition();

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
    * Called at the ending of each solution step
    */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
    
   /**
    * Called at the end of each iteration
    */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod();
    
    IntegrationMethod GetIntegrationMethod(const unsigned int integration_order, );
    
    /**
    * Initialize System Matrices
    */
    template<unsigned int TMatrixSize>
    void InitializeSystemMatrices( 
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
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

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

    /** 
     * This data will be used to compute teh derivatives
     */ 
    struct DerivativeData
    {
    public:
        
        // Auxiliar types
        typedef bounded_matrix<double, TNumNodes, TTensor>   Type1;
        typedef bounded_matrix<double, TNumNodes, TNumNodes> Type2;
        
        // Master and element geometries
        GeometryType SlaveGeometry;
        GeometryType MasterGeometry;
        
        // The current Lagrange Multipliers
        Type1 LagrangeMultipliers;
        
        // DoF
        Type1 u1;
        Type1 u2;
        
        // Ae
        Type2 Ae;
        
        // Default destructor
        ~DerivativeData(){}
        
        // Initializer method 
        void Initialize(      
                const GeometryType& GeometryInput  // The geometry of the slave 
                )
        {
            SlaveGeometry  = GeometryInput;
            
            // The current Lagrange Multipliers
            LagrangeMultipliers = ZeroMatrix(TNumNodes, TTensor);
            
            // DoF of the slave       
            if (TTensor == ScalarValue)
            {
                Variable<double> VarType = this->GetValue(SCALAR_VARIABLE);
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    const double var  = SlaveGeometry[iNode].FastGetSolutionStepValue(VarType);

                    u1(iNode, 0) = var;
                }
            }
            else if (TTensor == VectorValue)
            {
                Variable<array_1d<double, 3>> VarType = this->GetValue(COMPONENTS_VARIABLE);
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    const array_1d<double, 3> var  = SlaveGeometry[iNode].FastGetSolutionStepValue(VarType);

                    for (unsigned int iDof = 0; iDof < TDim; iDof++)
                    {
                        u1(iNode, iDof) = var[iDof];
                    }
                }
            }
        }
        
        // Initialize the Ae components
        void InitializeAeComponents()
        {
            // Ae
            Ae = ZeroMatrix(TNumNodes, TNumNodes);
        }
    
        // Updating the Master pair
        void UpdateMasterPair(
    //         const GeometryType& GeometryInput,          // The geometry of the current master
            const Condition::Pointer& pCond          // The pointer of the current master
        )
        {
            const GeometryType GeometryInput =  pCond->GetGeometry();
            MasterGeometry = GeometryInput; // Updating the geometry
            
            // DoF of the master
            if (TTensor == ScalarValue)
            {
                Variable<double> VarType = pCond->GetValue(SCALAR_VARIABLE);
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    const double var  = MasterGeometry[iNode].FastGetSolutionStepValue(VarType);

                    u2(iNode, 0) = var;
                }
            }
            else if (TTensor == VectorValue)
            {
                Variable<array_1d<double, 3>> VarType = pCond->GetValue(COMPONENTS_VARIABLE);
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    const array_1d<double, 3> var  = MasterGeometry[iNode].FastGetSolutionStepValue(VarType);

                    for (unsigned int iDof = 0; iDof < TDim; iDof++)
                    {
                        u2(iNode, iDof) = var[iDof];
                    }
                }
            }
        }
    };
    
    /** 
     * This data will be used to compute the Ae matrix
     */ 
    struct AeData
    {
    public:
        
        // Auxiliar types
        typedef bounded_matrix<double, TNumNodes, TNumNodes> Type3;
        
        // Matrices
        Type3 Me;
        Type3 De;
        
        // Default destructor
        ~AeData(){}
        
        // Initialize the Ae components
        void Initialize()
        {
            // Matrices
            Me = ZeroMatrix(TNumNodes, TNumNodes);
            De = ZeroMatrix(TNumNodes, TNumNodes);
        }
    };
    
    
   /*
    * Mortar condition matrices
    */
    struct MortarConditionMatrices
    {
    private:

    public:
       /*
        * Struct Member Variables
        */
       
        // Mortar condition matrices - DOperator and MOperator
        bounded_matrix<double, TNumNodes, TNumNodes> DOperator;
        bounded_matrix<double, TNumNodes, TNumNodes> MOperator;

       /*
        * Struct Methods
        */
        // Initializer method
        void Initialize()
        {
            // We initialize the D and M operators
            DOperator = ZeroMatrix(TNumNodes, TNumNodes);
            MOperator = ZeroMatrix(TNumNodes, TNumNodes);
        }
        
        void print() 
        {
            KRATOS_WATCH(DOperator);
//             LOG_MATRIX_PRETTY(DOperator);
            KRATOS_WATCH(MOperator);
//             LOG_MATRIX_PRETTY(MOperator);
        }
    };

    ///@}
    ///@name Protected member Variables
    ///@{

    /* Pair info */
    unsigned int mPairSize;                                          // The number of contact pairs
    std::vector<Condition::Pointer> mThisMasterConditions;           // Vector which contains the pointers to the master conditions

    /* Integration points */
    std::vector<IntegrationPointsType> mThisSlaveIntegrationPoints;  // The arrays containing the integration points of the slave side
    /* Element info */
    Element::Pointer mThisSlaveElement;                              // The slave element from which derives everything
    std::vector<Element::Pointer> mThisMasterElements;               // Vector which contains the pointers to the master elements
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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
    template<unsigned int TMatrixSize>
    void CalculateConditionSystem( 
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
    void InitializeDerivativeData( 
        DerivativeData& rDerivativeData,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Calculate Ae matrix
     */
    void CalculateAe( 
        DerivativeData& rDerivativeData,
        GeneralVariables& rVariables,
//         const GeometryType::IntegrationPointsArrayType& integration_points,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Initialize General Variables
     */
    void UpdateDerivativeData( 
        DerivativeData& rDerivativeData,
        const unsigned int& rMasterElementIndex
        );
    
    /**
     * This function loops over all conditions and calculates the overall number of DOFs
     * total_dofs = SUM( master_u_dofs + 2 * slave_u_dofs) 
     */
    template<unsigned int TMatrixSize>
    const unsigned int CalculateConditionSize( );
    
    /**
     * Calculate condition kinematics
     */
    bool CalculateKinematics( 
        GeneralVariables& rVariables,
        const DerivativeData rDerivativeData,
        const double& rPointNumber,
        const GeometryType::IntegrationPointsArrayType& integration_points
        );

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /*
     * Calculation and addition of the matrices of the LHS of a contact pair
     */
    template<unsigned int TMatrixSize>
    void CalculateAndAddLHS( 
        LocalSystemComponents& rLocalSystem,
        const bounded_matrix<double, TMatrixSize, TMatrixSize>& LHS_contact_pair, 
        const unsigned int rPairIndex
        );

    /*
     * Assembles the contact pair LHS block into the condition's LHS
     */
    template<unsigned int TMatrixSize>
    void AssembleContactPairLHSToConditionSystem( 
        const bounded_matrix<double, TMatrixSize, TMatrixSize>& rPairLHS,
        MatrixType& rConditionLHS,
        const unsigned int rPairIndex
        );

    /*
     * Calculates the local contibution of the LHS
     */
    template<unsigned int TMatrixSize>
    bounded_matrix<double, TMatrixSize, TMatrixSize> CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
        );
    
    /*
     * Calculation and addition fo the vectors of the RHS of a contact pair
     */
    template<unsigned int TMatrixSize>
    void CalculateAndAddRHS( 
        LocalSystemComponents& rLocalSystem,
        const array_1d<double, TMatrixSize>& RHS_contact_pair, 
        const unsigned int rPairIndex
        );
    
    /*
     * Assembles the contact pair RHS block into the condition's RHS
     */
    template<unsigned int TMatrixSize>
    void AssembleContactPairRHSToConditionSystem( 
        const array_1d<double, TMatrixSize>& rPairRHS,
        VectorType& rConditionRHS,
        const unsigned int rPairIndex
        );
    
    /*
     * Calculates the local contibution of the LHS
     */
    template<unsigned int TMatrixSize>
    array_1d<double, TMatrixSize> CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
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

    void CalculateMortarOperators(
        MortarConditionMatrices& rThisMortarConditionMatrices,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight
    );
    /*
     * Calculates the matrix De
     */
    bounded_matrix<double, TNumNodes, TNumNodes> ComputeDe(        
        const array_1d<double, TNumNodes> N1, 
        const double detJ 
        );
    
    /*
     * Calculates the matrix De
     */
    bounded_matrix<double, TNumNodes, TNumNodes> ComputeMe(        
        const array_1d<double, TNumNodes> N1, 
        const double detJ 
        );
    
    /*
     * Calculates the matrix Ae
     */
    void CalculateAe(
        DerivativeData& rDerivativeData,
        AeData& rAeData
        );
    
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/
    
    /*
     * Computes the selective integration method
     */
    void ComputeSelectiveIntegrationMethod(const unsigned int rPairIndex);
    
    /*
     * Initializes the integration method
     */
    void InitializeIntegrationMethod();
    
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

    // Serialization 
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

    ///@}

}; // Class MeshTyingMortarCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MESH_TYING_MORTAR_CONDITION_H_INCLUDED  defined 
