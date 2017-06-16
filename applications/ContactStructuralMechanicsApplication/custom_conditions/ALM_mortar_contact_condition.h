// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_ALM_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_ALM_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"
#include "boost/smart_ptr.hpp"
#include <vector>

// Project includes
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"

/* Utilities */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/logging_settings.hpp"

/* Geometries */
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"

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
    typedef Geometry<PointType>               GeometryPointType;
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
  
/** \brief DerivativeData
 * This data will be used to compute the derivatives
 */
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
class DerivativeData
{
public:
    
    // Auxiliar types
//     typedef int zero[0]; // NOTE: Problems in Windows
    typedef array_1d<double, TNumNodes> Type1;
    typedef bounded_matrix<double, TNumNodes, TDim> Type2;
    typedef bounded_matrix<double, TNumNodes, TNumNodes> Type3;
    
    // Auxiliar sizes
    static const unsigned int Size1 =     (TNumNodes * TDim);
    static const unsigned int Size2 = 2 * (TNumNodes * TDim);
    
    // Master and element geometries
    GeometryType SlaveGeometry;
    GeometryType MasterGeometry;
    
    // The ALM parameters
    double PenaltyParameter;
    double ScaleFactor;
    
    typename std::conditional< TFrictional,double,int >::type TangentFactor;
    
    // The normals of the nodes
    Type2 NormalMaster;
    Type2 NormalSlave;
    
    // Displacements and velocities
    Type2 X1;
    Type2 X2;
    Type2 u1;
    Type2 u2;
    
    typename std::conditional< TFrictional,Type2,int >::type u1old;
    typename std::conditional< TFrictional,Type2,int >::type u2old;
    
    // Derivatives    
    array_1d<double, Size1> DeltaJSlave;
    array_1d<Type1,  Size1> DeltaPhi;
    array_1d<Type1,  Size2> DeltaN1;
    array_1d<Type1,  Size2> DeltaN2;
    array_1d<Type2,  Size1> DeltaNormalSlave;
    array_1d<Type2,  Size1> DeltaNormalMaster;
    
    // Ae
    Type3 Ae;
    
    // Derivatives Ae
    array_1d<Type3, Size1> DeltaAe;
    
    // Default destructor
    ~DerivativeData(){}
    
    /**
     * Initializer method 
     * @param  GeometryInput: The geometry of the slave 
     */
    
    virtual void Initialize(
        const GeometryType& GeometryInput,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        SlaveGeometry  = GeometryInput;
        
        // The normals of the nodes
        NormalSlave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry,  NORMAL);
        
        // Displacements and velocities of the slave       
        u1 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 0);
        X1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(SlaveGeometry, false);
        
        // We get the ALM variables
        PenaltyParameter = rCurrentProcessInfo[PENALTY_PARAMETER];
        ScaleFactor = rCurrentProcessInfo[SCALE_FACTOR];
        
        // Derivatives 
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaPhi[i] = ZeroVector(TNumNodes);
            DeltaN1[i] = ZeroVector(TNumNodes);
            DeltaN1[i + TNumNodes * TDim] = ZeroVector(TNumNodes);
            DeltaN2[i] = ZeroVector(TNumNodes);
            DeltaN2[i + TNumNodes * TDim] = ZeroVector(TNumNodes);
            DeltaNormalSlave[i]      = ZeroMatrix(TNumNodes, TDim);
        }
        
        #if (TFrictional == true)
            TangentFactor = rCurrentProcessInfo[TANGENT_FACTOR];
            u1old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1);
        #endif
    }
    
    /**
     * Initialize the DeltaAe components
     */
    
    void InitializeDeltaAeComponents()
    {
        // Ae
        Ae = ZeroMatrix(TNumNodes, TNumNodes);
        
        // Derivatives Ae
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaAe[i] = ZeroMatrix(TNumNodes, TNumNodes);
        }
    }

    /**
     * Updating the Master pair
     * @param  pCond: The pointer of the current master
     */
    
    virtual void UpdateMasterPair(const Condition::Pointer& pCond)
    {
        MasterGeometry =  pCond->GetGeometry();
        
        NormalMaster = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry,  NORMAL);
        
        // Displacements, coordinates and normals of the master
        u2 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 0);
        X2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(MasterGeometry, false);

        // Derivative of master's normal
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaNormalMaster[i] = ZeroMatrix(TNumNodes, TDim);
        }
        
        #if (TFrictional == true)
            u2old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 1);
        #endif
    }
};  // Class DerivativeData
    
/** \brief AugmentedLagrangianMethodMortarContactCondition
 * This is a contact condition which employes the mortar method with dual lagrange multiplier 
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 */
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
class AugmentedLagrangianMethodMortarContactCondition: public Condition 
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Counted pointer of AugmentedLagrangianMethodMortarContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( AugmentedLagrangianMethodMortarContactCondition );

    typedef Condition                                                                    BaseType;
    
    typedef typename BaseType::VectorType                                              VectorType;

    typedef typename BaseType::MatrixType                                              MatrixType;

    typedef typename BaseType::IndexType                                                IndexType;

    typedef typename BaseType::GeometryType::Pointer                          GeometryPointerType;

    typedef typename BaseType::NodesArrayType                                      NodesArrayType;

    typedef typename BaseType::PropertiesType::Pointer                      PropertiesPointerType;
    
    typedef typename std::vector<array_1d<PointType,TDim>>                 ConditionArrayListType;
    
    typedef Line2D2<Point<3>>                                                            LineType;
    
    typedef Triangle3D3<Point<3>>                                                    TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompositionType;
    
    typedef DerivativeData<TDim, TNumNodes, TFrictional>                       DerivativeDataType;
    
    static constexpr unsigned int MatrixSize = TFrictional == true ? TDim * (TNumNodes + TNumNodes + TNumNodes) : TDim * (TNumNodes + TNumNodes) + TNumNodes;
    
    /**
     * Parameters to be used in the Condition as they are.
     */
    struct GeneralVariables
    {
    private:
        // Contact pair information
        Condition::Pointer pMasterElement;   // Pointer to the master contact segment only

    public:
        // Shape functions for contact pair
        VectorType NMaster;
        VectorType NSlave;
        VectorType PhiLagrangeMultipliers;

        // Shape functions local derivatives for contact pair
        MatrixType DNDeMaster;
        MatrixType DNDeSlave;
//         MatrixType DPhi_De_LagrangeMultipliers;

        // Determinant of slave cell's jacobian
        double DetjSlave;
        
        /*
         * Jacobians in current configuration on all integration points of slave segment
         * Only those two variables contain info on all GP
         * other variables contain info only on the currently-calculated GP
         */
        MatrixType jSlave;
        
        /********************************************************/
        /******************** STRUCT METHODS ********************/
        /********************************************************/

        // Initializer method 
        void Initialize()
        {
            pMasterElement = nullptr;

            // Shape functions
            NMaster                = ZeroVector(TNumNodes);
            NSlave                 = ZeroVector(TNumNodes);
            PhiLagrangeMultipliers = ZeroVector(TNumNodes);

            // Shape functions local derivatives
            DNDeMaster                 = ZeroMatrix(TNumNodes, TDim - 1);
            DNDeSlave                  = ZeroMatrix(TNumNodes, TDim - 1);
//             DPhi_De_LagrangeMultipliers  = ZeroMatrix(TNumNodes, TDim - 1);
            
            // Jacobian of slave
            DetjSlave = 0.0;
           
            // Jacobians on all integration points
            jSlave = ZeroMatrix(TDim, TDim - 1);
        }

        /* Setters and getters for the master element */
        void SetMasterElement( Condition::Pointer rMasterElement ) 
        { 
            pMasterElement = rMasterElement;
        }
        
        Condition::Pointer GetMasterElement( ) 
        { 
            return pMasterElement;
        }
        
        void print( )
        {
            KRATOS_WATCH( NSlave );
            KRATOS_WATCH( NMaster );
            KRATOS_WATCH( PhiLagrangeMultipliers );
            KRATOS_WATCH( jSlave );
            KRATOS_WATCH( DetjSlave );
        }
    };
         
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AugmentedLagrangianMethodMortarContactCondition(): Condition() 
    {
        mIntegrationOrder = 2; // Default value
    }
    
    // Constructor 1
    AugmentedLagrangianMethodMortarContactCondition(IndexType NewId, GeometryType::Pointer pGeometry):Condition(NewId, pGeometry)
    {
        mIntegrationOrder = 2; // Default value
    }
    
    // Constructor 2
    AugmentedLagrangianMethodMortarContactCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):Condition( NewId, pGeometry, pProperties )
    {
        mIntegrationOrder = 2; // Default value
    }

    ///Copy constructor
    AugmentedLagrangianMethodMortarContactCondition( AugmentedLagrangianMethodMortarContactCondition const& rOther){}

    /// Destructor.
    virtual ~AugmentedLagrangianMethodMortarContactCondition();

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
    * Initialize System Matrices
    */
    
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
    
    virtual Condition::Pointer Create( 
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
    
    virtual Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const;
        
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @return rResult: The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo: the current process info instance
     */
    
    virtual void EquationIdVector( 
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        );

    /**
     * Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @return rConditionalDofList
     * @param rCurrentProcessInfo: the current process info instance
     */
    
    virtual void GetDofList( 
        DofsVectorType& rConditionalDofList,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * Get on rVariable a double Value
     * @param rVariable: Internal values
     * @param rCurrentProcessInfo: The current process information
     * @return rValues: The values of interest (doubles)
     */
    
    void GetValueOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo 
        );
    
    /**
     * Get on rVariable a array_1d Value
     * @param rVariable: Internal values
     * @param rCurrentProcessInfo: The current process information
     * @return rValues: The values of interest (array_1d)
     */
    
    void GetValueOnIntegrationPoints( 
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Get on rVariable a Vector Value
     * @param rVariable: Internal values
     * @param rCurrentProcessInfo: The current process information
     * @return rValues: The values of interest (vector)
     */
    
    void GetValueOnIntegrationPoints( 
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * Calculate a double Variable
     * @param rVariable: Internal values
     * @param rCurrentProcessInfo: The current process information
     * @return rOutput: The values of interest (doubles)
     */
    
    void CalculateOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo 
        );
    
    /**
     * Calculate a array_1d Variable
     * @param rVariable: Internal values
     * @param rCurrentProcessInfo: The current process information
     * @return rOutput: The values of interest (array_1d)
     */
    
    void CalculateOnIntegrationPoints( 
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector< array_1d<double, 3 > >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Calculate a Vector Variable
     * @param rVariable: Internal values
     * @param rCurrentProcessInfo: The current process information
     * @return rOutput: The values of interest (vector)
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
     * This data will be used to compute the Ae and DeltaAe matrices
     */ 
    struct AeData
    {
    public:
        
        // Auxiliar types
        typedef bounded_matrix<double, TNumNodes, TNumNodes> Type3;
        
        // Auxiliar sizes
        static const unsigned int Size1 = (TNumNodes * TDim);
        
        // Matrices
        Type3 Me;
        Type3 De;

        // Derivatives matrices
        array_1d<Type3, Size1> DeltaMe;
        array_1d<Type3, Size1> DeltaDe;
        
        // Default destructor
        ~AeData(){}
        
        // Initialize the DeltaAe components
        void Initialize()
        {
            // Matrices
            Me = ZeroMatrix(TNumNodes, TNumNodes);
            De = ZeroMatrix(TNumNodes, TNumNodes);
            // Derivatives matrices
            for (unsigned int i = 0; i < TNumNodes * TDim; i++)
            {
                DeltaMe[i] = ZeroMatrix(TNumNodes, TNumNodes);
                DeltaDe[i] = ZeroMatrix(TNumNodes, TNumNodes);
            }
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
       
        static const unsigned int Size2 = 2 * (TNumNodes * TDim);
       
        // Mortar condition matrices - DOperator and MOperator
        bounded_matrix<double, TNumNodes, TNumNodes> DOperator;
        bounded_matrix<double, TNumNodes, TNumNodes> MOperator;
        
        // D and M directional derivatives
        array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2> DeltaDOperator;
        array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2> DeltaMOperator;

       /*
        * Struct Methods
        */
        // Initializer method
        void Initialize()
        {
            // We initialize the D and M operators
            DOperator = ZeroMatrix(TNumNodes, TNumNodes);
            MOperator = ZeroMatrix(TNumNodes, TNumNodes);
            
            for (unsigned int i = 0; i < TNumNodes * TDim; i++)
            {
                DeltaDOperator[i] = ZeroMatrix(TNumNodes, TNumNodes);
                DeltaDOperator[i + TNumNodes * TDim] = ZeroMatrix(TNumNodes, TNumNodes);
                DeltaMOperator[i] = ZeroMatrix(TNumNodes, TNumNodes);
                DeltaMOperator[i + TNumNodes * TDim] = ZeroMatrix(TNumNodes, TNumNodes);
            }
        }
        
        void print() 
        {
            KRATOS_WATCH(DOperator);
            KRATOS_WATCH(MOperator);
            
//             for (unsigned int i = 0; i < TNumNodes * TDim; i++)
//             {
//                 KRATOS_WATCH(DeltaDOperator[i]);
//                 KRATOS_WATCH(DeltaMOperator[i]);
//             }
        }
    };

    ///@}
    ///@name Protected member Variables
    ///@{

    IntegrationMethod mThisIntegrationMethod;            // Integration order of the element
    unsigned int mPairSize;                              // The number of contact pairs
    std::vector<Condition::Pointer> mThisMasterElements; // Vector which contains the pointers to the master elements
   
    unsigned int mIntegrationOrder;                      // The integration order to consider
    
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
     * Calculate Ae and DeltaAe matrices
     */
    
    bool CalculateAeAndDeltaAe( 
        DerivativeDataType& rDerivativeData,
        GeneralVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * This function loops over all conditions and calculates the overall number of DOFs
     * total_dofs = SUM( master_u_dofs + 2 * slave_u_dofs) 
     */
    
    const unsigned int CalculateConditionSize( );
    
    /**
     * Calculate condition kinematics
     */
    
    void CalculateKinematics( 
        GeneralVariables& rVariables,
        const DerivativeDataType rDerivativeData,
        const array_1d<double, 3> MasterNormal,
        const PointType& LocalPointDecomp,
        const PointType& LocalPointParent,
        GeometryPointType& GeometryDecomp,
        const bool DualLM = true
        );

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /**
     * Calculation and addition of the matrices of the LHS of a contact pair
     */

    void CalculateAndAddLHS( 
        LocalSystemComponents& rLocalSystem,
        const bounded_matrix<double, MatrixSize, MatrixSize>& LHS_contact_pair, 
        const unsigned int rPairIndex
        );

    /**
     * Assembles the contact pair LHS block into the condition's LHS
     */
    
    void AssembleContactPairLHSToConditionSystem( 
        const bounded_matrix<double, MatrixSize, MatrixSize>& rPairLHS,
        MatrixType& rConditionLHS,
        const unsigned int rPairIndex
        );

    /**
     * Calculates the local contibution of the LHS
     */
    
    virtual bounded_matrix<double, MatrixSize, MatrixSize> CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        );
    
    /**
     * Calculation and addition fo the vectors of the RHS of a contact pair
     */
    
    void CalculateAndAddRHS( 
        LocalSystemComponents& rLocalSystem,
        const array_1d<double, MatrixSize>& RHS_contact_pair, 
        const unsigned int rPairIndex
        );
    
    /**
     * Assembles the contact pair RHS block into the condition's RHS
     */
    
    void AssembleContactPairRHSToConditionSystem( 
        const array_1d<double, MatrixSize>& rPairRHS,
        VectorType& rConditionRHS,
        const unsigned int rPairIndex
        );
    
    /**
     * Calculates the local contibution of the LHS
     */
    
    virtual array_1d<double, MatrixSize> CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        );
    
    /********************************************************************************/
    /*************** METHODS TO CALCULATE MORTAR CONDITION DERIVATIVES **************/
    /********************************************************************************/
    
    void CalculateDeltaDetjSlave(
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData
        );
    
    bounded_matrix<double, TDim, TDim> LocalDeltaNormal(
        const GeometryType& CondGeometry,
        const unsigned int node_index
        );

    /**
     * Calculates the increment of the normal in the slave condition
     */
    
    void CalculateDeltaNormalSlave(DerivativeDataType& rDerivativeData);
    
    /**
     * Calculates the increment of the normal and in the master condition
     */
    
    void CalculateDeltaNormalMaster(DerivativeDataType& rDerivativeData);
    
    /**
     * Calculates the increment of the shape functions and the gap
     */
    
    void CalculateDeltaN(
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData
        );
    
    /**
     * Calculates the increment of Phi
     */
    
    void CalculateDeltaPhi(
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData
        );
    
    /***********************************************************************************/
    /**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
    /***********************************************************************************/

    /**
     * Calculates the values of the shape functions for the master element
     */
    
    void MasterShapeFunctionValue(
        GeneralVariables& rVariables,
        const array_1d<double, 3> MasterNormal,
        const PointType& local_point 
    );
    
    /**
     * Calculates the componets necessaries to compute the mortar operators and its derivatives
     * @param rVariables: Internal values
     * @param rDerivativeData: The derivative data
     * @param rIntegrationWeight: The integration weight
     * @return rThisMortarConditionMatrices: The mortar operators
     */
    
    virtual void CalculateMortarOperators(
        MortarConditionMatrices& rThisMortarConditionMatrices,
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        const double& rIntegrationWeight
        );

    /**
     * Calculates the componets necessaries to compute the mortar operators
     * @param rVariables: Internal values
     * @param rIntegrationWeight: The integration weight
     * @return rThisMortarConditionMatrices: The mortar operators
     */
        
    virtual void CalculateMortarOperators(
        MortarConditionMatrices& rThisMortarConditionMatrices,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight
        );
    
    /**
     * Calculates the componets necessaries to compute the derivatives of Phi
     * @param rVariables: Internal values
     * @param rDerivativeData: The derivative data
     * @param rIntegrationWeight: The integration weight
     * @return AeData: The Ae matrix and derivatives
     */
    
    virtual void CalculateDeltaAeComponents(
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        AeData& rAeData,
        const double& rIntegrationWeight
        );
    
    /**
     * Calculates the matrix De
     */
    
    bounded_matrix<double, TNumNodes, TNumNodes> ComputeDe(        
        const array_1d<double, TNumNodes> N1, 
        const double detJ 
        );
    
    /**
     * Calculates the matrix DeltaAe
     */
    
    bool CalculateDeltaAe(
        DerivativeDataType& rDerivativeData,
        AeData& rAeData
        );
    
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/
    
    /**
     * Returns a value depending of the active/inactive set
     */
    
    virtual unsigned int GetActiveInactiveValue(GeometryType& CurrentGeometry) const
    {
        KRATOS_ERROR << "You are calling to the base class method GetActiveInactiveValue, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
        
        return 0;
    }

    /**
     * It returns theintegration method considered
     */
    
    IntegrationMethod GetIntegrationMethod()
    {
        if (mIntegrationOrder == 1)
        {
            return GeometryData::GI_GAUSS_1;
        }
        else if (mIntegrationOrder == 2)
        {
            return GeometryData::GI_GAUSS_2;
        }
        else if (mIntegrationOrder == 3)
        {
            return GeometryData::GI_GAUSS_3;
        }
        else if (mIntegrationOrder == 4)
        {
            return GeometryData::GI_GAUSS_4;
        }
        else if (mIntegrationOrder == 5)
        {
            return GeometryData::GI_GAUSS_5;
        }
        else
        {
            return GeometryData::GI_GAUSS_2;
        }
    }
    
    /**
     * Returns a matrix with the increment of displacements, that can be used for compute the Jacobian "perturbed"
     */
    
    Matrix CalculateDeltaPosition();
    
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

}; // Class AugmentedLagrangianMethodMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_ALM_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined 
