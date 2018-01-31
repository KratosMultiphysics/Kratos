// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MESH_TYING_MORTAR_CONDITION_H_INCLUDED )
#define  KRATOS_MESH_TYING_MORTAR_CONDITION_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include <vector>

// Project includes
// #include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/condition.h"
#include "includes/kratos_flags.h"
#include "includes/mortar_classes.h"

/* Utilities */
#include "utilities/math_utils.h"

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
    
    typedef Point                                             PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                GeometryType;
    
    // Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsType;
    
    // Type definition of the components of an array_1d 
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;
    
///@}
///@name  Enum's
///@{
    
#if !defined(TENSOR_VALUE)
#define TENSOR_VALUE
    enum TensorValue {ScalarValue = 1, Vector2DValue = 2, Vector2DPScalarValue = 3, Vector3DValue = 3, Vector3DPScalarValue = 4 };
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

template< const unsigned int TDim, const unsigned int TNumNodesElem, TensorValue TTensor>
class MeshTyingMortarCondition: public Condition 
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Counted pointer of MeshTyingMortarCondition
    KRATOS_CLASS_POINTER_DEFINITION( MeshTyingMortarCondition );
    typedef Condition                                                                    BaseType;
    
    typedef typename BaseType::VectorType                                              VectorType;

    typedef typename BaseType::MatrixType                                              MatrixType;

    typedef typename BaseType::IndexType                                                IndexType;

    typedef typename BaseType::GeometryType::Pointer                          GeometryPointerType;

    typedef typename BaseType::NodesArrayType                                      NodesArrayType;

    typedef typename BaseType::PropertiesType::Pointer                      PropertiesPointerType;
    
    typedef typename std::vector<array_1d<PointType,TDim>>                 ConditionArrayListType;
    
    typedef Line2D2<Point>                                                            LineType;
    
    typedef Triangle3D3<Point>                                                    TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompositionType;

    static constexpr unsigned int NumNodes = (TNumNodesElem == 3 || (TDim == 2 && TNumNodesElem == 4)) ? 2 : TNumNodesElem == 4 ? 3 : 4;

//     static constexpr unsigned int MatrixSize = TTensor * (2* TNumNodesElem - NumNodes);
    static constexpr unsigned int MatrixSize = TTensor * (3 * NumNodes);
    
//     static constexpr unsigned int DimensionLocalElem = TTensor * TNumNodesElem;
    
    typedef MortarKinematicVariables<NumNodes>                                   GeneralVariables;
    
    typedef DualLagrangeMultiplierOperators<NumNodes>                                      AeData;
    
    typedef MortarOperator<NumNodes>                                      MortarConditionMatrices;
         
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MeshTyingMortarCondition(): Condition() 
    {
        mIntegrationOrder = 2; // Default value
    }
    
    // Constructor 1
    MeshTyingMortarCondition(IndexType NewId, GeometryType::Pointer pGeometry):Condition(NewId, pGeometry)
    {
        mIntegrationOrder = 2; // Default value
    }
    
    // Constructor 2
    MeshTyingMortarCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):Condition( NewId, pGeometry, pProperties )
    {
        mIntegrationOrder = 2; // Default value
    }

    ///Copy constructor
    MeshTyingMortarCondition( MeshTyingMortarCondition const& rOther){}

    /// Destructor.
    ~MeshTyingMortarCondition() override;

    /**
     * Flags related to the element computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
   /**
    * Called at the beginning of each solution step
    */
    void Initialize() override;

   /**
    * Called at the beginning of each solution step
    */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

   /**
    * Called at the beginning of each iteration
    */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
    * Called at the ending of each solution step
    */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;
    
   /**
    * Called at the end of each iteration
    */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;
    
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
        ) override;
    
    /**
    * Initialize Damping Matrix
    */
    void CalculateDampingMatrix( 
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;
    
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
        ) const override;
    
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
        PropertiesType::Pointer pProperties
        ) const override;
        
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
        ) override;

    /**
     * Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @return rConditionalDofList
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList( 
        DofsVectorType& rConditionalDofList,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
     * Get on rVariable a double Value
     */
    void GetValueOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo 
        ) override;
    
    /**
     * Get on rVariable a array_1d Value
     */
    void GetValueOnIntegrationPoints( 
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;
    
    /**
     * Get on rVariable a Vector Value
     */
    void GetValueOnIntegrationPoints( 
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Calculate a double Variable
     */
    void CalculateOnIntegrationPoints( 
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo 
        ) override;
    
    /**
     * Calculate a array_1d Variable
     */
    void CalculateOnIntegrationPoints( 
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector< array_1d<double, 3 > >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;
    
    /**
     * Calculate a Vector Variable
     */
    void CalculateOnIntegrationPoints( 
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

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
    struct DofData
    {
    public:
        
        // Auxiliar types
        typedef boost::numeric::ublas::bounded_matrix<double, NumNodes, TTensor>  Type1;
        typedef boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes> Type2;
        
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
        ~DofData()= default;
        
        // Initializer method 
        void Initialize(      
                const GeometryType& GeometryInput  // The geometry of the slave 
                )
        {
            SlaveGeometry  = GeometryInput;
            
            // The current Lagrange Multipliers
            u1 = ZeroMatrix(NumNodes, TTensor);
            u2 = ZeroMatrix(NumNodes, TTensor);
            LagrangeMultipliers = ZeroMatrix(NumNodes, TTensor);
        }
        
        // Initialize the Ae components
        void InitializeAeComponents()
        {
            // Ae
            Ae = ZeroMatrix(NumNodes, NumNodes);
        }
        
        // Updating the Master pair
        void UpdateMasterPair(
//         const GeometryType& GeometryInput,          // The geometry of the current master
            const Condition::Pointer& pCond          // The pointer of the current master
        )
        {
            const GeometryType GeometryInput =  pCond->GetGeometry();
            MasterGeometry = GeometryInput; // Updating the geometry
            
            /* DoF */
            if (TTensor == 1)
            {
                for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
                {
                    const double Value = MasterGeometry[iNode].FastGetSolutionStepValue(TEMPERATURE);
                    u2(iNode, 0) = Value;
                }
            }
            else
            {
                for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
                {
                    const array_1d<double, 3> Value = MasterGeometry[iNode].FastGetSolutionStepValue(DISPLACEMENT);
                    for (unsigned int iDof = 0; iDof < TTensor; iDof++)
                    {
                        u2(iNode, iDof) = Value[iDof];
                    }
                }
            }
        }

    };
    
    ///@}
    ///@name Protected member Variables
    ///@{

    /* Integration order */
    unsigned int mIntegrationOrder;                                      // The integration order considered
    
    /* Pair info */
    unsigned int mPairSize;                                              // The number of contact pairs
    std::vector<Condition::Pointer> mThisMasterConditions;               // Vector which contains the pointers to the master conditions

    std::vector<IntegrationPointsType> mIntegrationPointsVector;         // This vector contains all the integration nodes
    
//     /* The list of variable to compute */
//     array_1d<Variable< array_1d_component_type>, TDim> mTyingVarVector;  // Variable considered in the mesh tying
//     Variable<double> mTyingVarScalar;                                    // Variable considered in the mesh tying
    
//     /* Element info */
//     Element::Pointer mThisSlaveElement;                                  // The slave element from which derives everything
//     std::vector<Element::Pointer> mThisMasterElements;                   // Vector which contains the pointers to the master elements
    
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
        ) override;

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
        ) override;

    /**
     * This is called during the assembling process in order
     * to calculate the condition right hand side vector only
     * @param rRightHandSideVector: the condition right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

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
        ) override;

    /**
     * This is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide( 
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

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
        ) override;

    /**
     * Calculates the condition contribution
     */
    void CalculateConditionSystem( 
        LocalSystemComponents& rLocalSystem,
        const ProcessInfo& CurrentProcessInfo 
        );
    
    /**
     * Initialize Contact data
     */
    void InitializeDofData( 
        DofData& rDofData,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /**
     * Calculate Ae matrix
     */
    void CalculateAe( 
        DofData& rDofData,
        GeneralVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const unsigned int PairIndex
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
        const DofData rDofData,
        const array_1d<double, 3> MasterNormal,
        const double& rPointNumber,
        const IntegrationPointsType& IntegrationPointsSlave,
        const unsigned int PairIndex
        );

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /*
     * Calculation and addition of the matrices of the LHS of a contact pair
     */
    
    void CalculateAndAddLHS( 
        LocalSystemComponents& rLocalSystem,
        const boost::numeric::ublas::bounded_matrix<double, MatrixSize, MatrixSize>& LHS_contact_pair, 
        const unsigned int rPairIndex
        )
    {
        if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
        {
            /* COMPONENT-WISE LHS MATRIX */
            const std::vector<Variable<MatrixType> >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables( );

            for ( unsigned int i = 0; i < rLeftHandSideVariables.size( ); i++ )
            {
                MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrices( )[i];
                
                // Assemble in the correct position
                this->AssembleContactPairLHSToConditionSystem(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
            }
        }
        else 
        {   
            /* SINGLE LHS MATRIX */
            MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix( );      
            
            // Assemble in the correct position
            this->AssembleContactPairLHSToConditionSystem(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
        }
    }

    /*
     * Assembles the contact pair LHS block into the condition's LHS
     */
    
    void AssembleContactPairLHSToConditionSystem( 
        const boost::numeric::ublas::bounded_matrix<double, MatrixSize, MatrixSize>& rPairLHS,
        MatrixType& rConditionLHS,
        const unsigned int rPairIndex
        )
    {
        // Find location of the pair's master DOFs in ConditionRHS
        const unsigned int index_begin = rPairIndex * MatrixSize;
        const unsigned int index_end  = index_begin + MatrixSize;
        
        subrange( rConditionLHS, index_begin, index_end, index_begin, index_end) += rPairLHS;
    }

    /*
     * Calculates the local contibution of the LHS
     */
    
    template< unsigned int MatrixSize >
    boost::numeric::ublas::bounded_matrix<double, MatrixSize, MatrixSize> CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        DofData& rDofData,
//         const boost::numeric::ublas::bounded_matrix<double, DimensionLocalElem, DimensionLocalElem> LHS_SlaveElem_Contribution,
//         const Element::EquationIdVectorType& EquationIdSlaveElem,
        const unsigned int rMasterElementIndex,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /*
     * Calculation and addition fo the vectors of the RHS of a contact pair
     */

    void CalculateAndAddRHS( 
        LocalSystemComponents& rLocalSystem,
        const array_1d<double, MatrixSize>& RHS_contact_pair, 
        const unsigned int rPairIndex
        )
    {
        if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
        {
            /* COMPONENT-WISE RHS VECTOR */
            const std::vector<Variable<VectorType> >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables( );

            for ( unsigned int i = 0; i < rRightHandSideVariables.size( ); i++ )
            {
                VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVectors()[i];

                // Assemble
                this->AssembleContactPairRHSToConditionSystem( RHS_contact_pair, rRightHandSideVector, rPairIndex );
            }
        }
        else 
        {
            /* SINGLE RHS VECTOR */
            VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
            
            // Assemble
            this->AssembleContactPairRHSToConditionSystem( RHS_contact_pair, rRightHandSideVector, rPairIndex );
        }
    }
    
    /*
     * Assembles the contact pair RHS block into the condition's RHS
     */
    void AssembleContactPairRHSToConditionSystem( 
        const array_1d<double, MatrixSize>& rPairRHS,
        VectorType& rConditionRHS,
        const unsigned int rPairIndex
        )
    {
        // Find location of the pair's master DOFs in ConditionRHS
        const unsigned int index_begin = rPairIndex * MatrixSize;
        const unsigned int index_end  = index_begin + MatrixSize;
        
        subrange( rConditionRHS, index_begin, index_end) += rPairRHS;
    }
    
    /*
     * Calculates the local contibution of the LHS
     */
    template< unsigned int MatrixSize >
    array_1d<double, MatrixSize> CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        DofData& rDofData,
//         array_1d<double, DimensionLocalElem> RHS_SlaveElem_Contribution,
//         const Element::EquationIdVectorType& EquationIdSlaveElem,
        const unsigned int rMasterElementIndex,
        const ProcessInfo& rCurrentProcessInfo
        );
    
    /***********************************************************************************/
    /**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
    /***********************************************************************************/

    /*
     * Calculates the values of the shape functions for the master element
     */
    void MasterShapeFunctionValue(
        GeneralVariables& rVariables,
        const array_1d<double, 3> MasterNormal,
        const PointType& LocalPoint,
        const unsigned int PairIndex
        );
    
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/
    
    /**
     * It returns theintegration method considered
     */
    
    IntegrationMethod GetIntegrationMethod() override
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
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
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
