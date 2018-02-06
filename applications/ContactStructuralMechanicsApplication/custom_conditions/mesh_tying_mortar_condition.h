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

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "custom_conditions/paired_condition.h"
#include "includes/mortar_classes.h"

/* Utilities */
#include "utilities/math_utils.h"
#include "utilities/exact_mortar_segmentation_utility.h"

namespace Kratos 
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    typedef Point                                                PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                GeometryType;
    
    // Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsType;
    
    // Type definition of the components of an array_1d 
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;
    
///@}
///@name  Enum's
///@{
    
    enum TensorValue {ScalarValue = 1, Vector2DValue = 2, Vector2DPScalarValue = 3, Vector3DValue = 3, Vector3DPScalarValue = 4 };
    
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
class MeshTyingMortarCondition: public PairedCondition 
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Counted pointer of MeshTyingMortarCondition
    KRATOS_CLASS_POINTER_DEFINITION( MeshTyingMortarCondition );
    typedef PairedCondition                                                              BaseType;
    
    typedef typename BaseType::VectorType                                              VectorType;

    typedef typename BaseType::MatrixType                                              MatrixType;

    typedef typename BaseType::IndexType                                                IndexType;

    typedef typename BaseType::GeometryType::Pointer                          GeometryPointerType;

    typedef typename BaseType::NodesArrayType                                      NodesArrayType;

    typedef typename BaseType::PropertiesType::Pointer                      PropertiesPointerType;
    
    typedef typename std::vector<array_1d<PointType,TDim>>                 ConditionArrayListType;
    
    typedef Line2D2<Point>                                                               LineType;
    
    typedef Triangle3D3<Point>                                                       TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompositionType;

    static constexpr unsigned int NumNodes = (TNumNodesElem == 3 || (TDim == 2 && TNumNodesElem == 4)) ? 2 : TNumNodesElem == 4 ? 3 : 4;

    static constexpr unsigned int MatrixSize = TTensor * (3 * NumNodes);
    
    typedef MortarKinematicVariables<NumNodes>                                   GeneralVariables;
    
    typedef DualLagrangeMultiplierOperators<NumNodes>                                      AeData;
    
    typedef MortarOperator<NumNodes>                                      MortarConditionMatrices;
    
    typedef ExactMortarIntegrationUtility<TDim, NumNodes, false>               IntegrationUtility;
         
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MeshTyingMortarCondition()
        : PairedCondition(),
          mIntegrationOrder(2)
    {}
    
    // Constructor 1
    MeshTyingMortarCondition(
        IndexType NewId, 
        GeometryType::Pointer pGeometry
        ) :PairedCondition(NewId, pGeometry),
           mIntegrationOrder(2)
    {}
    
    // Constructor 2
    MeshTyingMortarCondition(
        IndexType NewId, 
        GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties
        ) :PairedCondition( NewId, pGeometry, pProperties ),
           mIntegrationOrder(2)
    {}
    
    // Constructor 3
    MeshTyingMortarCondition(
        IndexType NewId, 
        GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties, 
        GeometryType::Pointer pMasterGeometry
        )
        :PairedCondition( NewId, pGeometry, pProperties, pMasterGeometry),
         mIntegrationOrder(2)
    {}

    ///Copy constructor
    MeshTyingMortarCondition( MeshTyingMortarCondition const& rOther){}

    /// Destructor.
    ~MeshTyingMortarCondition() override;

    /**
     * Flags related to the element computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );

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
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create( 
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties 
        ) const override;
    
    /**
     * Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeom the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;
        
    /**
     * Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeom the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @param pMasterGeom the paired geometry
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties,
        GeometryType::Pointer pMasterGeom
        ) const override;
        
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector( 
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @param rConditionalDofList
     * @param rCurrentProcessInfo the current process info instance
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

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process information
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;
        
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
     * This data will be used to compute teh derivatives
     */ 
    struct DofData
    {
    public:
        
        // Auxiliar types
        typedef bounded_matrix<double, NumNodes, TTensor>  Type1;
        typedef bounded_matrix<double, NumNodes, NumNodes> Type2;
        
        // The DoF
        Type1 LagrangeMultipliers, u1, u2;
        
        // Ae
        Type2 Ae;
        
        // Default destructor
        ~DofData()= default;
        
        /** 
         * Updating the Slave pair
         * @param GeometryInput The pointer of the current master
         */
        void Initialize(const GeometryType& GeometryInput)
        {            
            // The current Lagrange Multipliers
            u1 = ZeroMatrix(NumNodes, TTensor);
            u2 = ZeroMatrix(NumNodes, TTensor);
            LagrangeMultipliers = ZeroMatrix(NumNodes, TTensor);
        }
        
        // Initialize the Ae components
        void InitializeAeComponents()
        {
            Ae = ZeroMatrix(NumNodes, NumNodes);
        }
        
        /** 
         * Updating the Master pair
         * @param GeometryInput The pointer of the current master
         */
        void UpdateMasterPair(const GeometryType& GeometryInput)
        { 
            /* DoF */
            if (TTensor == 1)
            {
                for (unsigned int i_node = 0; i_node < NumNodes; ++i_node)
                {
                    const double value = GeometryInput[i_node].FastGetSolutionStepValue(TEMPERATURE);
                    u2(i_node, 0) = value;
                }
            }
            else
            {
                for (unsigned int i_node = 0; i_node < NumNodes; ++i_node)
                {
                    const array_1d<double, 3>& value = GeometryInput[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                    for (unsigned int i_dof = 0; i_dof < TTensor; ++i_dof)
                    {
                        u2(i_node, i_dof) = value[i_dof];
                    }
                }
            }
        }

    };
    
    ///@}
    ///@name Protected member Variables
    ///@{

    Flags  mCalculationFlags;                              // Calculation flags
    
    MortarConditionMatrices mrThisMortarConditionMatrices; // The mortar operators
   
    unsigned int mIntegrationOrder;                        // The integration order to consider
    
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
     * @param rLeftHandSideMatrix the condition left hand side matrix
     * @param rRightHandSideVector the condition right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem( 
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
     * This is called during the assembling process in order
     * to calculate the condition right hand side vector only
     * @param rRightHandSideVector the condition right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
     * This is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix the condition left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide( 
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo 
        ) override;
        
    /**
     * Calculates the condition contribution
     */
    void CalculateConditionSystem( 
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& CurrentProcessInfo 
        );
    
    /**
     * Initialize Contact data
     */
    void InitializeDofData(DofData& rDofData);
    
    /**
     * Calculate Ae matrix
     */
    bool CalculateAe( 
        const array_1d<double, 3>& NormalMaster,
        DofData& rDofData,
        GeneralVariables& rVariables,
        ConditionArrayListType& ConditionsPointsSlave,
        IntegrationMethod ThisIntegrationMethod
        );
    
    /**
     * Calculate condition kinematics
     */
    void CalculateKinematics( 
        GeneralVariables& rVariables,
        const DofData& rDofData,
        const array_1d<double, 3>& NormalMaster,
        const PointType& LocalPointDecomp,
        const PointType& LocalPointParent,
        GeometryPointType& GeometryDecomp,
        const bool DualLM = false
        );

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /*
     * Calculates the local contibution of the LHS
     */
    
    bounded_matrix<double, MatrixSize, MatrixSize> CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DofData& rDofData
        );
    
    /*
     * Calculates the local contibution of the LHS
     */
    array_1d<double, MatrixSize> CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DofData& rDofData
        );
    
    /***********************************************************************************/
    /**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
    /***********************************************************************************/

    /*
     * Calculates the values of the shape functions for the master element
     */
    void MasterShapeFunctionValue(
        GeneralVariables& rVariables,
        const array_1d<double, 3>& NormalMaster,
        const PointType& LocalPoint
        );
    
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/
    
    /**
     * It returns theintegration method considered
     */
    
    IntegrationMethod GetIntegrationMethod() override
    {
        // Setting the auxiliar integration points
        switch (mIntegrationOrder) {
        case 1: return GeometryData::GI_GAUSS_1;
        case 2: return GeometryData::GI_GAUSS_2;
        case 3: return GeometryData::GI_GAUSS_3;
        case 4: return GeometryData::GI_GAUSS_4;
        case 5: return GeometryData::GI_GAUSS_5;
        default: return GeometryData::GI_GAUSS_2;
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
