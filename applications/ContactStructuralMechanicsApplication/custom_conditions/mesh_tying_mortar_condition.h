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

    /// The definition of the size type
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class MeshTyingMortarCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief MeshTyingMortarCondition
 * @details This is a mesh tying condition which employes the mortar method with dual lagrange multiplier
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodesElem The number of nodes of the slave
 * @tparam TNumNodesElemMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodesElem, const SizeType TNumNodesElemMaster = TNumNodesElem>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) MeshTyingMortarCondition
    : public PairedCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MeshTyingMortarCondition
    KRATOS_CLASS_POINTER_DEFINITION( MeshTyingMortarCondition );

    /// Base class definitions
    typedef PairedCondition                                                               BaseType;

    /// Vector type definition
    typedef typename BaseType::VectorType                                               VectorType;

    /// Matrix type definition
    typedef typename BaseType::MatrixType                                               MatrixType;

    /// Index type definition
    typedef typename BaseType::IndexType                                                 IndexType;

    /// Geometry pointer definition
    typedef typename BaseType::GeometryType::Pointer                           GeometryPointerType;

    /// Nodes array type definition
    typedef typename BaseType::NodesArrayType                                       NodesArrayType;

    /// Properties pointer definition
    typedef typename BaseType::PropertiesType::Pointer                       PropertiesPointerType;

    /// Point definition
    typedef Point                                                                        PointType;

    /// Node type definition
    typedef Node<3>                                                                       NodeType;

    /// Geoemtry type definition
    typedef Geometry<NodeType>                                                        GeometryType;

    // Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType                         IntegrationPointsType;

    // Type definition of the components of an array_1d
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > Array1DComponentsType;

    typedef typename std::vector<array_1d<PointType,TDim>>                  ConditionArrayListType;

    typedef Line2D2<Point>                                                                LineType;

    typedef Triangle3D3<Point>                                                        TriangleType;

    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type  DecompositionType;

    static constexpr SizeType NumNodes = (TNumNodesElem == 3 || (TDim == 2 && TNumNodesElem == 4)) ? 2 : TNumNodesElem == 4 ? 3 : 4;

    static constexpr SizeType NumNodesMaster = (TNumNodesElemMaster == 3 || (TDim == 2 && TNumNodesElemMaster == 4)) ? 2 : TNumNodesElemMaster == 4 ? 3 : 4;

    typedef BoundedMatrix<double, NumNodes, NumNodes>                                  MatrixDualLM;

    typedef MortarKinematicVariables<NumNodes, NumNodesMaster>                     GeneralVariables;

    typedef DualLagrangeMultiplierOperators<NumNodes, NumNodesMaster>                        AeData;

    typedef MortarOperator<NumNodes, NumNodesMaster>                        MortarConditionMatrices;

    typedef ExactMortarIntegrationUtility<TDim, NumNodes, false, NumNodesMaster> IntegrationUtility;

    ///@}
    ///@name  Enum's
    ///@{

    enum TensorValue {ScalarValue = 1, Vector2DValue = 2, Vector3DValue = 3};

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
    * @param rCurrentProcessInfo The current process info instance
    */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

   /**
    * Called at the beginning of each iteration
    * @param rCurrentProcessInfo The current process info instance
    */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
    * Called at the ending of each solution step
    * @param rCurrentProcessInfo The current process info instance
    */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

   /**
    * Called at the end of each iteration
    * @param rCurrentProcessInfo The current process info instance
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
     * @param rConditionalDofList The list of DOFs
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
    template< const TensorValue TTensor >
    struct DofData
    {
    public:

        // Auxiliar types
        typedef BoundedMatrix<double, NumNodes, TTensor>  MatrixUnknownSlave;
        typedef BoundedMatrix<double, NumNodesMaster, TTensor>  MatrixUnknownMaster;

        // The DoF
        MatrixUnknownSlave LagrangeMultipliers, u1;
        MatrixUnknownMaster u2;

        // Default destructor
        ~DofData()= default;

        /**
         * Updating the Slave pair
         * @param rGeometryInput The pointer of the current master
         */
        void Initialize(const GeometryType& rGeometryInput)
        {
            // The current Lagrange Multipliers
            u1 = ZeroMatrix(NumNodes, TTensor);
            u2 = ZeroMatrix(NumNodesMaster, TTensor);
            LagrangeMultipliers = ZeroMatrix(NumNodes, TTensor);
        }

        /**
         * @brief Updating the Master pair
         * @param rGeometryInput The pointer of the current master
         * @param rDoubleVariables The list of double variables
         * @param rArray1DVariables The list of components array1d
         */
        void UpdateMasterPair(
            const GeometryType& rGeometryInput,
            std::vector<Variable<double>>& rDoubleVariables,
            std::vector<Variable<array_1d<double, 3>>>& rArray1DVariables
            )
        {
            /* DoF */
            if (TTensor == 1) {
                for (IndexType i_node = 0; i_node < NumNodesMaster; ++i_node) {
                    u2(i_node, 0) = rGeometryInput[i_node].FastGetSolutionStepValue(rDoubleVariables[0]);
                }
            } else {
                for (IndexType i_node = 0; i_node < NumNodesMaster; ++i_node) {
                    const array_1d<double, 3>& value = rGeometryInput[i_node].FastGetSolutionStepValue(rArray1DVariables[0]);
                    for (IndexType i_dof = 0; i_dof < TTensor; ++i_dof) {
                        u2(i_node, i_dof) = value[i_dof];
                    }
                }
            }
        }

    };

    ///@}
    ///@name Protected member Variables
    ///@{

    Flags  mCalculationFlags;                                      /// Calculation flags

    MortarConditionMatrices mrThisMortarConditionMatrices;         /// The mortar operators

    IndexType mIntegrationOrder;                                   /// The integration order to consider

    std::vector<Variable<double>> mDoubleVariables;                /// The list of double variables

    std::vector<Variable<array_1d<double, 3>>> mArray1DVariables;  /// The list of components array1d

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
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * Initialize Contact data
     */
    template< const TensorValue TTensor >
    void InitializeDofData(DofData<TTensor>& rDofData)
    {
        // Slave element info
        rDofData.Initialize(GetGeometry());

        if (TTensor == ScalarValue) {
            for (IndexType i_node = 0; i_node < NumNodes; i_node++) {
                const double value = GetGeometry()[i_node].FastGetSolutionStepValue(mDoubleVariables[0]);
                const double lm = GetGeometry()[i_node].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
                rDofData.u1(i_node, 0) = value;
                rDofData.LagrangeMultipliers(i_node, 0) = lm;
            }
        } else {
            for (IndexType i_node = 0; i_node < NumNodes; i_node++) {
                const array_1d<double, 3>& value = GetGeometry()[i_node].FastGetSolutionStepValue(mArray1DVariables[0]);
                const array_1d<double, 3>& lm = GetGeometry()[i_node].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                for (IndexType i_dof = 0; i_dof < TDim; i_dof++) {
                    rDofData.u1(i_node, i_dof) = value[i_dof];
                    rDofData.LagrangeMultipliers(i_node, i_dof) = lm[i_dof];
                }
            }
        }
    }

    /**
     * Calculate Ae matrix
     */
    bool CalculateAe(
        const array_1d<double, 3>& rNormalMaster,
        MatrixDualLM& rAe,
        GeneralVariables& rVariables,
        ConditionArrayListType& rConditionsPointsSlave,
        IntegrationMethod ThisIntegrationMethod
        );

    /**
     * Calculate condition kinematics
     */
    void CalculateKinematics(
        GeneralVariables& rVariables,
        const MatrixDualLM& rAe,
        const array_1d<double, 3>& rNormalMaster,
        const PointType& rLocalPointDecomp,
        const PointType& rLocalPointParent,
        GeometryPointType& rGeometryDecomp,
        const bool DualLM = false
        );

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /**
     * @brief Calculates the local contibution of the LHS
     * @param rLocalLHS The local LHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDofData The class containing all the information needed in order to compute the jacobian
     */
    template< const TensorValue TTensor >
    void CalculateLocalLHS(
        Matrix& rLocalLHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DofData<TTensor>& rDofData
        );

    /**
     * @brief Calculates the local contibution of the LHS
     * @param rLocalRHS The local RHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDofData The class containing all the information needed in order to compute the jacobian
     */
    template< const TensorValue TTensor >
    void CalculateLocalRHS(
        Vector& rLocalRHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DofData<TTensor>& rDofData
        );

    /***********************************************************************************/
    /**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
    /***********************************************************************************/

    /*
     * Calculates the values of the shape functions for the master element
     */
    void MasterShapeFunctionValue(
        GeneralVariables& rVariables,
        const array_1d<double, 3>& rNormalMaster,
        const PointType& rLocalPoint
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PairedCondition );
        rSerializer.save("CalculationFlags", mCalculationFlags);
        rSerializer.save("IntegrationOrder", mIntegrationOrder);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PairedCondition );
        rSerializer.load("CalculationFlags", mCalculationFlags);
        rSerializer.load("IntegrationOrder", mIntegrationOrder);
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
