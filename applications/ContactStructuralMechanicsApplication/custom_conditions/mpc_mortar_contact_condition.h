// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_MPC_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_MPC_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "custom_conditions/paired_condition.h"
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"
#include "includes/checks.h"
#include "includes/mortar_classes.h"

/* Utilities */
#include "utilities/exact_mortar_segmentation_utility.h"
#include "custom_utilities/derivatives_utilities.h"

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
 * @class MPCMortarContactCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief MPCMortarContactCondition
 * @details This is a contact condition which employes the mortar method with dual lagrange multiplier in a explicit manner in order to compute the gap/slip and the weights
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes,const SizeType TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) MPCMortarContactCondition
    : public PairedCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPCMortarContactCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPCMortarContactCondition );

    /// Base class definitions
    typedef PairedCondition                                                                      BaseType;

    /// Vector type definition
    typedef typename BaseType::VectorType                                                      VectorType;

    /// Matrix type definition
    typedef typename BaseType::MatrixType                                                      MatrixType;

    /// Index type definition
    typedef typename BaseType::IndexType                                                        IndexType;

    /// Geometry pointer definition
    typedef typename BaseType::GeometryType::Pointer                                  GeometryPointerType;

    /// Nodes array type definition
    typedef typename BaseType::NodesArrayType                                              NodesArrayType;

    /// Properties pointer definition
    typedef typename BaseType::PropertiesType::Pointer                              PropertiesPointerType;

    /// Point definition
    typedef Point                                                                               PointType;

    typedef array_1d<Point,TDim>                                                       ConditionArrayType;

    typedef typename std::vector<ConditionArrayType>                               ConditionArrayListType;

    /// Node type definition
    typedef Node<3>                                                                              NodeType;

    /// Geoemtry type definition
    typedef Geometry<NodeType>                                                               GeometryType;

    // Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType                                IntegrationPointsType;

    typedef Line2D2<PointType>                                                                   LineType;

    typedef Triangle3D3<PointType>                                                           TriangleType;

    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type         DecompositionType;

    typedef MortarKinematicVariables<TNumNodes, TNumNodesMaster>                         GeneralVariables;

    typedef DualLagrangeMultiplierOperators<TNumNodes, TNumNodesMaster>                            AeData;

    typedef MortarOperator<TNumNodes, TNumNodesMaster>                            MortarConditionMatrices;

    typedef ExactMortarIntegrationUtility<TDim, TNumNodes, false, TNumNodesMaster>     IntegrationUtility;

    typedef DerivativesUtilities<TDim, TNumNodes, false, false, TNumNodesMaster> DerivativesUtilitiesType;

    static constexpr IndexType MatrixSize = TDim * (TNumNodes + TNumNodesMaster);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MPCMortarContactCondition()
        : PairedCondition()
    {}

    // Constructor 1
    MPCMortarContactCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        ) :PairedCondition(NewId, pGeometry)
    {}

    // Constructor 2
    MPCMortarContactCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        ) :PairedCondition( NewId, pGeometry, pProperties )
    {}

    // Constructor 3
    MPCMortarContactCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        GeometryType::Pointer pMasterGeometry
        )
        :PairedCondition( NewId, pGeometry, pProperties, pMasterGeometry)
    {}

    ///Copy constructor
    MPCMortarContactCondition( MPCMortarContactCondition const& rOther){}

    /// Destructor.
    ~MPCMortarContactCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element pointer from an arry of nodes
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
     * @brief Creates a new element pointer from an existing geometry
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
     * @brief Creates a new element pointer from an existing geometry
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

   /**
    * @brief Called at the beginning of each solution step
    */
    void Initialize() override;

   /**
    * @brief Called at the beginning of each solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

   /**
    * @brief Called at the beginning of each iteration
    * @param rCurrentProcessInfo the current process info instance
    */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
    * @brief Called at the ending of each solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

   /**
    * @brief Called at the end of each iteration
    * @param rCurrentProcessInfo the current process info instance
    */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @param rConditionalDofList The list of DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rConditionalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;


    /**
    * @brief This method computes the mass matrix
    * @param rMassMatrix The mass matrix to be computed
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Initialize Damping Matrix
     * @param rDampingMatrix The damping matrix to be computed
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the condition contribution in explicit calculation. NodalData is modified Inside the function, so the
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH A CONDITION IS ALLOWED TO WRITE ON ITS NODES. The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo) override;

    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * @brief Get on rVariable a double Value
     * @param rVariable Internal values
     * @param rCurrentProcessInfo The current process information
     * @param rValues The values of interest (doubles)
     */
    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a array_1d Value
     * @param rVariable Internal values
     * @param rCurrentProcessInfo The current process information
     * @param rValues The values of interest (array_1d)
     */
    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a Vector Value
     * @param rVariable Internal values
     * @param rCurrentProcessInfo The current process information
     * @param rValues The values of interest (vector)
     */
    void GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a double Variable
     * @param rVariable Internal values
     * @param rCurrentProcessInfo The current process information
     * @param rOutput The values of interest (doubles)
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a array_1d Variable
     * @param rVariable Internal values
     * @param rCurrentProcessInfo The current process information
     * @param rOutput The values of interest (array_1d)
     */
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector< array_1d<double, 3 > >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Vector Variable
     * @param rVariable Internal values
     * @param rCurrentProcessInfo The current process information
     * @param rOutput The values of interest (vector)
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MPCMortarContactCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MPCMortarContactCondition #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        PrintInfo(rOStream);
        this->GetParentGeometry().PrintData(rOStream);
        this->GetPairedGeometry().PrintData(rOStream);
    }

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

    bool mPreviousMortarOperatorsInitialized = false; /// In order to know iw we need to initialize the previous operators

    MortarConditionMatrices mPreviousMortarOperators; /// These are the mortar operators from the previous converged step, necessary for a consistent definition of the r_gt

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
     * @brief This is called during the assembling process in order
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
     * @brief This is called during the assembling process in order
     * to calculate the condition right hand side vector only
     * @param rRightHandSideVector the condition right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix the condition left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief It returns theintegration method considered
     */
    IntegrationMethod GetIntegrationMethod() override
    {
        // Setting the auxiliar integration points
        const IndexType integration_order = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
        switch (integration_order) {
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

    /**
     * @brief This method updates the database of the dofs of the constraint
     * @param rRelationMatrix The relation matrix between the master/slave DoF
     * @param rConstantVector The vector containing the additional kinematic relationship
     * @param rCurrentProcessInfo The current process instance
     */
    void ConstraintDofDatabaseUpdate(
        Matrix& rRelationMatrix,
        Vector& rConstantVector,
        ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method updates the realtion matrix and the constant vector (frictionless)
     * @param rMortarConditionMatrices Already computed mortar matrices
     * @param rRelationMatrix The relation matrix between the master/slave DoF
     * @param rConstantVector The vector containing the additional kinematic relationship
     * @param rCurrentProcessInfo The current process instance
     * @param DualLM If considering dual LM
     */
    void UpdateConstraintFrictionless(
        MortarConditionMatrices& rMortarConditionMatrices,
        Matrix& rRelationMatrix,
        Vector& rConstantVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool DualLM = true
        );

    /**
     * @brief This method updates the realtion matrix and the constant vector (frictional)
     * @param rMortarConditionMatrices Already computed mortar matrices
     * @param rRelationMatrix The relation matrix between the master/slave DoF
     * @param rConstantVector The vector containing the additional kinematic relationship
     * @param rCurrentProcessInfo The current process instance
     * @param DualLM If considering dual LM
     */
    void UpdateConstraintFrictional(
        MortarConditionMatrices& rMortarConditionMatrices,
        Matrix& rRelationMatrix,
        Vector& rConstantVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool DualLM = true
        );

    /**
     * @brief This method updates the realtion matrix and the constant vector (tying)
     * @param rMortarConditionMatrices Already computed mortar matrices
     * @param rRelationMatrix The relation matrix between the master/slave DoF
     * @param rConstantVector The vector containing the additional kinematic relationship
     * @param rCurrentProcessInfo The current process instance
     * @param DualLM If considering dual LM
     */
    void UpdateConstraintTying(
        MortarConditionMatrices& rMortarConditionMatrices,
        Matrix& rRelationMatrix,
        Vector& rConstantVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool DualLM = true
        );

    /**
     * @brief It computes the previous mortar operators
     * @param rCurrentProcessInfo The current process instance
     */
    void ComputePreviousMortarOperators(ProcessInfo& rCurrentProcessInfo);

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
        rSerializer.save("PreviousMortarOperators", mPreviousMortarOperators);
        rSerializer.save("PreviousMortarOperatorsInitialized", mPreviousMortarOperatorsInitialized);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PairedCondition );
        rSerializer.load("PreviousMortarOperators", mPreviousMortarOperators);
        rSerializer.load("PreviousMortarOperatorsInitialized", mPreviousMortarOperatorsInitialized);
    }

    ///@}

}; // Class MPCMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MPC_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined
