// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ /
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

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
    using SizeType = std::size_t;

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
    using BaseType = PairedCondition;

    /// Vector type definition
    using VectorType = typename BaseType::VectorType;

    /// Matrix type definition
    using MatrixType = typename BaseType::MatrixType;

    /// Index type definition
    using IndexType = typename BaseType::IndexType;

    /// Geometry pointer definition
    using GeometryPointerType = typename BaseType::GeometryType::Pointer;

    /// Nodes array type definition
    using NodesArrayType = typename BaseType::NodesArrayType;

    /// Properties pointer definition
    using PropertiesPointerType = typename BaseType::PropertiesType::Pointer;

    /// Point definition
    using PointType = Point;

    /// Array type for condition with points
    using ConditionArrayType = array_1d<Point, TDim>;

    /// Type definition for a list of condition arrays
    using ConditionArrayListType = typename std::vector<ConditionArrayType>;

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    // Type definition for integration methods
    using IntegrationPointsType = typename GeometryType::IntegrationPointsArrayType;

    /// Line type definition
    using LineType = Line2D2<PointType>;

    /// Triangle type definition
    using TriangleType = Triangle3D3<PointType>;

    /// Type definition for decomposition based on dimension
    using DecompositionType = typename std::conditional<TDim == 2, LineType, TriangleType>::type;

    /// Type definition for general variables
    using GeneralVariables = MortarKinematicVariables<TNumNodes, TNumNodesMaster>;

    /// Type definition for AE data
    using AeData = DualLagrangeMultiplierOperators<TNumNodes, TNumNodesMaster>;

    /// Type definition for mortar condition matrices
    using MortarConditionMatrices = MortarOperator<TNumNodes, TNumNodesMaster>;

    /// Type definition for integration utility
    using IntegrationUtility = ExactMortarIntegrationUtility<TDim, TNumNodes, false, TNumNodesMaster>;

    /// Type definition for derivatives utilities
    using DerivativesUtilitiesType = DerivativesUtilities<TDim, TNumNodes, false, false, TNumNodesMaster>;

    /// Constant expression for matrix size
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
     * @brief Creates a new element pointer from an array of nodes
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
    * @param rCurrentProcessInfo the current process info instance
    */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

   /**
    * @brief Called at the beginning of each solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

   /**
    * @brief Called at the beginning of each iteration
    * @param rCurrentProcessInfo the current process info instance
    */
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * @brief Called at the ending of each solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

   /**
    * @brief Called at the end of each iteration
    * @param rCurrentProcessInfo the current process info instance
    */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @param rConditionalDofList The list of DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rConditionalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;


    /**
    * @brief This method computes the mass matrix
    * @param rMassMatrix The mass matrix to be computed
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Initialize Damping Matrix
     * @param rDampingMatrix The damping matrix to be computed
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the condition contribution in explicit calculation. NodalData is modified Inside the function, so the
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH A CONDITION IS ALLOWED TO WRITE ON ITS NODES. The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo) override;

    /******************************************************************/
    /********** AUXILIARY METHODS FOR GENERAL CALCULATIONS ************/
    /******************************************************************/

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
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order
     * to calculate the condition right hand side vector only
     * @param rRightHandSideVector the condition right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix the condition left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief It returns theintegration method considered
     */
    IntegrationMethod GetIntegrationMethod() const override
    {
        // Setting the auxiliary integration points
        const IndexType integration_order = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
        switch (integration_order) {
        case 1: return GeometryData::IntegrationMethod::GI_GAUSS_1;
        case 2: return GeometryData::IntegrationMethod::GI_GAUSS_2;
        case 3: return GeometryData::IntegrationMethod::GI_GAUSS_3;
        case 4: return GeometryData::IntegrationMethod::GI_GAUSS_4;
        case 5: return GeometryData::IntegrationMethod::GI_GAUSS_5;
        default: return GeometryData::IntegrationMethod::GI_GAUSS_2;
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
        const ProcessInfo& rCurrentProcessInfo
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
        const ProcessInfo& rCurrentProcessInfo,
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
        const ProcessInfo& rCurrentProcessInfo,
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
        const ProcessInfo& rCurrentProcessInfo,
        const bool DualLM = true
        );

    /**
     * @brief It computes the previous mortar operators
     * @param rCurrentProcessInfo The current process instance
     */
    void ComputePreviousMortarOperators(const ProcessInfo& rCurrentProcessInfo);

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
