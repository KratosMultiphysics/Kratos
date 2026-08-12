// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ /
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
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
// #include "custom_utilities/logging_settings.hpp"

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
 * @class ALM3dMortarFrictionlessCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief ALM3dMortarFrictionlessCondition
 * @details Augmented Lagrangian Mortar Frictionless Condition in 3D.
 * This condition has the parent (slave) and the paired (master) geometries, it then
 * performs the geometrical segmentation and the mortar integration. It is able to compute the local system, the right hand side and the left hand side.
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ALM3dMortarFrictionlessCondition
    : public PairedCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ALM3dMortarFrictionlessCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ALM3dMortarFrictionlessCondition);

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

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    // Type definition for integration methods
    using IntegrationPointsType = typename GeometryType::IntegrationPointsArrayType;

    /// The type of points belongs to be considered
    using BelongType = typename PointBelongsTriangle3D3N;

    /// The definition of the point with belonging
    using PointBelongType = PointBelong<3, 3>;

    /// Type definition for the geometry with point belonging
    // using GeometryPointBelongType = Geometry<PointBelongType>;

    /// Type definition for an array of points with belonging
    using ConditionArrayType = array_1d<PointBelongType, 3>;

    /// Type definition for a list of arrays of points with belonging
    using ConditionArrayListType = typename std::vector<ConditionArrayType>;

    /// Type definition for a line in 2D
    // using LineType = Line2D2<PointType>;

    /// Type definition for a triangle in 3D
    using TriangleType = Triangle3D3<PointType>;

    /// The decomposition type
    // using DecompositionType = typename std::conditional<TDim == 2, LineType, TriangleType >::type;
    using DecompositionType = typename TriangleType;

    /// The derivative data type
    // using DerivativeDataType = typename std::conditional<TFrictional == FrictionalCase::FRICTIONAL || TFrictional == FrictionalCase::FRICTIONAL_PENALTY, DerivativeDataFrictional<TDim, TNumNodes, TNumNodesMaster>, DerivativeData<TDim, TNumNodes, TNumNodesMaster> >::type;

    /// The matrix size definition
    // static constexpr IndexType MatrixSize = (TFrictional == FrictionalCase::FRICTIONLESS) ? TDim * (TNumNodesMaster + TNumNodes) + TNumNodes : (TFrictional == FrictionalCase::FRICTIONLESS_COMPONENTS || TFrictional == FrictionalCase::FRICTIONAL) ? TDim * (TNumNodesMaster + TNumNodes + TNumNodes) :  TDim * (TNumNodesMaster + TNumNodes);

    /// The definition of the frictional flag
    // static constexpr bool IsFrictional  = (TFrictional == FrictionalCase::FRICTIONAL || TFrictional == FrictionalCase::FRICTIONAL_PENALTY) ? true: false;

    /// Type definition for general variables with derivatives
    // using GeneralVariables = MortarKinematicVariablesWithDerivatives<TDim, TNumNodes, TNumNodesMaster>;

    // /// Type definition for AE data with derivatives
    // using AeData = DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, IsFrictional, TNumNodesMaster>;

    // /// Type definition for mortar condition matrices with derivatives
    // using MortarConditionMatrices = MortarOperatorWithDerivatives<TDim, TNumNodes, IsFrictional, TNumNodesMaster>;

    /// Type definition for integration utility with derivatives
    using IntegrationUtility = ExactMortarIntegrationUtility<3, 3, true, 3>;

    // /// Type definition for derivatives utilities with derivatives
    using DerivativesUtilitiesType = DerivativesUtilities<3, 3, false, false, 3>;

    // The threshold coefficient considered for checking
    static constexpr double CheckThresholdCoefficient = 1.0e-12;
    static constexpr double MinIntegrationAreaRatioTolerance = 1.0e-8;
    static constexpr double AugmentedLMtolerance = 1.0e-12;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ALM3dMortarFrictionlessCondition()
        : PairedCondition()
    {}

    // Constructor 1
    ALM3dMortarFrictionlessCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        ) :PairedCondition(NewId, pGeometry)
    {}

    // Constructor 2
    ALM3dMortarFrictionlessCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        ) :PairedCondition( NewId, pGeometry, pProperties )
    {}

    // Constructor 3
    ALM3dMortarFrictionlessCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        GeometryType::Pointer pMasterGeometry
        )
        :PairedCondition( NewId, pGeometry, pProperties, pMasterGeometry)
    {}

    ///Copy constructor
    ALM3dMortarFrictionlessCondition( ALM3dMortarFrictionlessCondition const& rOther){}

    /// Destructor.
    ~ALM3dMortarFrictionlessCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

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

    /******************************************************************/
    /********** AUXILIARY METHODS FOR GENERAL CALCULATIONS ************/
    /******************************************************************/

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
     * @param rConditionalDofList The list of DoF
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& rConditionalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

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

    MatrixType CalculateSlaveNormalDerivatives(const MatrixType& J, const MatrixType& DN_De)
    {
        MatrixType dN_dAs(3, 3 * 3);
        array_1d<array_1d<double, 3>, 3*3> dN_dAs_array;
        dN_dAs_array =  DerivativesUtilitiesType::GPDeltaNormalSlave(J, DN_De);
        for (IndexType j = 0; j < 9; ++j) { // cols
            for (IndexType i = 0; i < 3; ++i) { // rows
                dN_dAs(i, j) = dN_dAs_array[j][i];
            }
        }
        return dN_dAs;
    }

    /**
     * @brief Project a point over a line/plane following an arbitrary direction
     * @tparam TGeometryType The type of the geometry
     * @param rGeom The geometry where to be projected
     * @param rPointToProject The point to be projected
     * @param rPointProjected The point pojected over the plane
     * @param rNormal The normal of the geometry
     * @param rVector The direction to project
     * @param EchoLevel If we want debugging info we should consider greater than 0
     * @return if the projection falls inside the destination geometry or not.
     */
    bool FastProjectDirection(
        const GeometryType &rGeometryToProject,
        const PointType &rPointToProject,
        PointType &rPointProjected,
        const array_1d<double, 3> &rNormal);

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
        buffer << "ALM3dMortarFrictionlessCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ALM3dMortarFrictionlessCondition #" << this->Id();
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
     * @brief Calculates the condition contribution
     */
    void CalculateConditionSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& CurrentProcessInfo,
        const bool ComputeLHS = true,
        const bool ComputeRHS = true
        );

    /**
     * @brief It checks if the element is isolated or not
     * @param DeltaTime The increment of time in each time step
     * @param HalfJump If the increment of time considered is just half or the whole time step
     */
    bool CheckIsolatedElement(
        const double DeltaTime,
        const bool HalfJump = true
        );

    IndexType GetSystemSize() const
    {
        IndexType slave_nodes = this->GetParentGeometry().PointsNumber();
        IndexType master_nodes = this->GetPairedGeometry().PointsNumber();
        return 3 * (master_nodes + slave_nodes) + slave_nodes;
    }

    void CalculateInterpolationMatrices(
        const Vector &rN_slave,
        const Vector &rN_master,
        const Vector &rN_LMgeom, // coming from the geometry
        Matrix &rNs,
        Matrix &rNm,
        Vector &rN_LM);

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need the values of the time derivatives of any of the dof
     * set by the element. If the derivatives do not exist can set to zero
     * these methods are: MANDATORY ( when compatibility with dynamics is required )
     */

    /**
     * Getting method to obtain the variable which defines the degrees of freedom
     */
    void GetValuesVector(Vector &rValues, int Step = 0) const;

    void AddRightHandSideContribution(
        VectorType &rRightHandSideVector,
        const double k,
        const double penalty,
        const double gap,
        const double LM,
        const double integration_weight,
        const Matrix &rNs,
        const Matrix &rNm,
        const Vector &rN_LM,  // Phi
        const Vector &rNormal, // slave normal vector
        const bool active_contact
    );

    void AddLeftHandSideContribution(
        MatrixType &rLeftHandSideMatrix,
        const double k,
        const double penalty,
        const double integration_weight,
        const double AugmentedLM,
        const Matrix &rNs,
        const Matrix &rNm,
        const Vector &rN_LM,  // Phi
        const Vector &rNormal, // slave normal vector
        const Vector& rDeltaX, // Xs - Xm
        const bool active_contact,
        const Matrix& rDnda_slave
    );

    /**
     * @brief It returns theintegration method considered
     */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return GeometryData::IntegrationMethod::GI_GAUSS_5;
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

    Vector mSlaveNormal; /// The normal of the slave surface

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PairedCondition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PairedCondition);
    }

    ///@}

}; // Class ALM3dMortarFrictionlessCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.