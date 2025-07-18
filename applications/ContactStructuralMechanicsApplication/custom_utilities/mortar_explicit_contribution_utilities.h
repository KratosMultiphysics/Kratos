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
#include "utilities/math_utils.h"
#include "custom_conditions/paired_condition.h"
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
 * @class MortarExplicitContributionUtilities
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This namespace includes several utilities necessaries for the computation of the explicit contribution of the mortar conditions
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If we are solving a frictional or frictionless problem
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template<std::size_t TDim, std::size_t TNumNodes, const FrictionalCase TFrictional, const bool TNormalVariation, std::size_t TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) MortarExplicitContributionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// The size type definition
    using SizeType = std::size_t;

    /// The index type definition
    using IndexType = std::size_t;

    /// Point definition
    using PointType = Point;

    /// Geometry point type definition
    using GeometryPointType = Geometry<PointType>;

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    // Type definition for integration methods
    using IntegrationPointsType = GeometryType::IntegrationPointsArrayType;

    /// The type of points belongs to be considered
    using BelongType = typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, typename std::conditional<TNumNodesMaster == 3, PointBelongsTriangle3D3N, PointBelongsTriangle3D3NQuadrilateral3D4N>::type, typename std::conditional<TNumNodesMaster == 3, PointBelongsQuadrilateral3D4NTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type>::type;

    /// The definition of the point with belonging
    using PointBelongType = PointBelong<TNumNodes, TNumNodesMaster>;

    /// Type definition for the geometry of the point belong
    using GeometryPointBelongType = Geometry<PointBelongType>;

    /// Array type for the conditions
    using ConditionArrayType = array_1d<PointBelongType, TDim>;

    /// Type definition for the array list of conditions
    using ConditionArrayListType = typename std::vector<ConditionArrayType>;

    /// Type definition for the line
    using LineType = Line2D2<PointType>;

    /// Type definition for the triangle
    using TriangleType = Triangle3D3<PointType>;

    /// Type definition for the decomposition based on dimension
    using DecompositionType = typename std::conditional<TDim == 2, LineType, TriangleType>::type;

    /// Type definition for the derivative data based on frictional case
    using DerivativeDataType = typename std::conditional<TFrictional == FrictionalCase::FRICTIONAL || TFrictional == FrictionalCase::FRICTIONAL_PENALTY, DerivativeDataFrictional<TDim, TNumNodes, TNumNodesMaster>, DerivativeData<TDim, TNumNodes, TNumNodesMaster>>::type;

    /// The size of the matrix
    static constexpr IndexType MatrixSize = (TFrictional == FrictionalCase::FRICTIONLESS) ? TDim * (TNumNodesMaster + TNumNodes) + TNumNodes : (TFrictional == FrictionalCase::FRICTIONLESS_COMPONENTS || TFrictional == FrictionalCase::FRICTIONAL) ? TDim * (TNumNodesMaster + TNumNodes + TNumNodes) : TDim * (TNumNodesMaster + TNumNodes);

    /// Flag indicating if the case is frictional
    static constexpr bool IsFrictional = (TFrictional == FrictionalCase::FRICTIONAL || TFrictional == FrictionalCase::FRICTIONAL_PENALTY) ? true : false;

    /// Type definition for the general variables
    using GeneralVariables = MortarKinematicVariablesWithDerivatives<TDim, TNumNodes, TNumNodesMaster>;

    /// Type definition for the Ae data
    using AeData = DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, IsFrictional, TNumNodesMaster>;

    /// Type definition for the mortar condition matrices
    using MortarConditionMatrices = MortarOperatorWithDerivatives<TDim, TNumNodes, IsFrictional, TNumNodesMaster>;

    /// Type definition for the integration utility
    using IntegrationUtility = ExactMortarIntegrationUtility<TDim, TNumNodes, true, TNumNodesMaster>;

    /// Type definition for the derivatives utilities
    using DerivativesUtilitiesType = DerivativesUtilities<TDim, TNumNodes, IsFrictional, TNormalVariation, TNumNodesMaster>;

    /// Type definition for the integration method
    using IntegrationMethod = GeometryData::IntegrationMethod;

    // The threshold coefficient considered for checking
    static constexpr double CheckThresholdCoefficient = 1.0e-12;



    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the explicit contributions of the mortar contact conditions
     * @details This method is created in order to avoid duplicated code
     * @param pCondition The condition pointer to compute the explicit contribution
     * @param rCurrentProcessInfo The current instance process info
     * @param IntegrationOrder The integration order of the utility
     * @param AxisymmetricCase If consider the axisymmetric coefficient
     * @param ComputeNodalArea If the contribution of the nodal are must be computed
     * @param ComputeDualLM If condider dual LM to begin with
     * @param rAreaVariable The nodal area variable
     * @return The mortar operators
     */
    static MortarConditionMatrices AddExplicitContributionOfMortarCondition(
        PairedCondition* pCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const IndexType IntegrationOrder = 2,
        const bool AxisymmetricCase = false,
        const bool ComputeNodalArea = false,
        const bool ComputeDualLM = true,
        const Variable<double>& rAreaVariable = NODAL_AREA
        );
    /**
     * @brief This method computes the explicit contributions of the mortar contact conditions
     * @details This method is created in order to avoid duplicated code
     * @param pCondition The condition pointer to compute the explicit contribution
     * @param rCurrentProcessInfo The current instance process info
     * @param rPreviousMortarOperators The previous mortar operators
     * @param IntegrationOrder The integration order of the utility
     * @param AxisymmetricCase If consider the axisymmetric coefficient
     * @param ComputeNodalArea If the contribution of the nodal are must be computed
     * @param ComputeDualLM If condider dual LM to begin with
     * @param rAreaVariable The nodal area variable
     * @param ConsiderObjetiveFormulation If the objetive formulation is considered always
     * @return The mortar operators
     */
    static MortarConditionMatrices AddExplicitContributionOfMortarFrictionalCondition(
        PairedCondition* pCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const MortarOperator<TNumNodes, TNumNodesMaster>& rPreviousMortarOperators,
        const IndexType IntegrationOrder = 2,
        const bool AxisymmetricCase = false,
        const bool ComputeNodalArea = false,
        const bool ComputeDualLM = true,
        const Variable<double>& rAreaVariable = NODAL_AREA,
        const bool ConsiderObjetiveFormulation = false
        );

    /**
     * @brief Calculate the operator Ae for the dual LM, without taking into account derivatives
     * @details This can be used in the mortar conditions
     * @param rSlaveGeometry The geometry of the slave side
     * @param rVariables Container of the jacobians, shape functions, etc...
     * @param rConditionsPointsSlave Container of the jacobians, shape functions, etc...
     * @param rAe The dual LM operator
     * @param rIntegrationMethod The integration method considered
     * @param AxiSymCoeff The coefficient of axisymmetry
     */
    static bool ExplicitCalculateAe(
        const GeometryType& rSlaveGeometry,
        GeneralVariables& rVariables,
        const ConditionArrayListType& rConditionsPointsSlave,
        BoundedMatrix<double, TNumNodes, TNumNodes>& rAe,
        const IntegrationMethod& rIntegrationMethod,
        const double AxiSymCoeff = 1.0
        );

    /**
     * @brief Calculate condition kinematics (shape functions, jacobians, ...), without taking into account derivatives
     * @details This can be used in the mortar conditions
     * @param pCondition The pointer to the condition
     * @param rVariables Container of the jacobians, shape functions, etc...
     * @param rAe The dual LM operator
     * @param rNormalMaster The normal of the master side
     * @param rLocalPointDecomp The local points of the decomposed geometries using the mortar segmentation
     * @param rLocalPointParent The local points of the parent geometry
     * @param rGeometryDecomp The geometry decomposed
     * @param DualLM If the dual LM is considered or not
     */
    static void ExplicitCalculateKinematics(
        const PairedCondition* pCondition,
        GeneralVariables& rVariables,
        const BoundedMatrix<double, TNumNodes, TNumNodes>& rAe,
        const array_1d<double, 3>& rNormalMaster,
        const PointType& rLocalPointDecomp,
        const PointType& rLocalPointParent,
        const GeometryPointType& rGeometryDecomp,
        const bool DualLM = true
        );

    /**
     * @brief This method computes the nodal area
     * @details This method is created in order to avoid duplicated code
     * @param pCondition The condition pointer to compute the explicit contribution
     * @param rCurrentProcessInfo The current instance process info
     * @param rAreaVariable The nodal area variable
     * @param IntegrationOrder The integration order of the utility
     * @param AxisymmetricCase If consider the axisymmetric coefficient
     * @return True is dual LM, false otherwise
     */
    static void ComputeNodalArea(
        PairedCondition* pCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const Variable<double>& rAreaVariable = NODAL_AREA,
        const IndexType IntegrationOrder = 2,
        const bool AxisymmetricCase = false
        );

    /**
     * @brief This method computes the previous mortar operators
     * @details This method is created in order to avoid duplicated code
     * @param pCondition The condition pointer to compute the explicit contribution
     * @param rCurrentProcessInfo The current instance process info
     * @param rPreviousMortarOperators The previous mortar operators
     * @param IntegrationOrder The integration order of the utility
     * @param AxisymmetricCase If consider the axisymmetric coefficient
     * @param ComputeNodalArea If the contribution of the nodal are must be computed
     * @param ComputeDualLM If condider dual LM to begin with
     * @param rAreaVariable The nodal area variable
     * @return True is dual LM, false otherwise
     */
    static bool ComputePreviousMortarOperators(
        PairedCondition* pCondition,
        const ProcessInfo& rCurrentProcessInfo,
        MortarOperator<TNumNodes, TNumNodesMaster>& rPreviousMortarOperators,
        const IndexType IntegrationOrder = 2,
        const bool AxisymmetricCase = false,
        const bool ComputeNodalArea = false,
        const bool ComputeDualLM = true,
        const Variable<double>& rAreaVariable = NODAL_AREA
        );

    /**
     * @brief Calculate condition kinematics (shape functions, jacobians, ...)
     * @details This can be used in the mortar conditions
     * @param pCondition The pointer to the condition
     * @param rVariables Container of the jacobians, shape functions, etc...
     * @param rDerivativeData The container of the directional derivatives of the mortar parameters
     * @param rNormalMaster The normal of the master side
     * @param rLocalPointDecomp The local points of the decomposed geometries using the mortar segmentation
     * @param rLocalPointParent The local points of the parent geometry
     * @param rGeometryDecomp The geometry decomposed
     * @param DualLM If the dual LM is considered or not
     */
    static void CalculateKinematics(
        const PairedCondition* pCondition,
        GeneralVariables& rVariables,
        const DerivativeDataType& rDerivativeData,
        const array_1d<double, 3>& rNormalMaster,
        const PointType& rLocalPointDecomp,
        const PointType& rLocalPointParent,
        const GeometryPointType& rGeometryDecomp,
        const bool DualLM = true
        );

    /**
     * @brief Calculates the values of the shape functions for the master element
     * @param pCondition The pointer to the condition
     * @param rVariables Container of the jacobians, shape functions, etc...
     * @param rNormalMaster The normal of the master side
     * @param rLocalPoint The current local point
     */
    static void MasterShapeFunctionValue(
        const PairedCondition* pCondition,
        GeneralVariables& rVariables,
        const array_1d<double, 3>& rNormalMaster,
        const PointType& rLocalPoint
        );

}; // class MortarExplicitContributionUtilities

namespace AuxiliaryOperationsUtilities
{
    /**
     * @brief This functions computes the integration weight to consider
     * @param pCondition The condition pointer to compute the explicit contribution
     * @param rNSlave The shape functions of the slave side
     */
    double KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) GetAxisymmetricCoefficient(
        const PairedCondition* pCondition,
        const Vector& rNSlave
        );

    /**
     * @brief Calculates the radius of axisymmetry
     * @param pCondition The condition pointer to compute the explicit contribution
     * @param rNSlave The shape functions of the slave side
     * @return Radius: The radius of axisymmetry
     */
    double KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) CalculateRadius(
        const PairedCondition* pCondition,
        const Vector& rNSlave
        );
}

}  // namespace Kratos
