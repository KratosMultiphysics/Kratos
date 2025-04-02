//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/plane_3d.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "modified_shape_functions/modified_shape_functions.h"

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

namespace ShiftedBoundaryUtilityInternals {

    //TODO
    template <std::size_t TDim>
    double CalculatePointDistance(
        const Geometry<Node>& rObjectGeometry,
        const Point& rPoint);

    template <std::size_t TDim>
    Plane3D CreateIntersectionPlane(const std::vector<array_1d<double,3>>& rIntPtsVector);

    void VoigtTransformForProduct(
        const array_1d<double,3>& rVector,
        BoundedMatrix<double,2,3>& rVoigtMatrix);

    void VoigtTransformForProduct(
        const array_1d<double,3>& rVector,
        BoundedMatrix<double,3,6>& rVoigtMatrix);

    template <std::size_t TDim>
    void CalculateStrainMatrix(
        const Matrix& rDN_DX,
        const std::size_t& NumNodes,
        Matrix& rB);

}  // namespace ShiftedBoundaryUtilityInternals

///@}
///@name Kratos Classes
///@{

/**
 * @brief Utilities for the SBM-WTE extension operator calculation
 * This class provides the utilities for the calculation of the extension operator
 * in the Shifted Boundary Method Without Taylor Expansions allowing for a discontinuous interface (thin-walled structure),
 * which is defined by points on the surface of the interface (e.g. IGA integration points).
 */
class KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility
{
public:

    ///@name Type Definitions
    ///@{

    enum class ExtensionOperator
    {
        MLS,
        RBF
    };

    using PointDistanceFunctionType = std::function<double(const Geometry<Node>&, const Point&)>;
    using IntersectionPlaneConstructorType = std::function<Plane3D(const std::vector<array_1d<double,3>>&)>;

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using ElementType = Element;

    using GeometryType = ModelPart::GeometryType;

    using ShapeFunctionsGradientsType = GeometryType::ShapeFunctionsGradientsType;

    using MeshlessShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

    using ElementSizeFunctionType = std::function<double(const GeometryType&)>;

    using SkinPointsDataVectorType = DenseVector<std::tuple< array_1d<double,3>, array_1d<double,3>, std::size_t >>; // vector of position, area normal and ID of skin points

    using SkinPointsToElementsMapType = std::unordered_map<ElementType::Pointer, SkinPointsDataVectorType, SharedPointerHasher<ElementType::Pointer>, SharedPointerComparator<ElementType::Pointer>>;

    using SidesVectorToElementsMapType = std::unordered_map<ElementType::Pointer, Vector, SharedPointerHasher<ElementType::Pointer>, SharedPointerComparator<ElementType::Pointer>>;

    using AverageSkinToElementsMapType = std::unordered_map<ElementType::Pointer, std::pair<array_1d<double,3>, array_1d<double,3>>, SharedPointerHasher<ElementType::Pointer>, SharedPointerComparator<ElementType::Pointer>>;

    using NodesCloudSetType = std::unordered_set<NodeType::Pointer, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    using CloudDataVectorType = DenseVector<std::pair<NodeType::Pointer, double>>;

    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ShiftedBoundaryPointBasedInterfaceUtility
    KRATOS_CLASS_POINTER_DEFINITION(ShiftedBoundaryPointBasedInterfaceUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Standard constructor
    /// @param rModel Model container
    /// @param ThisParameters Parameters object encapsulating the settings
    ShiftedBoundaryPointBasedInterfaceUtility(
        Model& rModel,
        Parameters ThisParameters);

    /// Copy constructor.
    ShiftedBoundaryPointBasedInterfaceUtility(ShiftedBoundaryPointBasedInterfaceUtility const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ShiftedBoundaryPointBasedInterfaceUtility& operator=(ShiftedBoundaryPointBasedInterfaceUtility const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    //TODO
    void CalculateAndAddPointBasedInterface();

    //TODO
    void ResetFlags();

    /** TODO
     * @brief Set BOUNDARY flags for elements that are intersected by the tessellated skin geometry.
       Calculate DISTANCE of nodes of intersected elements to the skin geometry and relocate those with a very small DISTANCE value in direction of the normal.
     * element BOUNDARY : elements in which a part of the surface is located (split/ intersected elements)
     */
    void SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes();

    //TODO
    // To be done after SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes() to locate skin points after node relocation in the updated volume part element geometries.
    void LocateSkinPoints();

    // TODO Set the corresponding flags at the interface elements
    // element INTERFACE : elements owning the surrogate boundary nodes adjacent to an deactivated BOUNDARY element
    // To be done after all necessary elements have been declared SBM_BOUNDARY (--> after SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes and LocateSkinPoints)
    void SetInterfaceFlags();

    // TODO node ACTIVE : nodes that belong to the elements to be assembled (all nodes as the interface is discontinuous)
    // element ACTIVE : elements which are not BOUNDARY (the ones to be assembled)
    // To be done after all necessary elements have been declared SBM_BOUNDARY (--> after SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes and LocateSkinPoints)
    void DeactivateElementsAndNodes();

    //TODO
    void CalculateAndAddSkinIntegrationPointConditions();

    //TODO
    void FixEnclosedVolumesPressure();

    //TODO
    // Calculate positive and negative side velocity and pressure at the nodes of the skin model part
    // Result is stored in POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE
    template <std::size_t TDim>
    void CalculatePressureAtSkinNodesTemplated();
    void CalculatePressureAtSkinNodes();
    template <std::size_t TDim>
    void CalculateVelocityAtSkinNodesTemplated();
    void CalculateVelocityAtSkinNodes();
    //void CalculateVelocityAtSkinNodes();
    //void CalculateDragForceAtSkinPoints();
    template <std::size_t TDim>
    void CalculateSkinDragTemplated();
    void CalculateSkinDrag();
    //     const Variable<array_1d<double, 3>>& rVariable,
    //     array_1d<double, 3>& rOutput);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    const Parameters GetDefaultParameters() const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "ShiftedBoundaryPointBasedInterfaceUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ShiftedBoundaryPointBasedInterfaceUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Static Member Variables
    ///@{

    //static constexpr std::size_t VoigtSize = 3 * (TDim-1);
    //static constexpr std::size_t BlockSize = TDim + 1;

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;
    //TODO needs to be discretized to be used for FindIntersectedGeometricalObjectsProcess.
    // Also right now integration points of mpSkinModelPart (1 GP) are taken as skin points with area and normal.
    ModelPart* mpSkinModelPart = nullptr;
    ModelPart* mpBoundarySubModelPart = nullptr;
    ModelPart* mpSkinPointsModelPart = nullptr;

    std::string mSkinModelPartName = "";

    SkinPointsToElementsMapType mSkinPointsMap;
    //TODO too big to store - (only) necessary for calculation of values at skin
    SidesVectorToElementsMapType mSidesVectorMap;
    NodesCloudMapType mExtensionOperatorMap;

    bool mConformingBasis;
    ExtensionOperator mExtensionOperator;
    std::size_t mMLSExtensionOperatorOrder;

    const Condition* mpConditionPrototype;

    bool mFindEnclosedVolumes = false;
    bool mPositiveSideIsActive = true;
    bool mNegativeSideIsActive = true;
    bool mPositiveSideIsEnclosed = false;
    bool mNegativeSideIsEnclosed = false;

    bool mCrossBoundaryNeighbors = false;
    bool mUseTessellatedBoundary = true;

    bool mAddedLeftEdge = false;
    bool mAddedRightEdge = false;

    /// @brief Protected empty constructor for derived classes
    ShiftedBoundaryPointBasedInterfaceUtility() {}

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Calculate the extension basis coefficients
     * This method calculates the extension operator coefficients in a meshless fashion
     * First, a cloud of points is created for each one of the "inner" nodes. This cloud of points is then used in the
     * calculation of the MLS approximants at these points. Then, the basis is made conformant by using the meshless
     * approximants to calculate the interpolation at each one of the point conditions location.
     */
    //void CalculateMeshlessBasedConformingExtensionBasis();

    // TODO
    double CalculateSkinDistanceToNode(
        const NodeType& rNode,
        const PointerVector<GeometricalObject>& rIntersectingObjects,
        PointDistanceFunctionType pPointDistanceFunction,
        IntersectionPlaneConstructorType pIntersectionPlaneConstructor,
        const double& DistanceThreshold,
        const double& ThresholdForSignedness);

    /**
     * @brief TODO. This method requires the skin to be stored in mpSkinModelPart as a Kratos model part with elements, integration points and area normals.
     *
     * @tparam TDim Working space dimension
     * @param rSkinPointsMap
     */
    template <std::size_t TDim>
    void MapSkinPointsToElements(SkinPointsToElementsMapType& rSkinPointsMap);

    /**TODO*/
    void SetSidesVectorsAndSkinNormalsForSplitElements(
        const SkinPointsToElementsMapType& rSkinPointsMap,
        SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap);

    /** TODO
     * @brief Set the Extension Operators For Split Element Nodes object selecting a support cloud for each node of each split element and using MLS shape functions
     *
     * @param rExtensionOperatorMap
     */
    void SetExtensionOperatorsForSplitElementNodes(
        const SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap,
        NodesCloudMapType& rExtensionOperatorMap);

    //TODO
    /**
     * @brief Set the support cloud for the same side of the boundary of which the first nodes nodes are given.
     * For given nodes on one side of the boundary, this function creates a cloud of nodes around them on the same side of the boundary.
     * These first nodes could be the positive nodes of a BOUNDARY element in order to get a support cloud and extension basis for one of the element's negative nodes and vice versa.
     * The support cloud created by this function is to be used for calculating the MLS-based extension operator.
     * @param rSameSideNodes Pointer to the first nodes of the support cloud
     * TODO
     * @param rCloudNodes Vector containing pointers to the nodes of the cloud
     * @param rCloudCoordinates Matrix containing the coordinates of the nodes of the cloud
     */
    void SetLateralSupportCloud(
        const NodeType::Pointer pOtherSideNode,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates,
        const Kratos::Flags& rSearchSideFlag);

    //TODO without dot product
    void AddLateralSupportLayer(
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesCloudSetType& SupportNodesSet);

    //TODO with dot product
    void AddLateralSupportLayer(
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesCloudSetType& SupportNodesSet);

    /**
     * @brief Create a pointer vector of pointers to all the nodes affecting the respective side of a split element's boundary.
     * Create vectors to pointers of cloud nodes for a split element using all the nodes affecting the element.
     * Pointer Vectors are sorted by ID to properly get the extension operator data
     *
     * @param rElement
     * @param rSidesVector
     * @param rExtensionOperatorMap
     * @param rCloudNodeVectorPositiveSide
     * @param rCloudNodeVectorNegativeSide
     */
    void CreateCloudNodeVectorsForSplitElement(
        const ElementType& rElement,
        const Vector& rSidesVector,
        NodesCloudMapType& rExtensionOperatorMap,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide);

    /* TODO */
    void GetDataForSplitElementSkinPoint(
        const ElementType& rElement,
        const array_1d<double,3>& rIntPtCoordinates,
        Vector& rIntPtShapeFunctionValues,
        Matrix& rIntPtShapeFunctionDerivatives);

    /* TODO */
    bool AddIntegrationPointCondition(
        const ElementType& rElement,
        const Vector& rSidesVector,
        const double ElementSize,
        const array_1d<double,3>& rIntPtCoordinates,
        const array_1d<double,3>& rIntPtAreaNormal,
        NodesCloudMapType& rExtensionOperatorMap,
        const PointerVector<NodeType>& rCloudNodeVector,
        const Vector& rIntPtShapeFunctionValues,
        const Matrix& rIntPtShapeFunctionDerivatives,
        std::size_t& r_ConditionId,
        const bool ConsiderPositiveSide);

    bool SetEnclosedNodesPressure(
        ElementType& rElement,
        const Vector& rSidesVector,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal);

    template <std::size_t TDim>
    bool CalculateUnknownsForBothSidesOfSplitElement(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns);

    //TODO
    // Calculate positive and negative side pressure inside a given SBM_BOUNDARY element using given shape function values.
    // returns true if pressure of point was calculated successfully
    bool CalculatePressureAtSplitElementSkinPoint(
        const ElementType::Pointer pElement,
        const Vector& rPointShapeFunctionValues,
        double& rPositiveSidePressure,
        double& rNegativeSidePressure);

    bool CalculateVelocityAtSplitElementSkinPoint(
        const ElementType::Pointer pElement,
        const Vector& rPointShapeFunctionValues,
        array_1d<double,3>& rPositiveSideVelocity,
        array_1d<double,3>& rNegativeSideVelocity);

    /**
     * @brief Get the MLS shape functions factory object
     * This function returns a prototype for the MLS shape functions calculation
     * @return MLSShapeFunctionsFunctionType MLS shape functions call prototype
     */
    MeshlessShapeFunctionsFunctionType GetMLSShapeFunctionsFunction() const;

    /**
     * @brief Get the RBF shape functions factory object
     * This function returns a prototype for the RBF shape functions calculation
     * @return RBFShapeFunctionsFunctionType RBF shape functions call prototype
     */
    MeshlessShapeFunctionsFunctionType GetRBFShapeFunctionsFunction() const;

    /**
     * @brief Get the element size function object
     * This function returns a prototype for the element size calculation call
     * @param rGeometry Geometry to calculate the element size
     * @return ElementSizeFunctionType Element size calculation call
     */
    ElementSizeFunctionType GetElementSizeFunction(const GeometryType& rGeometry);

    /**
     * @brief Calculates the kernel function radius
     * For the given cloud of points coordinates, this function calculates the maximum distance between
     * the center of the kernel (origin) and the points. This is supposed to be used in the MLS
     * approximation calculation.
     * @param rCloudCoordinates Matrix containing the coordinates of the nodes of the cloud
     * @param rOrigin Coordinates of the point on which the kernel function is centered
     * @return double Kernel function radius
     */
    double CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin);

    /**
     * @brief Get the number of required points
     * For the MLS approximation, this function returns the minimum number of required points according to
     * dimension and order. If the cloud of points has less points than the value returned by this function
     * the case is ill-conditioned, meaning that the cloud needs to be enlarged.
     * @return std::size_t Number of required points
     */
    std::size_t GetRequiredNumberOfPoints();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class ShiftedBoundaryPointBasedInterfaceUtility

}  // namespace Kratos.
