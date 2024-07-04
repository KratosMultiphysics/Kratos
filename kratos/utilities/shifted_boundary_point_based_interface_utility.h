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
        MLS
    };

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using ElementType = Element;

    using GeometryType = ModelPart::GeometryType;

    using ShapeFunctionsGradientsType = GeometryType::ShapeFunctionsGradientsType;

    using MeshlessShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

    using ElementSizeFunctionType = std::function<double(const GeometryType&)>;

    using SkinPointsDataVectorType = DenseVector<std::pair< array_1d<double,3>, array_1d<double,3> >>; // vector of position and area normal of skin points

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

    void AddSkinIntegrationPointConditions();

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

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;
    ModelPart* mpSkinModelPart = nullptr;
    ModelPart* mpBoundarySubModelPart = nullptr;

    bool mConformingBasis;
    ExtensionOperator mExtensionOperator;
    std::size_t mMLSExtensionOperatorOrder;

    const Condition* mpConditionPrototype;

    bool mInterpolateBoundary;

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
    void CalculateMeshlessBasedConformingExtensionBasis();

    /**
     * @brief TODO
     *
     * @tparam TDim Working space dimension
     * @param rSkinPointsMap
     */
    template <std::size_t TDim>
    void MapSkinPointsToElements(SkinPointsToElementsMapType& rSkinPointsMap);

    /**
     * @brief Set the corresponding flags at the interface
     * This methods sets the corresponding flags at the interface nodes and elements. These are
     * node ACTIVE : nodes that belong to the elements to be assembled (all nodes as the interface is discontinuous)
     * node INTERFACE : nodes that belong to the cloud of points (multiple layers surrounding the BOUNDARY elements)
     * element ACTIVE : elements which are not BOUNDARY (the ones to be assembled)
     * element BOUNDARY : elements in which a part of the surface is located (split elements)
     * element INTERFACE : elements owning the surrogate boundary nodes adjacent to an deactivated BOUNDARY element
     */
    void SetInterfaceFlags(const SkinPointsToElementsMapType& rSkinPointsMap);

    /**TODO*/
    void SetSidesVectorsAndSkinNormalsForSplitElements(
        const SkinPointsToElementsMapType& rSkinPointsMap,
        SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap
    );

    /** TODO
     * @brief Set the Extension Operators For Split Element Nodes object selecting a support cloud for each node of each split element and using MLS shape functions
     *
     * @param rExtensionOperatorMap
     */
    void SetExtensionOperatorsForSplitElementNodes(
        const SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap,
        NodesCloudMapType& rExtensionOperatorMap
    );

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
        const std::vector<NodeType::Pointer>& rSameSideNodes,
        const array_1d<double,3>& rSkinPosition,
        const array_1d<double,3>& rSkinNormal,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates,
        bool ConsiderPositiveSide);

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
    void GetDataForSplitElementIntegrationPoint(
        const ElementType& rElement,
        const array_1d<double,3>& rIntPtCoordinates,
        Vector& rIntPtShapeFunctionValues,
        Matrix& rIntPtShapeFunctionDerivatives);

    /* TODO */
    void AddIntegrationPointCondition(
        const ElementType& rElement,
        const Vector& rSidesVector,
        const double ElementSize,
        const array_1d<double,3>& rIntPtCoordinates,
        const array_1d<double,3>& rIntPtAreaNormal,
        NodesCloudMapType& rExtensionOperatorMap,
        const PointerVector<NodeType>& rCloudNodeVector,
        const Vector& rIntPtShapeFunctionValues,
        const Matrix& rIntPtShapeFunctionDerivatives,
        const std::size_t ConditionId,
        bool ConsiderPositiveSide);

    /**
     * @brief Get the MLS shape functions factory object
     * This function returns a prototype for the MLS shape functions calculation
     * @return MLSShapeFunctionsFunctionType MLS shape functions call prototype
     */
    MeshlessShapeFunctionsFunctionType GetMLSShapeFunctionsFunction() const;

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
     * the center of the kernel (origin) and the points. This is supposed to be used in the meshless
     * approximants calculation.
     * @param rCloudCoordinates Matrix containing the coordinates of the nodes of the cloud
     * @param rOrigin Coordinates of the point on which the kernel function is centered
     * @return double Kernel function radius
     */
    double CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin);

    /**
     * @brief Get the number of required points
     * For the MLS approximants, this function returns the minimum number of required points according to
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

} // namespace Kratos.
