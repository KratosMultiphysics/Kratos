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
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/shifted_boundary_meshless_interface_utility.h"

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
 * This class provides the utilities for calculating a meshless conforming extension operator for a Shifted Boundary Method
 * without Taylor expansions allowing for a discontinuous interface, which is defined by a level set (e.g. ELEMENTAL_DISTANCES)
 * from the calculation of intersections between the volume mesh (e.g. fluid) and the mesh of the interface (e.g. thin-walled structure).
 */
class KRATOS_API(KRATOS_CORE) ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility : public ShiftedBoundaryMeshlessInterfaceUtility
{
public:

    ///@name Type Definitions
    ///@{

    enum class ExtensionOperator
    {
        MLS,
        RBF,
        GradientBased
    };

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using ShapeFunctionsGradientsType = GeometryType::ShapeFunctionsGradientsType;

    using ModifiedShapeFunctionsFactoryType = std::function<ModifiedShapeFunctions::UniquePointer(const GeometryType::Pointer, const Vector&)>;

    using MeshlessShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

    using MLSShapeFunctionsAndGradientsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&, Matrix&)>;

    using ElementSizeFunctionType = std::function<double(const GeometryType&)>;

    using NodesCloudSetType = std::unordered_set<NodeType::Pointer, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    using CloudDataVectorType = DenseVector<std::pair<NodeType::Pointer, double>>;

    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility
    KRATOS_CLASS_POINTER_DEFINITION(ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Standard constructor
    /// @param rModel Model container
    /// @param ThisParameters Parameters object encapsulating the settings
    ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility(
        Model& rModel,
        Parameters ThisParameters);

    /// Copy constructor.
    ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility(ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility& operator=(ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

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
        return "ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const Variable<Vector>* mpDiscontinuousLevelSetVariable;

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
    void CalculateMeshlessBasedConformingExtensionBasis() override;

    /**
     * @brief Set the corresponding flags at the interface
     * This methods sets the corresponding flags at the interface nodes and elements. These are
     * node ACTIVE : nodes that belong to the elements to be assembled (all nodes as the interface is discontinuous)
     * node INTERFACE : nodes that belong to the cloud of points (multiple layers surrounding the BOUNDARY elements)
     * element ACTIVE : elements which are not intersected by the boundary (the ones to be assembled)
     * element BOUNDARY : intersected elements, as well as incised elements if ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED were calculated
     * element INTERFACE : elements owning the surrogate boundary nodes adjacent to an deactivated BOUNDARY element
     * Incised elements will be ACTIVE, BOUNDARY and INTERFACE at the same time, to be assembled and contributing to the boundary flux as if they weren't incised,
     * but still part of BOUNDARY, so that the point cloud for one side can not cross over to the other side via an incised element, separating both sides better.
     */
    void SetInterfaceFlags() override;

    /** TODO
     * @brief Set the Extension Operators For Split Element Nodes object selecting a support cloud for each node of each split element and using MLS shape functions
     *
     * @param rExtensionOperatorMap
     */
    void SetExtensionOperatorsForSplitElementNodes(
        NodesCloudMapType& rExtensionOperatorMap
    );

    /**
     * @brief Set the support cloud for the same side of the boundary of which the first nodes nodes are given.
     * For given nodes on one side of the boundary, this function creates a cloud of nodes around them on the same side of the boundary.
     * These first nodes could be the positive nodes of a BOUNDARY element in order to get a support cloud and extension basis for one of the element's negative nodes and vice versa.
     * The support cloud created by this function is to be used for calculating the MLS-based extension operator.
     * @param rSameSideNodes Pointer to the first nodes of the support cloud
     * @param rCloudNodes Vector containing pointers to the nodes of the cloud
     * @param rCloudCoordinates Matrix containing the coordinates of the nodes of the cloud
     */
    void SetLateralSupportCloud(
        const std::vector<NodeType::Pointer>& rSameSideNodes,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates);

    /** TODO
     * @brief Get the Data For Split Element Boundary object using modified shape functions
     *
     * @param rElement
     * @param pModifiedShapeFunctionsFactory
     * @param rBoundaryShapeFunctionValues
     * @param rBoundaryShapeFunctionDerivatives
     * @param rBoundaryWeights
     * @param rBoundaryAreaNormals
     */
    void GetDataForSplitElementBoundary(
        const ElementType& rElement,
        ModifiedShapeFunctionsFactoryType pModifiedShapeFunctionsFactory,
        Matrix& rBoundaryShapeFunctionValues,
        ModifiedShapeFunctions::ShapeFunctionsGradientsType& rBoundaryShapeFunctionDerivatives,
        Vector& rBoundaryWeights,
        std::vector<array_1d<double,3>>& rBoundaryAreaNormals);

    /**
     * @brief Create a pointer vector of pointers to all the nodes affecting the respective side of a split element's boundary.
     * Create vectors to pointers of cloud nodes for a split element using all the nodes affecting the element.
     * Pointer Vectors are sorted by ID to properly get the extension operator data
     *
     * @param rElement
     * @param rExtensionOperatorMap
     * @param rCloudNodeVectorPositiveSide
     * @param rCloudNodeVectorNegativeSide
     */
    void CreateCloudNodeVectorsForSplitElement(
        const ElementType& rElement,
        NodesCloudMapType& rExtensionOperatorMap,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide);

    /**
     * @brief
     *
     * @param rElement
     * @param ElementSize
     * @param rExtensionOperatorMap
     * @param rCloudNodeVector
     * @param rIntPtCoordinates
     * @param rIntPtShapeFunctionValues
     * @param rIntPtShapeFunctionDerivatives
     * @param rIntPtWeight
     * @param rIntPtNormal
     * @param conditionId
     * @param ConsiderPositiveSide
     */
    void AddIntegrationPointCondition(
        const ElementType& rElement,
        const double ElementSize,
        NodesCloudMapType& rExtensionOperatorMap,
        const PointerVector<NodeType>& rCloudNodeVector,
        const array_1d<double,3>& rIntPtCoordinates,
        const DenseVector<double>& rIntPtShapeFunctionValues,
        const DenseMatrix<double>& rIntPtShapeFunctionDerivatives,
        const double rIntPtWeight,
        const array_1d<double,3>& rIntPtNormal,
        const std::size_t conditionId,
        bool ConsiderPositiveSide);

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
}; // Class ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility

} // namespace Kratos.
