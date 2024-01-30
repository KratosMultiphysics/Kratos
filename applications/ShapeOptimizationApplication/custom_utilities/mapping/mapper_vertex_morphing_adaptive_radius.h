// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    David Schmölz
//                   Suneth Warnakulasuriya
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_ADAPTIVE_RADIUS_H
#define MAPPER_VERTEX_MORPHING_ADAPTIVE_RADIUS_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

*/

template<class TBaseVertexMorphingMapper>
class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) MapperVertexMorphingAdaptiveRadius : public TBaseVertexMorphingMapper
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = TBaseVertexMorphingMapper;

    // Type definitions for better reading later
    using NodeType = Node;

    using ElementType = Element;

    using ConditionType = Element;

    using IndexType = std::size_t;

    using NodeTypePointer = NodeType::Pointer;

    using NodeVector = std::vector<NodeTypePointer>;

    using DoubleVectorIterator = std::vector<double>::iterator ;

    using NodeIterator = std::vector<NodeType::Pointer>::iterator;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;

    // Type definitions for tree-search
    using BucketType = Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator >;
    using KDTree = Tree< KDTreePartition<BucketType> >;

    /// Pointer definition of MapperVertexMorphingAdaptiveRadius
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingAdaptiveRadius);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingAdaptiveRadius(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Parameters MapperSettings);

    /// Destructor.
    virtual ~MapperVertexMorphingAdaptiveRadius() = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    void Update() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    // Initialized by class constructor
    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    std::string mRadiusFunctionType;
    double mRadiusFunctionParameter;
    double mMinimumFilterRadius;
    double mCurvatureLimit;
    IndexType mNumberOfSmoothingIterations;
    IndexType mMaxNumberOfNeighbors;

    IndexType mBucketSize = 100;
    Kratos::unique_ptr<KDTree> mpSearchTree;
    NodeVector mListOfNodesInOriginModelPart;

    ///@}
    ///@name Private Operations
    ///@{

    double CurvatureFunction(const double& rCurvature, const double& rElementSize);

    void CalculateCurvatureBasedFilterRadius();

    void SmoothenCurvatureBasedFilterRadius();

    void CalculateAdaptiveVertexMorphingRadius();

    double GetVertexMorphingRadius(const NodeType& rNode) const override;

    void CreateSearchTreeWithAllNodesInOriginModelPart();

    void CreateListOfNodesInOriginModelPart();

    virtual void ComputeWeightForAllNeighbors(
        const ModelPart::NodeType& destination_node,
        const NodeVector& neighbor_nodes,
        const unsigned int number_of_neighbors,
        std::vector<double>& list_of_weights,
        double& sum_of_weights) override;

    void AssignMappingIds();

    ///@}

}; // Class MapperVertexMorphingAdaptiveRadius

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_ADAPTIVE_RADIUS_H
