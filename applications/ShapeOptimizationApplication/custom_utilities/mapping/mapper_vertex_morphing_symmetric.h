// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_SYMMETRIC_H
#define MAPPER_VERTEX_MORPHING_SYMMETRIC_H

// System includes
#include <iostream>
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "custom_utilities/filter_function.h"
#include "symmetry_base.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) MapperVertexMorphingSymmetric : public Mapper
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef Node NodeType;
    typedef NodeType::Pointer NodePointerType;
    typedef std::vector<NodePointerType> NodeVectorType;
    typedef std::vector<NodePointerType>::iterator NodeIteratorType;
    typedef std::vector<double>::iterator DoubleVectorIteratorType;
    typedef array_1d<double,3> array_3d;

    // Type definitions for linear algebra including sparse systems
    typedef TUblasSparseSpace<double> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVectorType, NodePointerType, NodeIteratorType, DoubleVectorIteratorType > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of MapperVertexMorphingSymmetric
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingSymmetric);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingSymmetric( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : mrOriginModelPart(rOriginModelPart),
          mrDestinationModelPart(rDestinationModelPart),
          mMapperSettings(MapperSettings),
          mFilterRadius(MapperSettings["filter_radius"].GetDouble())
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingSymmetric()
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override;

    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override;

    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override;

    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override;

    void Update() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "MapperVertexMorphingSymmetric";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingSymmetric";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    // Initialized by class constructor
    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    Parameters mMapperSettings;
    Kratos::unique_ptr<FilterFunction> mpFilterFunction;
    bool mIsMappingInitialized = false;

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void InitializeComputationOfMappingMatrix();

    ///@}

private:
    ///@name Member Variables
    ///@{

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    Kratos::unique_ptr<KDTree> mpSearchTree;

    // Variables for mapping
    SparseMatrixType mMappingMatrix;
    SparseMatrixType mMappingMatrixScalar;
    double mFilterRadius;

    Kratos::unique_ptr<SymmetryBase> mpSymmetry;

    ///@}
    ///@name Private Operations
    ///@{

    void CreateFilterFunction();

    void InitializeMappingVariables();

    void AssignMappingIds();

    void CreateSearchTreeWithAllNodesInOriginModelPart();

    void ComputeMappingMatrix();

    void AllocateMatrix();

    virtual void ComputeWeightForAllNeighbors(  const ModelPart::NodeType& destination_node,
                                                const NodeVectorType& neighbor_nodes,
                                                const unsigned int number_of_neighbors,
                                                std::vector<double>& list_of_weights,
                                                double& sum_of_weights );

    void FillMappingMatrixWithWeights(  const ModelPart::NodeType& destination_node,
                                        const NodeVectorType& neighbor_nodes,
                                        const unsigned int number_of_neighbors,
                                        const std::vector<double>& list_of_weights,
                                        const std::vector<bool>& transform,
                                        const double& sum_of_weights );

    double GetVertexMorphingRadius(const NodeType& rNode) const override
    {
        return mFilterRadius;
    }

}; // Class MapperVertexMorphingSymmetric

}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_SYMMETRIC_H
