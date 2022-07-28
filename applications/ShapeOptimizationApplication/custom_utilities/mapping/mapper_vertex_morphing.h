// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_H
#define MAPPER_VERTEX_MORPHING_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "custom_utilities/filter_function.h"

// ==============================================================================

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

/// Short class definition.
/** Detail class definition.
*/

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) MapperVertexMorphing : public Mapper
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef array_1d<double,3> array_3d;

    // Type definitions for linear algebra including sparse systems
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of MapperVertexMorphing
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphing( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : mrOriginModelPart(rOriginModelPart),
          mrDestinationModelPart(rDestinationModelPart),
          mMapperSettings(MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphing()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override;

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override;

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override;

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override;

    // --------------------------------------------------------------------------
    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override;

    // --------------------------------------------------------------------------
    void Update() override;

    // --------------------------------------------------------------------------

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
    virtual std::string Info() const override
    {
        return "MapperVertexMorphing";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphing";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
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

    // Initialized by class constructor
    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    Parameters mMapperSettings;
    FilterFunction::UniquePointer mpFilterFunction;
    bool mIsMappingInitialized = false;
    SparseMatrixType mMappingMatrix;
    KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual void InitializeComputationOfMappingMatrix();

    // --------------------------------------------------------------------------
    virtual void ComputeMappingMatrix()
    {
        InitializeComputationOfMappingMatrix();
        CreateSearchTreeWithAllNodesInOriginModelPart();

        double filter_radius = mMapperSettings["filter_radius"].GetDouble();
        unsigned int max_number_of_neighbors = mMapperSettings["max_nodes_in_filter_radius"].GetInt();

        // working vectors
        std::vector<double> resulting_squared_distances( max_number_of_neighbors );
        NodeVector neighbor_nodes( max_number_of_neighbors );
        std::vector<double> list_of_weights( max_number_of_neighbors );

        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             filter_radius,
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(),
                                                                             max_number_of_neighbors );

            list_of_weights.resize(number_of_neighbors);
            std::fill(list_of_weights.begin(), list_of_weights.end(), 0.0);
            double sum_of_weights = 0.0;

            if(number_of_neighbors >= max_number_of_neighbors)
                KRATOS_WARNING("ShapeOpt::MapperVertexMorphing") << "For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << max_number_of_neighbors << " nodes) reached!" << std::endl;

            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            FillMappingMatrixWithWeights( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
        }
    }

    // --------------------------------------------------------------------------
    virtual void ComputeWeightForAllNeighbors(  ModelPart::NodeType& origin_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double weight = mpFilterFunction->ComputeWeight( origin_node.Coordinates(), neighbor_node.Coordinates() );

            list_of_weights[neighbor_itr] = weight;
            sum_of_weights += weight;
        }
    }

    // --------------------------------------------------------------------------
    virtual void InitializeMappingVariables()
    {
        const unsigned int origin_node_number = mrOriginModelPart.Nodes().size();
        mValuesOrigin.resize(3);
        mValuesOrigin[0] = ZeroVector(origin_node_number);
        mValuesOrigin[1] = ZeroVector(origin_node_number);
        mValuesOrigin[2] = ZeroVector(origin_node_number);

        const unsigned int destination_node_number = mrDestinationModelPart.Nodes().size();
        mValuesDestination.resize(3);
        mValuesDestination[0] = ZeroVector(destination_node_number);
        mValuesDestination[1] = ZeroVector(destination_node_number);
        mValuesDestination[2] = ZeroVector(destination_node_number);

        mMappingMatrix.resize(destination_node_number,origin_node_number,false);
    }

    // --------------------------------------------------------------------------
    virtual void FillMappingMatrixWithWeights(  ModelPart::NodeType& destination_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        unsigned int row_id = destination_node.GetValue(MAPPING_ID);

        std::vector<std::pair<int, double>> sorted_neighbors;
        sorted_neighbors.reserve(number_of_neighbors);

        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            sorted_neighbors.push_back(std::make_pair(neighbor_node.GetValue(MAPPING_ID), list_of_weights[neighbor_itr]));
        }

        // Sort the vector of pairs according the column id - speeds up matrix insertion
        std::sort(std::begin(sorted_neighbors), std::end(sorted_neighbors),
            [&](const std::pair<int, double>& a, const std::pair<int, double>& b)
            {
                return a.first < b.first;
            });

        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            int collumn_id = sorted_neighbors[neighbor_itr].first;
            double weight = sorted_neighbors[neighbor_itr].second / sum_of_weights;
            mMappingMatrix.insert_element(row_id,collumn_id,weight);
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

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesInOriginModelPart;

    // Variables for mapping
    std::vector<Vector> mValuesOrigin;
    std::vector<Vector> mValuesDestination;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CreateListOfNodesInOriginModelPart();

    // --------------------------------------------------------------------------
    void CreateFilterFunction();

    // --------------------------------------------------------------------------
    void InitializeMappingVariables();

    // --------------------------------------------------------------------------
    void AssignMappingIds();

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesInOriginModelPart();

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix();

    // --------------------------------------------------------------------------
    virtual void ComputeWeightForAllNeighbors(  ModelPart::NodeType& origin_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights );

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& origin_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights );

    // --------------------------------------------------------------------------

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      MapperVertexMorphing& operator=(MapperVertexMorphing const& rOther);

    /// Copy constructor.
//      MapperVertexMorphing(MapperVertexMorphing const& rOther);


    ///@}

}; // Class MapperVertexMorphing

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_H
