// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "custom_utilities/filter_function.h"
#include "symmetry_base.h"
#include "symmetry_plane.h"
#include "symmetry_revolution.h"
#include "mapper_vertex_morphing_symmetric.h"

// ==============================================================================

namespace Kratos
{

void MapperVertexMorphingSymmetric::Initialize()
{
    BuiltinTimer timer;
    KRATOS_INFO("ShapeOpt") << "Starting initialization of mapper..." << std::endl;

    CreateFilterFunction();

    mIsMappingInitialized = true;

    Update();

    KRATOS_INFO("ShapeOpt") << "Finished initialization of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphingSymmetric::Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

    Vector values_origin(mrOriginModelPart.Nodes().size()*3);
    Vector values_destination(mrDestinationModelPart.Nodes().size()*3);

    // Prepare vectors for mapping
    values_origin.clear();
    values_destination.clear();

    block_for_each(mrOriginModelPart.Nodes(), [&](const ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        const array_3d& r_nodal_variable = rNode.FastGetSolutionStepValue(rOriginVariable);
        values_origin[i*3+0] = r_nodal_variable[0];
        values_origin[i*3+1] = r_nodal_variable[1];
        values_origin[i*3+2] = r_nodal_variable[2];
    });

    // Perform mapping
    SparseSpaceType::Mult(mMappingMatrix, values_origin, values_destination);

    // Assign results to nodal variable
    block_for_each(mrDestinationModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        array_3d& r_node_vector = rNode.FastGetSolutionStepValue(rDestinationVariable);
        r_node_vector(0) = values_destination[i*3+0];
        r_node_vector(1) = values_destination[i*3+1];
        r_node_vector(2) = values_destination[i*3+2];
    });

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphingSymmetric::Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

    Vector values_origin(mrOriginModelPart.Nodes().size());
    Vector values_destination(mrDestinationModelPart.Nodes().size());

    // Prepare vectors for mapping
    values_origin.clear();
    values_destination.clear();

    block_for_each(mrOriginModelPart.Nodes(), [&](const ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        const double& r_nodal_variable = rNode.FastGetSolutionStepValue(rOriginVariable);
        values_origin[i] = r_nodal_variable;
    });

    // Perform mapping
    SparseSpaceType::Mult(mMappingMatrixScalar, values_origin, values_destination);

    // Assign results to nodal variable
    block_for_each(mrDestinationModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        double& r_node_double = rNode.FastGetSolutionStepValue(rDestinationVariable);
        r_node_double = values_destination[i];
    });

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;

 }

void MapperVertexMorphingSymmetric::InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("ShapeOpt") << "\nStarting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;

    Vector values_origin(mrOriginModelPart.Nodes().size()*3);
    Vector values_destination(mrDestinationModelPart.Nodes().size()*3);

    // Prepare vectors for mapping
    values_origin.clear();
    values_destination.clear();

    block_for_each(mrDestinationModelPart.Nodes(), [&](const ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        const array_3d& r_nodal_variable = rNode.FastGetSolutionStepValue(rDestinationVariable);
        values_destination[i*3+0] = r_nodal_variable[0];
        values_destination[i*3+1] = r_nodal_variable[1];
        values_destination[i*3+2] = r_nodal_variable[2];
    });

    SparseSpaceType::TransposeMult(mMappingMatrix,values_destination,values_origin);

    // Assign results to nodal variable
    block_for_each(mrOriginModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        array_3d& r_node_vector = rNode.FastGetSolutionStepValue(rOriginVariable);
        r_node_vector(0) = values_origin[i*3+0];
        r_node_vector(1) = values_origin[i*3+1];
        r_node_vector(2) = values_origin[i*3+2];
    });

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphingSymmetric::InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

    Vector values_origin(mrOriginModelPart.Nodes().size());
    Vector values_destination(mrDestinationModelPart.Nodes().size());

    // Prepare vectors for mapping
    values_origin.clear();
    values_destination.clear();

    block_for_each(mrDestinationModelPart.Nodes(), [&](const ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        const double& r_nodal_value = rNode.FastGetSolutionStepValue(rDestinationVariable);
        values_destination[i] = r_nodal_value;
    });

    SparseSpaceType::TransposeMult(mMappingMatrixScalar,values_destination,values_origin);

    // Assign results to nodal variable
    block_for_each(mrOriginModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        const int i = rNode.GetValue(MAPPING_ID);
        double& r_nodal_value = rNode.FastGetSolutionStepValue(rOriginVariable);
        r_nodal_value = values_origin[i];
    });

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphingSymmetric::Update()
{
    if (mIsMappingInitialized == false) {
        KRATOS_ERROR << "Mapping has to be initialized before calling the Update-function!";
    }

    BuiltinTimer timer;
    KRATOS_INFO("ShapeOpt") << "Starting to update mapper..." << std::endl;

    InitializeMappingVariables();
    AssignMappingIds();

    if (mMapperSettings["plane_symmetry"].GetBool()) {
        mpSymmetry = Kratos::make_unique<SymmetryPlane>(mrOriginModelPart, mrDestinationModelPart, mMapperSettings["plane_symmetry_settings"]);
    } else if (mMapperSettings["revolution"].GetBool()) {
        mpSymmetry = Kratos::make_unique<SymmetryRevolution>(mrOriginModelPart, mrDestinationModelPart, mMapperSettings["revolution_settings"]);
    } else {
        KRATOS_ERROR << "No symmetry type specified" << std::endl;
    }

    ComputeMappingMatrix();

    KRATOS_INFO("ShapeOpt") << "Finished updating of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphingSymmetric::InitializeComputationOfMappingMatrix()
{
    mpSearchTree.reset();
    mMappingMatrix.clear();
    mMappingMatrixScalar.clear();
}

void MapperVertexMorphingSymmetric::CreateFilterFunction()
{
    const std::string filter_type = mMapperSettings["filter_function_type"].GetString();

    mpFilterFunction = Kratos::make_unique<FilterFunction>(filter_type);
}

void MapperVertexMorphingSymmetric::InitializeMappingVariables()
{
    const unsigned int origin_node_number = mrOriginModelPart.Nodes().size();
    const unsigned int destination_node_number = mrDestinationModelPart.Nodes().size();
    mMappingMatrix.resize(destination_node_number*3,origin_node_number*3,false);
    mMappingMatrixScalar.resize(destination_node_number,origin_node_number,false);
}

void MapperVertexMorphingSymmetric::AssignMappingIds()
{
    // Note: loop in the same order as in AllocateMatrix(), to avoid reallocations of the matrix.
    IndexPartition<int>(mrOriginModelPart.NumberOfNodes()).for_each([&](const int Index) {
        (mrOriginModelPart.NodesBegin() + Index)->SetValue(MAPPING_ID, Index);
    });

    IndexPartition<int>(mrDestinationModelPart.NumberOfNodes()).for_each([&](const int Index) {
        (mrDestinationModelPart.NodesBegin() + Index)->SetValue(MAPPING_ID, Index);
    });
}

void MapperVertexMorphingSymmetric::CreateSearchTreeWithAllNodesInOriginModelPart()
{
    mpSearchTree = Kratos::make_unique<KDTree>(mpSymmetry->GetOriginSearchNodes().begin(), mpSymmetry->GetOriginSearchNodes().end(), mBucketSize);
}

void MapperVertexMorphingSymmetric::ComputeMappingMatrix()
{
    InitializeComputationOfMappingMatrix();
    CreateSearchTreeWithAllNodesInOriginModelPart();
    AllocateMatrix();
    const double filter_radius = mMapperSettings["filter_radius"].GetDouble();
    const unsigned int max_number_of_neighbors = mMapperSettings["max_nodes_in_filter_radius"].GetInt();

    struct tls_vecs {
        // variable size, fixed capacity vecs
        std::vector<bool> transform;
        NodeVectorType total_neighbor_nodes;
        std::vector<double> total_list_of_weights;
        std::vector<double> list_of_weights;

        // fixed size vec
        NodeVectorType neighbor_nodes;

        // constructor
        explicit tls_vecs(const unsigned int MaxNumberOfNeighbors) {
            transform.reserve(MaxNumberOfNeighbors);
            total_neighbor_nodes.reserve( MaxNumberOfNeighbors );
            total_list_of_weights.reserve( MaxNumberOfNeighbors );
            list_of_weights.reserve( MaxNumberOfNeighbors );
            neighbor_nodes.resize( MaxNumberOfNeighbors );
        }
    };

    block_for_each(mrDestinationModelPart.Nodes(), tls_vecs(max_number_of_neighbors), [&](const NodeType& rNode_i, tls_vecs& rTLS){
        rTLS.transform.clear();
        rTLS.total_neighbor_nodes.clear();
        rTLS.total_list_of_weights.clear();

        auto search_nodes = mpSymmetry->GetDestinationSearchNodes(rNode_i.GetValue(MAPPING_ID));
        unsigned int total_number_of_neighbors = 0;
        double total_sum_of_weights = 0;
        for (auto& pair_i : search_nodes) {
            const NodeType search_node(0, pair_i.first);
            const unsigned int number_of_neighbors = mpSearchTree->SearchInRadius(search_node,
                                                                            filter_radius,
                                                                            rTLS.neighbor_nodes.begin(),
                                                                            max_number_of_neighbors );

            total_number_of_neighbors += number_of_neighbors;
            rTLS.transform.resize(total_number_of_neighbors, pair_i.second);

            rTLS.list_of_weights.clear();
            rTLS.list_of_weights.resize( number_of_neighbors, 0.0 );
            ComputeWeightForAllNeighbors( search_node, rTLS.neighbor_nodes, number_of_neighbors, rTLS.list_of_weights, total_sum_of_weights );

            rTLS.total_list_of_weights.insert(rTLS.total_list_of_weights.end(), rTLS.list_of_weights.begin(), rTLS.list_of_weights.begin()+number_of_neighbors);
            rTLS.total_neighbor_nodes.insert(rTLS.total_neighbor_nodes.end(), rTLS.neighbor_nodes.begin(), rTLS.neighbor_nodes.begin()+number_of_neighbors);
        }

        if(total_number_of_neighbors >= max_number_of_neighbors)
            KRATOS_WARNING("ShapeOpt::MapperVertexMorphingSymmetric") << "For node " << rNode_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << max_number_of_neighbors << " nodes) reached!" << std::endl;

        // each thread assembles to one row only (and only once) thats why no explicit synchronisation is needed in FillMappingMatrixWithWeights
        FillMappingMatrixWithWeights( rNode_i, rTLS.total_neighbor_nodes, total_number_of_neighbors, rTLS.total_list_of_weights, rTLS.transform, total_sum_of_weights );
    });
}

void MapperVertexMorphingSymmetric::AllocateMatrix() {
    const double filter_radius = mMapperSettings["filter_radius"].GetDouble();
    const unsigned int max_number_of_neighbors = mMapperSettings["max_nodes_in_filter_radius"].GetInt();

    // fixed size vecs
    NodeVectorType neighbor_nodes( max_number_of_neighbors );
    for(auto& node_i : mrDestinationModelPart.Nodes())
    {
        auto search_nodes = mpSymmetry->GetDestinationSearchNodes(node_i.GetValue(MAPPING_ID));
        unsigned int total_number_of_neighbors = 0;
        for (auto& pair_i : search_nodes) {
            NodeType search_node(0, pair_i.first);
            const unsigned int number_of_neighbors = mpSearchTree->SearchInRadius(search_node,
                                                                            filter_radius,
                                                                            neighbor_nodes.begin() + total_number_of_neighbors,
                                                                            max_number_of_neighbors - total_number_of_neighbors);
            total_number_of_neighbors += number_of_neighbors;
        }

        if(total_number_of_neighbors >= max_number_of_neighbors)
            KRATOS_WARNING("ShapeOpt::MapperVertexMorphingSymmetric") << "For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << max_number_of_neighbors << " nodes) reached!" << std::endl;

        // allocate the mapping matrix row-by row, with sorted column entries
        const unsigned int row_id = node_i.GetValue(MAPPING_ID);

        std::vector<unsigned int> unique_col_ids(total_number_of_neighbors);
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<total_number_of_neighbors ; neighbor_itr++)
        {
            const NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            unique_col_ids[neighbor_itr] = neighbor_node.GetValue(MAPPING_ID);
        }
        std::sort (unique_col_ids.begin(), unique_col_ids.end());
        unique_col_ids.erase( unique( unique_col_ids.begin(), unique_col_ids.end() ), unique_col_ids.end() );

        // each node affects a 3x3 part of the mapping matrix
        for (unsigned int i: unique_col_ids) {
            mMappingMatrix.insert_element(row_id*3+0, i*3+0, 0.0);
            mMappingMatrix.insert_element(row_id*3+0, i*3+1, 0.0);
            mMappingMatrix.insert_element(row_id*3+0, i*3+2, 0.0);
        }
        for (unsigned int i: unique_col_ids) {
            mMappingMatrix.insert_element(row_id*3+1, i*3+0, 0.0);
            mMappingMatrix.insert_element(row_id*3+1, i*3+1, 0.0);
            mMappingMatrix.insert_element(row_id*3+1, i*3+2, 0.0);
        }
        for (unsigned int i: unique_col_ids) {
            mMappingMatrix.insert_element(row_id*3+2, i*3+0, 0.0);
            mMappingMatrix.insert_element(row_id*3+2, i*3+1, 0.0);
            mMappingMatrix.insert_element(row_id*3+2, i*3+2, 0.0);
        }
        // each node affects a 1x1 part of the scalar mapping matrix
        for (unsigned int i: unique_col_ids) {
            mMappingMatrixScalar.insert_element(row_id, i, 0.0);
        }
    }
}

void MapperVertexMorphingSymmetric::ComputeWeightForAllNeighbors(
    const ModelPart::NodeType& destination_node,
    const NodeVectorType& neighbor_nodes,
    const unsigned int number_of_neighbors,
    std::vector<double>& list_of_weights,
    double& sum_of_weights )
{
    for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
    {
        const NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
        const double weight = mpFilterFunction->ComputeWeight( destination_node.Coordinates(), neighbor_node.Coordinates(), GetVertexMorphingRadius(destination_node) );

        list_of_weights[neighbor_itr] = weight;
        sum_of_weights += weight;
    }
}

void MapperVertexMorphingSymmetric::FillMappingMatrixWithWeights(
    const ModelPart::NodeType& destination_node,
    const NodeVectorType& neighbor_nodes,
    const unsigned int number_of_neighbors,
    const std::vector<double>& list_of_weights,
    const std::vector<bool>& transform,
    const double& sum_of_weights )
{
    const unsigned int row_id = destination_node.GetValue(MAPPING_ID);

    BoundedMatrix<double, 3, 3> transformation_matrix;
    for(unsigned int neighbor_itr = 0; neighbor_itr<number_of_neighbors; neighbor_itr++)
    {
        const NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
        const unsigned int column_id = neighbor_node.GetValue(MAPPING_ID);

        if (transform[neighbor_itr]) {
            mpSymmetry->TransformationMatrix(row_id, column_id, transformation_matrix);
        } else {
            noalias(transformation_matrix) = IdentityMatrix(3);
        }
        const double weight = list_of_weights[neighbor_itr] / sum_of_weights;
        for (unsigned int i=0; i<3; ++i){
            for (unsigned int j=0; j<3; ++j) {
                mMappingMatrix(row_id*3+i, column_id*3+j) += weight*transformation_matrix(i,j);
            }
        }
        // no transformation of weights for scalar mapping necessary
        mMappingMatrixScalar(row_id, column_id) += weight;
    }
}

}  // namespace Kratos.
