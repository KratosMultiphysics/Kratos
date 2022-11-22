// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
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
#include "mapper_vertex_morphing.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/builtin_timer.h"
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "custom_utilities/filter_function.h"

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include "external_libraries/igl/include/igl/gaussian_curvature.h"
#include "external_libraries/igl/include/igl/massmatrix.h"
#include "external_libraries/igl/include/igl/invert_diag.h"
#include "external_libraries/igl/include/igl/writeOFF.h"

// ==============================================================================

namespace Kratos
{

void MapperVertexMorphing::Initialize()
{
    BuiltinTimer timer;
    KRATOS_INFO("ShapeOpt") << "Starting initialization of mapper..." << std::endl;

    CreateFilterFunction();
    mIsMappingInitialized = true;

    Update();

    KRATOS_INFO("ShapeOpt") << "Finished initialization of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphing::Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

    // Prepare vectors for mapping
    mValuesOrigin[0].clear();
    mValuesOrigin[1].clear();
    mValuesOrigin[2].clear();
    mValuesDestination[0].clear();
    mValuesDestination[1].clear();
    mValuesDestination[2].clear();

    for(auto& node_i : mrOriginModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);
        array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rOriginVariable);
        mValuesOrigin[0][i] = r_nodal_variable[0];
        mValuesOrigin[1][i] = r_nodal_variable[1];
        mValuesOrigin[2][i] = r_nodal_variable[2];
    }

    // Perform mapping
    noalias(mValuesDestination[0]) = prod(mMappingMatrix,mValuesOrigin[0]);
    noalias(mValuesDestination[1]) = prod(mMappingMatrix,mValuesOrigin[1]);
    noalias(mValuesDestination[2]) = prod(mMappingMatrix,mValuesOrigin[2]);

    // Assign results to nodal variable
    for(auto& node_i : mrDestinationModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);

        array_3d& r_node_vector = node_i.FastGetSolutionStepValue(rDestinationVariable);
        r_node_vector(0) = mValuesDestination[0][i];
        r_node_vector(1) = mValuesDestination[1][i];
        r_node_vector(2) = mValuesDestination[2][i];
    }

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphing::Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

    // Prepare vectors for mapping
    mValuesOrigin[0].clear();
    mValuesDestination[0].clear();

    for(auto& node_i : mrOriginModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);
        mValuesOrigin[0][i] = node_i.FastGetSolutionStepValue(rOriginVariable);
    }

    // Perform mapping
    noalias(mValuesDestination[0]) = prod(mMappingMatrix,mValuesOrigin[0]);

    // Assign results to nodal variable
    for(auto& node_i : mrDestinationModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);
        node_i.FastGetSolutionStepValue(rDestinationVariable) = mValuesDestination[0][i];
    }

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphing::InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;

    // Prepare vectors for mapping
    mValuesOrigin[0].clear();
    mValuesOrigin[1].clear();
    mValuesOrigin[2].clear();
    mValuesDestination[0].clear();
    mValuesDestination[1].clear();
    mValuesDestination[2].clear();

    for(auto& node_i : mrDestinationModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);
        array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rDestinationVariable);
        mValuesDestination[0][i] = r_nodal_variable[0];
        mValuesDestination[1][i] = r_nodal_variable[1];
        mValuesDestination[2][i] = r_nodal_variable[2];
    }

    // Perform mapping
    if(mMapperSettings["consistent_mapping"].GetBool())
    {
        KRATOS_ERROR_IF(mrOriginModelPart.Nodes().size() != mrDestinationModelPart.Nodes().size()) << "Consistent mapping requires matching origin and destination model part.";

        noalias(mValuesOrigin[0]) = prod(mMappingMatrix,mValuesDestination[0]);
        noalias(mValuesOrigin[1]) = prod(mMappingMatrix,mValuesDestination[1]);
        noalias(mValuesOrigin[2]) = prod(mMappingMatrix,mValuesDestination[2]);
    }
    else
    {
        SparseSpaceType::TransposeMult(mMappingMatrix,mValuesDestination[0],mValuesOrigin[0]);
        SparseSpaceType::TransposeMult(mMappingMatrix,mValuesDestination[1],mValuesOrigin[1]);
        SparseSpaceType::TransposeMult(mMappingMatrix,mValuesDestination[2],mValuesOrigin[2]);
    }

    // Assign results to nodal variable
    for(auto& node_i : mrOriginModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);

        array_3d& r_node_vector = node_i.FastGetSolutionStepValue(rOriginVariable);
        r_node_vector(0) = mValuesOrigin[0][i];
        r_node_vector(1) = mValuesOrigin[1][i];
        r_node_vector(2) = mValuesOrigin[2][i];
    }

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphing::InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable)
{
    if (mIsMappingInitialized == false)
        Initialize();

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;

    // Prepare vectors for mapping
    mValuesOrigin[0].clear();
    mValuesDestination[0].clear();

    for(auto& node_i : mrDestinationModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);
        mValuesDestination[0][i] = node_i.FastGetSolutionStepValue(rDestinationVariable);
    }

    // Perform mapping
    if(mMapperSettings["consistent_mapping"].GetBool())
    {
        KRATOS_ERROR_IF(mrOriginModelPart.Nodes().size() != mrDestinationModelPart.Nodes().size()) << "Consistent mapping requires matching origin and destination model part.";
        noalias(mValuesOrigin[0]) = prod(mMappingMatrix,mValuesDestination[0]);
    }
    else
        SparseSpaceType::TransposeMult(mMappingMatrix,mValuesDestination[0],mValuesOrigin[0]);

    // Assign results to nodal variable
    for(auto& node_i : mrOriginModelPart.Nodes())
    {
        int i = node_i.GetValue(MAPPING_ID);
        node_i.FastGetSolutionStepValue(rOriginVariable) = mValuesOrigin[0][i];
    }

    KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphing::Update()
{
    if (mIsMappingInitialized == false)
        KRATOS_ERROR << "Mapping has to be initialized before calling the Update-function!";

    BuiltinTimer timer;
    KRATOS_INFO("ShapeOpt") << "Starting to update mapper..." << std::endl;

    CreateListOfNodesInOriginModelPart();
    InitializeMappingVariables();
    AssignMappingIds();

    ComputeMappingMatrix();

    // CalculateCurvature();

    KRATOS_INFO("ShapeOpt") << "Finished updating of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void MapperVertexMorphing::InitializeComputationOfMappingMatrix()
{
    mpSearchTree.reset();
    mMappingMatrix.clear();
}

void MapperVertexMorphing::CreateListOfNodesInOriginModelPart()
{
    mListOfNodesInOriginModelPart.resize(mrOriginModelPart.Nodes().size());
    int counter = 0;
    for (ModelPart::NodesContainerType::iterator node_it = mrOriginModelPart.NodesBegin(); node_it != mrOriginModelPart.NodesEnd(); ++node_it)
    {
        NodeTypePointer pnode = *(node_it.base());
        mListOfNodesInOriginModelPart[counter++] = pnode;
    }
}

void MapperVertexMorphing::CreateFilterFunction()
{
    std::string filter_type = mMapperSettings["filter_function_type"].GetString();

    mpFilterFunction = Kratos::make_unique<FilterFunction>(filter_type);
}

void MapperVertexMorphing::InitializeMappingVariables()
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

void MapperVertexMorphing::AssignMappingIds()
{
    unsigned int i = 0;
    for(auto& node_i : mrOriginModelPart.Nodes())
        node_i.SetValue(MAPPING_ID,i++);

    i = 0;
    for(auto& node_i : mrDestinationModelPart.Nodes())
        node_i.SetValue(MAPPING_ID,i++);
}

void MapperVertexMorphing::CreateSearchTreeWithAllNodesInOriginModelPart()
{
    mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInOriginModelPart.begin(), mListOfNodesInOriginModelPart.end(), mBucketSize));
}

void MapperVertexMorphing::ComputeMappingMatrix()
{
    InitializeComputationOfMappingMatrix();
    CreateSearchTreeWithAllNodesInOriginModelPart();

    unsigned int max_number_of_neighbors = mMapperSettings["max_nodes_in_filter_radius"].GetInt();

    for(auto& node_i : mrDestinationModelPart.Nodes())
    {
        NodeVector neighbor_nodes( max_number_of_neighbors );
        std::vector<double> resulting_squared_distances( max_number_of_neighbors );
        unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                            GetVertexMorphingRadius(node_i),
                                                                            neighbor_nodes.begin(),
                                                                            resulting_squared_distances.begin(),
                                                                            max_number_of_neighbors );



        std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
        double sum_of_weights = 0.0;

        if(number_of_neighbors >= max_number_of_neighbors)
            KRATOS_WARNING("ShapeOpt::MapperVertexMorphing") << "For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << max_number_of_neighbors << " nodes) reached!" << std::endl;

        ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
        FillMappingMatrixWithWeights( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
    }
}

void MapperVertexMorphing::ComputeWeightForAllNeighbors(  ModelPart::NodeType& origin_node,
                                    NodeVector& neighbor_nodes,
                                    unsigned int number_of_neighbors,
                                    std::vector<double>& list_of_weights,
                                    double& sum_of_weights )
{
    for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
    {
        ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
        double weight = mpFilterFunction->ComputeWeight( origin_node.Coordinates(), neighbor_node.Coordinates(), GetVertexMorphingRadius(origin_node) );

        list_of_weights[neighbor_itr] = weight;
        sum_of_weights += weight;
    }
}

void MapperVertexMorphing::FillMappingMatrixWithWeights(  ModelPart::NodeType& origin_node,
                                    NodeVector& neighbor_nodes,
                                    unsigned int number_of_neighbors,
                                    std::vector<double>& list_of_weights,
                                    double& sum_of_weights )
{


    unsigned int row_id = origin_node.GetValue(MAPPING_ID);
    for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
    {
        ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
        int collumn_id = neighbor_node.GetValue(MAPPING_ID);


        double weight = list_of_weights[neighbor_itr] / sum_of_weights;
        mMappingMatrix.insert_element(row_id,collumn_id,weight);
    }
}

// TODO: DELETE AGAIN
void MapperVertexMorphing::CalculateCurvature()
{
    // count number of quad elements
    // TODO: deal with quads 4, 8, 9 nodes
    // TODO: deal with triangles 6 nodes
    // int numQuadElms = 0;
    // for (int i = 0; i < mrDestinationModelPart.NumberOfElements; ++i) {
    //     ConditionType::Pointer condition_i = mrDestinationModelPart.mpConditions[i];
    //     if (condition_i -> NumberOfNodes() == 4) numQuadElms++;
    // }

    // Load a mesh in OFF format
    Eigen::MatrixXd V(mrDestinationModelPart.NumberOfNodes(), 3);
    // Eigen::MatrixXi F(mrDestinationModelPart->mNumElems + numQuadElms, 3);
    unsigned int i = 0;
    for (auto& node_i : mrDestinationModelPart.Nodes()) {
        V(i,0) = node_i.X();
        V(i,1) = node_i.Y();
        V(i,2) = node_i.Z();
        ++i;
    }
    KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "Size of V : " << V.size() << "..." << std::endl;
    // for (int i = 0; i < mrDestinationModelPart.NumberOfNodes(); ++i) {
    //     NodeType::Pointer node_i = mrDestinationModelPart.GetNode()[i];
    //     V(i,0) = node_i->X();
    //     V(i,1) = node_i->Y();
    //     V(i,2) = node_i->Z();
    // }

    // int counter = 0;
    Eigen::MatrixXi F(mrDestinationModelPart.NumberOfConditions(), 3);
    unsigned int j = 0;
    for (auto& condition_j : mrDestinationModelPart.Conditions()) {
        NodeType& r_node_1 = condition_j.GetGeometry()[0];
        NodeType& r_node_2 = condition_j.GetGeometry()[1];
        NodeType& r_node_3 = condition_j.GetGeometry()[2];
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_1 : " << r_node_1 << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_2 : " << r_node_2 << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_3 : " << r_node_3 << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_1_index : " << r_node_1.GetValue(MAPPING_ID) << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_2_index : " << r_node_2.GetValue(MAPPING_ID) << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_3_index : " << r_node_3.GetValue(MAPPING_ID) << "..." << std::endl;
        F(j,0) = r_node_1.GetValue(MAPPING_ID);
        F(j,1) = r_node_2.GetValue(MAPPING_ID);
        F(j,2) = r_node_3.GetValue(MAPPING_ID);
        ++j;
        // counter ++;
    }
    KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "Size of F : " << F.size() << "..." << std::endl;
    Eigen::VectorXd K;
    igl::gaussian_curvature(V,F,K);
    Eigen::SparseMatrix<double> M,Minv;
    // SparseMatrixType M,Minv;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
    igl::invert_diag(M,Minv);
    K = (Minv*K).eval();

    igl::writeOFF("test.off",V,F);

    i = 0;
    for (auto& node_i : mrDestinationModelPart.Nodes()) {
        double& gaussian_curvature_i = node_i.FastGetSolutionStepValue(GAUSSIAN_CURVATURE);
        gaussian_curvature_i = K[i];
        ++i;
    }
    // // if quad elements in mesh, add a second triangle
    // for (int j = 0; j < mpPart->mNumElems; ++j) {
    //     Element::ConstPointer element_j = mpPart->mElements[j];
    //     if (element_j -> NumberOfNodes() == 4) {
    //         F(counter,0) = mpPart->GetNodeIndex(element_j -> mNodes[0]->mId);
    //         F(counter,1) = mpPart->GetNodeIndex(element_j -> mNodes[2]->mId);
    //         F(counter,2) = mpPart->GetNodeIndex(element_j -> mNodes[3]->mId);
    //         counter ++;
    //     }
    // }
}
// END TODO: DELETE AGAIN

}  // namespace Kratos.
