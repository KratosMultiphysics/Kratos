// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_REVOLUTION_H
#define MAPPER_VERTEX_MORPHING_REVOLUTION_H

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
#include "utilities/math_utils.h"
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "filter_function.h"

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


class RevolutionModelPart
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(RevolutionModelPart);

    typedef Node<3> NodeType;
    typedef NodeType::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef array_1d<double,3> array_3d;


    RevolutionModelPart(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, array_3d Point,  array_3d Axis)
    : mrOriginModelPart(rOriginModelPart), mrDestinationModelPart(rDestinationModelPart), mPoint(Point), mAxis(Axis)
    {
        mPlaneVector1 = ZeroVector(3);
        int index = 0;
        double max_value = 0;
        for (int i=0; i<3; ++i){
            if (std::abs(mAxis[i]) > max_value){
                index = i;
                max_value = std::abs(mAxis[i]);
            }
        }

        const int index2 = (index == 2) ? 0 : index + 1;

        mPlaneVector1[index2] = mAxis[index];
        mPlaneVector1[index] = -mAxis[index2];

        // create transformed node vecs
        mOriginNodes.resize(mrOriginModelPart.Nodes().size());
        mTransformedOriginNodes.resize(mrOriginModelPart.Nodes().size());
        for (auto& r_node_i : mrOriginModelPart.Nodes()){
            const int mapping_id = r_node_i.GetValue(MAPPING_ID);
            mOriginNodes[mapping_id] = &r_node_i;
            mTransformedOriginNodes[mapping_id] = GetTransformedNode(r_node_i);
        }

        mDestinationNodes.resize(mrDestinationModelPart.Nodes().size());
        mTransformedDestinationNodes.resize(mrDestinationModelPart.Nodes().size());
        for (auto& r_node_i : mrDestinationModelPart.Nodes()){
            const int mapping_id = r_node_i.GetValue(MAPPING_ID);
            mDestinationNodes[mapping_id] = &r_node_i;
            mTransformedDestinationNodes[mapping_id] = GetTransformedNode(r_node_i);
        }
    }

    NodeVector& GetOriginSearchNodes(){
        return mTransformedOriginNodes;
    }

    std::vector<std::pair<array_3d, bool>> GetDestinationSearchNodes(const size_t MappingId){
        return {
            std::make_pair(mTransformedDestinationNodes[MappingId]->Coordinates(), true),
        };
    }

    BoundedMatrix<double, 3, 3> TransformationMatrix(const size_t DestinationMappingId, const size_t OriginMappingId) const
    {
        const NodeType& r_origin_node = *mOriginNodes[OriginMappingId];
        const NodeType& r_destination_node = *mDestinationNodes[DestinationMappingId];

        // find angle between ortho
        const array_3d origin_v1 = r_origin_node.Coordinates() - mPoint;
        const array_3d origin_along_axis = inner_prod(mAxis, origin_v1) * mAxis;
        array_3d origin_ortho = origin_v1 - origin_along_axis;

        const double norm_origin = norm_2(origin_ortho);
        if (norm_origin < std::numeric_limits<double>::epsilon()) {
            BoundedMatrix<double, 3, 3 > r = ZeroMatrix(3,3);
            r(0,0) = mAxis[0];
            r(1,1) = mAxis[1];
            r(2,2) = mAxis[2];
            return r;
        }
        origin_ortho /= norm_origin;

        const array_3d destination_v1 = r_destination_node.Coordinates() - mPoint;
        const array_3d destination_along_axis = inner_prod(mAxis, destination_v1) * mAxis;
        array_3d destination_ortho = destination_v1 - destination_along_axis;

        const double norm_destination = norm_2(destination_ortho);
        if (norm_destination < std::numeric_limits<double>::epsilon()) {
            BoundedMatrix<double, 3, 3 > r = ZeroMatrix(3,3);
            r(0,0) = mAxis[0];
            r(1,1) = mAxis[1];
            r(2,2) = mAxis[2];
            return r;
        }
        destination_ortho /= norm_destination;

        double dot_prod = inner_prod(origin_ortho, destination_ortho);
        // result of inner_prod can be slightly larger (smaller) than 1 (-1). This results in NAN for acos!
        if (dot_prod >= 1.0) dot_prod = 1.0;
        else if (dot_prod <= -1.0) dot_prod = -1.0;

        double angle = acos(dot_prod);
        if (inner_prod(mAxis, MathUtils<double>::CrossProduct(origin_ortho, destination_ortho)) < 0) { // Or > 0
            angle = -angle;
        }

        BoundedMatrix<double, 3, 3 > r;
		const double c=cos(angle);
		const double s=sin(angle);
        const double t = 1-c;

        r(0,0) = t*mAxis[0]*mAxis[0] + c;          r(0,1) = t*mAxis[0]*mAxis[1] - mAxis[2]*s; r(0,2) = t*mAxis[0]*mAxis[2] + mAxis[1]*s;
        r(1,0) = t*mAxis[0]*mAxis[1] + mAxis[2]*s; r(1,1) = t*mAxis[1]*mAxis[1] + c;          r(1,2) = t*mAxis[1]*mAxis[2] - mAxis[0]*s;
        r(2,0) = t*mAxis[0]*mAxis[2] - mAxis[1]*s; r(2,1) = t*mAxis[1]*mAxis[2] + mAxis[0]*s; r(2,2) = t*mAxis[2]*mAxis[2] + c;

        return r;
    }

    NodeTypePointer GetTransformedNode(NodeType& rNode) {
        NodeTypePointer p_new_node = Kratos::make_intrusive<NodeType>(rNode.Id(), rNode[0], rNode[1], rNode[2]);
        p_new_node->SetValue(MAPPING_ID, rNode.GetValue(MAPPING_ID));

        const array_3d v1 = p_new_node->Coordinates() - mPoint;
        const array_3d along_axis = inner_prod(mAxis, v1) * mAxis;
        const array_3d ortho = v1 - along_axis;

        p_new_node->Coordinates() = mPoint + along_axis + mPlaneVector1 * norm_2(ortho);

        return p_new_node;
    }

    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    array_3d mPoint;
    array_3d mAxis;
    array_3d mPlaneVector1; // vector in the plane orthogonal to the mAxis

    NodeVector mOriginNodes;
    NodeVector mDestinationNodes;
    NodeVector mTransformedOriginNodes;
    NodeVector mTransformedDestinationNodes;


}; // Class RevolutionModelPart


class MapperVertexMorphingRevolution : public Mapper
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

    /// Pointer definition of MapperVertexMorphingRevolution
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingRevolution);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingRevolution( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : mrOriginModelPart(rOriginModelPart),
          mrDestinationModelPart(rDestinationModelPart),
          mMapperSettings(MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingRevolution()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override
    {
        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting initialization of mapper..." << std::endl;

        CreateFilterFunction();

        mIsMappingInitialized = true;

        Update();

        KRATOS_INFO("ShapeOpt") << "Finished initialization of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
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

        for(auto& node_i : mrOriginModelPart.Nodes())
        {
            const int i = node_i.GetValue(MAPPING_ID);
            const array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rOriginVariable);
            values_origin[i*3+0] = r_nodal_variable[0];
            values_origin[i*3+1] = r_nodal_variable[1];
            values_origin[i*3+2] = r_nodal_variable[2];
        }

        // Perform mapping
        noalias(values_destination) = prod(mMappingMatrix,values_origin);

        // Assign results to nodal variable
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            const int i = node_i.GetValue(MAPPING_ID);

            array_3d& r_node_vector = node_i.FastGetSolutionStepValue(rDestinationVariable);
            r_node_vector(0) = values_destination[i*3+0];
            r_node_vector(1) = values_destination[i*3+1];
            r_node_vector(2) = values_destination[i*3+2];
        }

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
    {
        KRATOS_ERROR << "Scalar mapping not possible." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;

        Vector values_origin(mrOriginModelPart.Nodes().size()*3);
        Vector values_destination(mrDestinationModelPart.Nodes().size()*3);

        // Prepare vectors for mapping
        values_origin.clear();
        values_destination.clear();

        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            const int i = node_i.GetValue(MAPPING_ID);
            const array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rDestinationVariable);
            values_destination[i*3+0] = r_nodal_variable[0];
            values_destination[i*3+1] = r_nodal_variable[1];
            values_destination[i*3+2] = r_nodal_variable[2];
        }

        SparseSpaceType::TransposeMult(mMappingMatrix,values_destination,values_origin);

        // Assign results to nodal variable
        for(auto& node_i : mrOriginModelPart.Nodes())
        {
            const int i = node_i.GetValue(MAPPING_ID);

            array_3d& r_node_vector = node_i.FastGetSolutionStepValue(rOriginVariable);
            r_node_vector(0) = values_origin[i*3+0];
            r_node_vector(1) = values_origin[i*3+1];
            r_node_vector(2) = values_origin[i*3+2];
        }

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        KRATOS_ERROR << "Scalar mapping not possible." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Update() override
    {
        if (mIsMappingInitialized == false)
            KRATOS_ERROR << "Mapping has to be initialized before calling the Update-function!";

        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting to update mapper..." << std::endl;

        InitializeMappingVariables();
        AssignMappingIds();

        const array_3d point = mMapperSettings["revolution_settings"]["point"].GetVector();
        array_3d axis = mMapperSettings["revolution_settings"]["axis"].GetVector();
        axis /= norm_2(axis);
        mpRevolution = Kratos::make_shared<RevolutionModelPart>(mrOriginModelPart, mrDestinationModelPart, point, axis);

        InitializeComputationOfMappingMatrix();
        CreateSearchTreeWithAllNodesInOriginModelPart();
        ComputeMappingMatrix();

        KRATOS_INFO("ShapeOpt") << "Finished updating of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

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
        return "MapperVertexMorphingRevolution";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingRevolution";
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
    FilterFunction::Pointer mpFilterFunction;
    bool mIsMappingInitialized = false;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual void InitializeComputationOfMappingMatrix()
    {
        mpSearchTree.reset();
        mMappingMatrix.clear();
    }

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
    KDTree::Pointer mpSearchTree;

    // Variables for mapping
    SparseMatrixType mMappingMatrix;

    RevolutionModelPart::Pointer mpRevolution;

    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CreateFilterFunction()
    {
        std::string filter_type = mMapperSettings["filter_function_type"].GetString();
        double filter_radius = mMapperSettings["filter_radius"].GetDouble();

        mpFilterFunction = Kratos::shared_ptr<FilterFunction>(new FilterFunction(filter_type, filter_radius));
    }

    // --------------------------------------------------------------------------
    void InitializeMappingVariables()
    {
        const unsigned int origin_node_number = mrOriginModelPart.Nodes().size();
        const unsigned int destination_node_number = mrDestinationModelPart.Nodes().size();
        mMappingMatrix.resize(destination_node_number*3,origin_node_number*3,false);
    }

    // --------------------------------------------------------------------------
    void AssignMappingIds()
    {
        unsigned int i = 0;
        for(auto& node_i : mrOriginModelPart.Nodes())
            node_i.SetValue(MAPPING_ID,i++);

        i = 0;
        for(auto& node_i : mrDestinationModelPart.Nodes())
            node_i.SetValue(MAPPING_ID,i++);
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesInOriginModelPart()
    {
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mpRevolution->GetOriginSearchNodes().begin(), mpRevolution->GetOriginSearchNodes().end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix()
    {
        double filter_radius = mMapperSettings["filter_radius"].GetDouble();
        unsigned int max_number_of_neighbors = mMapperSettings["max_nodes_in_filter_radius"].GetInt();

        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            auto search_nodes = mpRevolution->GetDestinationSearchNodes(node_i.GetValue(MAPPING_ID));
            std::vector<bool> transform;
            unsigned int total_number_of_neighbors = 0;
            NodeVector total_neighbor_nodes;
            std::vector<double> total_list_of_weights;
            double total_sum_of_weights = 0;
            double sum_of_weights = 0.0;
            for (auto& pair_i : search_nodes) {
                NodeType search_node(0, pair_i.first);
                NodeVector neighbor_nodes( max_number_of_neighbors );
                std::vector<double> resulting_squared_distances( max_number_of_neighbors );
                unsigned int number_of_neighbors = mpSearchTree->SearchInRadius(search_node,
                                                                                filter_radius,
                                                                                neighbor_nodes.begin(),
                                                                                resulting_squared_distances.begin(),
                                                                                max_number_of_neighbors );
                total_number_of_neighbors += number_of_neighbors;
                transform.resize(total_number_of_neighbors, pair_i.second);

                std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
                ComputeWeightForAllNeighbors( search_node, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );

                total_list_of_weights.insert(total_list_of_weights.end(), list_of_weights.begin(), list_of_weights.begin()+number_of_neighbors);
                total_neighbor_nodes.insert(total_neighbor_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.begin()+number_of_neighbors);
                total_sum_of_weights += sum_of_weights;
            }

            if(total_number_of_neighbors >= max_number_of_neighbors)
                KRATOS_WARNING("ShapeOpt::MapperVertexMorphingRevolution") << "For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << max_number_of_neighbors << " nodes) reached!" << std::endl;

            FillMappingMatrixWithWeights( node_i, total_neighbor_nodes, total_number_of_neighbors, total_list_of_weights, transform, total_sum_of_weights );
        }
    }

    // --------------------------------------------------------------------------
    virtual void ComputeWeightForAllNeighbors(  ModelPart::NodeType& destination_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double weight = mpFilterFunction->compute_weight( destination_node.Coordinates(), neighbor_node.Coordinates() );

            list_of_weights[neighbor_itr] = weight;
            sum_of_weights += weight;
        }
    }

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& destination_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        std::vector<bool>& transform,
                                        double& sum_of_weights )
    {
        unsigned int row_id = destination_node.GetValue(MAPPING_ID);

        std::vector<double> col_ids(number_of_neighbors);
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            const int collumn_id = neighbor_node.GetValue(MAPPING_ID);
            col_ids[neighbor_itr] = collumn_id;
        }
        std::sort (col_ids.begin(), col_ids.end());

        // initialize non-zeros row by row - sorted access helps a lot
        for (int i: col_ids) {
            mMappingMatrix.insert_element(row_id*3+0, i*3+0, 0.0);
            mMappingMatrix.insert_element(row_id*3+0, i*3+1, 0.0);
            mMappingMatrix.insert_element(row_id*3+0, i*3+2, 0.0);
        }
        for (int i: col_ids) {
            mMappingMatrix.insert_element(row_id*3+1, i*3+1, 0.0);
            mMappingMatrix.insert_element(row_id*3+1, i*3+1, 0.0);
            mMappingMatrix.insert_element(row_id*3+1, i*3+2, 0.0);
        }
        for (int i: col_ids) {
            mMappingMatrix.insert_element(row_id*3+2, i*3+2, 0.0);
            mMappingMatrix.insert_element(row_id*3+2, i*3+1, 0.0);
            mMappingMatrix.insert_element(row_id*3+3, i*3+2, 0.0);
        }

        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            const int collumn_id = neighbor_node.GetValue(MAPPING_ID);

            auto mat = (transform[neighbor_itr]) ? mpRevolution->TransformationMatrix(row_id, collumn_id) : IdentityMatrix(3);
            const double weight = list_of_weights[neighbor_itr] / sum_of_weights;
            for (int i=0; i<3; ++i){
                for (int j=0; j<3; ++j) {
                    mMappingMatrix(row_id*3+i, collumn_id*3+j) += weight*mat(i,j);
                }
            }
        }
    }

    // --------------------------------------------------------------------------

}; // Class MapperVertexMorphingRevolution

}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_SYMMETRIC_H
