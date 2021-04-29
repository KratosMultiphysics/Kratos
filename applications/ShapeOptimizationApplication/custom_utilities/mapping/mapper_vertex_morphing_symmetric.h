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
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "filter_function.h"
#include "symmetry_base.h"
#include "symmetry_plane.h"
#include "symmetry_revolution.h"

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

class MapperVertexMorphingSymmetric : public Mapper
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

    /// Pointer definition of MapperVertexMorphingSymmetric
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingSymmetric);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingSymmetric( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : mrOriginModelPart(rOriginModelPart),
          mrDestinationModelPart(rDestinationModelPart),
          mMapperSettings(MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingSymmetric()
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
        if (mIsMappingInitialized == false) {
            KRATOS_ERROR << "Mapping has to be initialized before calling the Update-function!";
        }

        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting to update mapper..." << std::endl;

        InitializeMappingVariables();
        AssignMappingIds();

        if (mMapperSettings["plane_symmetry"].GetBool()) {
            mpSymmetry = Kratos::make_shared<SymmetryPlane>(mrOriginModelPart, mrDestinationModelPart, mMapperSettings["plane_symmetry_settings"]);
        } else if (mMapperSettings["revolution"].GetBool()) {
            mpSymmetry = Kratos::make_shared<SymmetryRevolution>(mrOriginModelPart, mrDestinationModelPart, mMapperSettings["revolution_settings"]);
        } else {
            KRATOS_ERROR << "No symmetry type specified" << std::endl;
        }

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

    SymmetryBase::Pointer mpSymmetry;

    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CreateFilterFunction()
    {
        const std::string filter_type = mMapperSettings["filter_function_type"].GetString();
        const double filter_radius = mMapperSettings["filter_radius"].GetDouble();

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
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mpSymmetry->GetOriginSearchNodes().begin(), mpSymmetry->GetOriginSearchNodes().end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix()
    {
        const double filter_radius = mMapperSettings["filter_radius"].GetDouble();
        const unsigned int max_number_of_neighbors = mMapperSettings["max_nodes_in_filter_radius"].GetInt();

        // variable size, fixed capacity vecs
        std::vector<bool> transform;
        NodeVector total_neighbor_nodes;
        std::vector<double> total_list_of_weights;
        std::vector<double> list_of_weights;
        transform.reserve(max_number_of_neighbors);
        total_neighbor_nodes.reserve( max_number_of_neighbors );
        total_list_of_weights.reserve( max_number_of_neighbors );
        list_of_weights.reserve( max_number_of_neighbors );

        // fixed size vecs
        std::vector<double> resulting_squared_distances( max_number_of_neighbors );
        NodeVector neighbor_nodes( max_number_of_neighbors );

        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            transform.clear();
            total_neighbor_nodes.clear();
            total_list_of_weights.clear();

            auto search_nodes = mpSymmetry->GetDestinationSearchNodes(node_i.GetValue(MAPPING_ID));
            unsigned int total_number_of_neighbors = 0;
            double total_sum_of_weights = 0;
            for (auto& pair_i : search_nodes) {
                NodeType search_node(0, pair_i.first);
                const unsigned int number_of_neighbors = mpSearchTree->SearchInRadius(search_node,
                                                                                filter_radius,
                                                                                neighbor_nodes.begin(),
                                                                                resulting_squared_distances.begin(),
                                                                                max_number_of_neighbors );
                total_number_of_neighbors += number_of_neighbors;
                transform.resize(total_number_of_neighbors, pair_i.second);

                list_of_weights.clear();
                list_of_weights.resize( number_of_neighbors, 0.0 );
                ComputeWeightForAllNeighbors( search_node, neighbor_nodes, number_of_neighbors, list_of_weights, total_sum_of_weights );

                total_list_of_weights.insert(total_list_of_weights.end(), list_of_weights.begin(), list_of_weights.begin()+number_of_neighbors);
                total_neighbor_nodes.insert(total_neighbor_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.begin()+number_of_neighbors);
            }

            if(total_number_of_neighbors >= max_number_of_neighbors)
                KRATOS_WARNING("ShapeOpt::MapperVertexMorphingSymmetric") << "For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << max_number_of_neighbors << " nodes) reached!" << std::endl;

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
            const NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            const double weight = mpFilterFunction->compute_weight( destination_node.Coordinates(), neighbor_node.Coordinates() );

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
        const unsigned int row_id = destination_node.GetValue(MAPPING_ID);

        std::vector<unsigned int> unique_col_ids(number_of_neighbors);
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            const NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            unique_col_ids[neighbor_itr] = neighbor_node.GetValue(MAPPING_ID);
        }
        std::sort (unique_col_ids.begin(), unique_col_ids.end());
        unique_col_ids.erase( unique( unique_col_ids.begin(), unique_col_ids.end() ), unique_col_ids.end() );

        // initialize non-zeros row by row - sorted access helps a lot
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
        }
    }

    // --------------------------------------------------------------------------

}; // Class MapperVertexMorphingSymmetric

}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_SYMMETRIC_H
