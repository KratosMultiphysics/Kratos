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

class MapperVertexMorphing : public Mapper
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
    void Initialize() override
    {
        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting initialization of mapper..." << std::endl;

        CreateListOfNodesInOriginModelPart();
        CreateFilterFunction();
        InitializeMappingVariables();
        AssignMappingIds();

        InitializeComputationOfMappingMatrix();
        CreateSearchTreeWithAllNodesInOriginModelPart();
        ComputeMappingMatrix();

        mIsMappingInitialized = true;

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

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
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

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override
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

    // --------------------------------------------------------------------------
    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
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

    // --------------------------------------------------------------------------
    void Update() override
    {
        if (mIsMappingInitialized == false)
            KRATOS_ERROR << "Mapping has to be initialized before calling the Update-function!";

        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting to update mapper..." << std::endl;

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
    KDTree::Pointer mpSearchTree;

    // Variables for mapping
    SparseMatrixType mMappingMatrix;
    std::vector<Vector> mValuesOrigin;
    std::vector<Vector> mValuesDestination;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CreateListOfNodesInOriginModelPart()
    {
        mListOfNodesInOriginModelPart.resize(mrOriginModelPart.Nodes().size());
        int counter = 0;
        for (ModelPart::NodesContainerType::iterator node_it = mrOriginModelPart.NodesBegin(); node_it != mrOriginModelPart.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesInOriginModelPart[counter++] = pnode;
        }
    }

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

        mMappingMatrix.resize(destination_node_number,origin_node_number,false);
        mMappingMatrix.clear();

        mValuesOrigin.resize(3,ZeroVector(origin_node_number));
        mValuesDestination.resize(3,ZeroVector(destination_node_number));
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
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInOriginModelPart.begin(), mListOfNodesInOriginModelPart.end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix()
    {
        double filter_radius = mMapperSettings["filter_radius"].GetDouble();
        unsigned int max_number_of_neighbors = mMapperSettings["max_nodes_in_filter_radius"].GetInt();

        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            NodeVector neighbor_nodes( max_number_of_neighbors );
            std::vector<double> resulting_squared_distances( max_number_of_neighbors );
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             filter_radius,
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
            double weight = mpFilterFunction->compute_weight( origin_node.Coordinates(), neighbor_node.Coordinates() );

            list_of_weights[neighbor_itr] = weight;
            sum_of_weights += weight;
        }
    }

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& origin_node,
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
