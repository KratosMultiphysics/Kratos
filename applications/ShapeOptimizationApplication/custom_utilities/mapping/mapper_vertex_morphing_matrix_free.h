// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_MATRIX_FREE_H
#define MAPPER_VERTEX_MORPHING_MATRIX_FREE_H

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
#include "shape_optimization_application.h"
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

class MapperVertexMorphingMatrixFree : public Mapper
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

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of MapperVertexMorphingMatrixFree
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingMatrixFree);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingMatrixFree( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : mrOriginModelPart( rOriginModelPart ),
          mrDestinationModelPart( rDestinationModelPart ),
          mMapperSettings( MapperSettings ),
          mFilterRadius( MapperSettings["filter_radius"].GetDouble() ),
          mMaxNumberOfNeighbors( MapperSettings["max_nodes_in_filter_radius"].GetInt())
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingMatrixFree()
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
        std::cout << "> Starting initialization of mapping..." << std::endl;

        if (mIsMappingInitialized == false)
        {
            CreateListOfNodesInOriginModelPart();
            CreateFilterFunction();
            InitializeMappingVariables();
            AssignMappingIds();
        }
        CreateSearchTreeWithAllNodesInOriginModelPart();

        mIsMappingInitialized = true;

        std::cout << "> Finished initialization of mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable ) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer mapping_time;
        std::cout << "\n> Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mXValuesDestination.clear();
        mYValuesDestination.clear();
        mZValuesDestination.clear();

        // Perform mapping
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            NodeVector neighbor_nodes(mMaxNumberOfNeighbors);
            std::vector<double> resulting_squared_distances(mMaxNumberOfNeighbors);
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             mFilterRadius,
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(),
                                                                             mMaxNumberOfNeighbors );

            ThrowWarningIfNumberOfNeighborsExceedsLimit(node_i, number_of_neighbors);

            std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;
            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );

            int node_i_mapping_id = node_i.GetValue(MAPPING_ID);
            for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
            {
                double weight = list_of_weights[neighbor_itr] / sum_of_weights;

                ModelPart::NodeType& node_j = *neighbor_nodes[neighbor_itr];
                array_3d& nodal_variable = node_j.FastGetSolutionStepValue(rOriginVariable);

                mXValuesDestination[node_i_mapping_id] += weight*nodal_variable[0];
                mYValuesDestination[node_i_mapping_id] += weight*nodal_variable[1];
                mZValuesDestination[node_i_mapping_id] += weight*nodal_variable[2];
            }
        }

        // Assign results to nodal variable
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);

            array_3d& r_node_vector = node_i.FastGetSolutionStepValue(rDestinationVariable);
            r_node_vector(0) = mXValuesDestination[i];
            r_node_vector(1) = mYValuesDestination[i];
            r_node_vector(2) = mZValuesDestination[i];
        }

        std::cout << "> Finished mapping in " << mapping_time.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable ) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer mapping_time;
        std::cout << "\n> Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mXValuesDestination.clear();

        // Perform mapping
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            NodeVector neighbor_nodes(mMaxNumberOfNeighbors);
            std::vector<double> resulting_squared_distances(mMaxNumberOfNeighbors);
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             mFilterRadius,
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(),
                                                                             mMaxNumberOfNeighbors );

            ThrowWarningIfNumberOfNeighborsExceedsLimit(node_i, number_of_neighbors);

            std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;
            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );

            int node_i_mapping_id = node_i.GetValue(MAPPING_ID);
            for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
            {
                double weight = list_of_weights[neighbor_itr] / sum_of_weights;
                ModelPart::NodeType& node_j = *neighbor_nodes[neighbor_itr];

                mXValuesDestination[node_i_mapping_id] += weight*node_j.FastGetSolutionStepValue(rOriginVariable);
            }
        }

        // Assign results to nodal variable
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            node_i.FastGetSolutionStepValue(rDestinationVariable) = mXValuesDestination[i];
        }

        std::cout << "> Finished mapping in " << mapping_time.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable ) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer mapping_time;
        std::cout << "\n> Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mXValuesOrigin.clear();
        mYValuesOrigin.clear();
        mZValuesOrigin.clear();

        // Perform mapping
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
            std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             mFilterRadius,
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(),
                                                                             mMaxNumberOfNeighbors );

            ThrowWarningIfNumberOfNeighborsExceedsLimit(node_i, number_of_neighbors);

            std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;
            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );

            array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rDestinationVariable);
            for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
            {
                ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
                int neighbor_node_mapping_id = neighbor_node.GetValue(MAPPING_ID);

                double weight = list_of_weights[neighbor_itr] / sum_of_weights;

                mXValuesOrigin[neighbor_node_mapping_id] += weight*nodal_variable[0];
                mYValuesOrigin[neighbor_node_mapping_id] += weight*nodal_variable[1];
                mZValuesOrigin[neighbor_node_mapping_id] += weight*nodal_variable[2];
            }
        }

        // Assign results to nodal variable
        for(auto& node_i : mrOriginModelPart.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);

            array_3d& r_node_vector = node_i.FastGetSolutionStepValue(rOriginVariable);
            r_node_vector(0) = mXValuesOrigin[i];
            r_node_vector(1) = mYValuesOrigin[i];
            r_node_vector(2) = mZValuesOrigin[i];
        }

        std::cout << "> Finished mapping in " << mapping_time.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable ) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer mapping_time;
        std::cout << "\n> Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mXValuesOrigin.clear();

        // Perform mapping
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
            std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             mFilterRadius,
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(),
                                                                             mMaxNumberOfNeighbors );

            ThrowWarningIfNumberOfNeighborsExceedsLimit(node_i, number_of_neighbors);

            std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;
            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );

            double variable_value = node_i.FastGetSolutionStepValue(rDestinationVariable);
            for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
            {
                ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
                int neighbor_node_mapping_id = neighbor_node.GetValue(MAPPING_ID);

                double weight = list_of_weights[neighbor_itr] / sum_of_weights;

                mXValuesOrigin[neighbor_node_mapping_id] += weight*variable_value;
            }
        }

        // Assign results to nodal variable
        for(auto& node_i : mrOriginModelPart.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            node_i.FastGetSolutionStepValue(rOriginVariable) = mXValuesOrigin[i];
        }

        std::cout << "> Finished mapping in " << mapping_time.ElapsedSeconds() << " s." << std::endl;
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
        return "MapperVertexMorphingMatrixFree";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingMatrixFree";
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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


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

    // Initialized by class constructor
    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    Parameters mMapperSettings;
    double mFilterRadius;
    unsigned int mMaxNumberOfNeighbors;
    FilterFunction::Pointer mpFilterFunction;

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesInOriginModelPart;
    KDTree::Pointer mpSearchTree;

    // Variables for mapping
    Vector mXValuesOrigin, mYValuesOrigin, mZValuesOrigin;
    Vector mXValuesDestination, mYValuesDestination, mZValuesDestination;
    bool mIsMappingInitialized = false;

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
        mXValuesOrigin.resize(origin_node_number,0.0);
        mYValuesOrigin.resize(origin_node_number,0.0);
        mZValuesOrigin.resize(origin_node_number,0.0);

        const unsigned int destination_node_number = mrDestinationModelPart.Nodes().size();
        mXValuesDestination.resize(destination_node_number,0.0);
        mYValuesDestination.resize(destination_node_number,0.0);
        mZValuesDestination.resize(destination_node_number,0.0);
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
        BuiltinTimer timer;
        std::cout << "> Creating search tree to perform mapping..." << std::endl;
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInOriginModelPart.begin(), mListOfNodesInOriginModelPart.end(), mBucketSize));
        std::cout << "> Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void ThrowWarningIfNumberOfNeighborsExceedsLimit(ModelPart::NodeType& given_node, unsigned int number_of_neighbors)
    {
        if(number_of_neighbors >= mMaxNumberOfNeighbors)
            std::cout << "\n> WARNING!!!!! For node " << given_node.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << mMaxNumberOfNeighbors << " nodes) reached!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void ComputeWeightForAllNeighbors(  ModelPart::NodeType& design_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double weight = mpFilterFunction->compute_weight( design_node.Coordinates(), neighbor_node.Coordinates() );

            list_of_weights[neighbor_itr] = weight;
            sum_of_weights += weight;
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
//      MapperVertexMorphingMatrixFree& operator=(MapperVertexMorphingMatrixFree const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingMatrixFree(MapperVertexMorphingMatrixFree const& rOther);


    ///@}

}; // Class MapperVertexMorphingMatrixFree

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_MATRIX_FREE_H
