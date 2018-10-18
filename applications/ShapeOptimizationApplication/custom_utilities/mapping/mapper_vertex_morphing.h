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
    MapperVertexMorphing( ModelPart& rModelPart, Parameters MapperSettings )
        : mrOriginMdpa( rModelPart ),
          mrDestinationMdpa( rModelPart ),
          mOriginNodeNumber(rModelPart.Nodes().size()),
          mDestinationNodeNumber(rModelPart.Nodes().size()),
          mFilterType( MapperSettings["filter_function_type"].GetString() ),
          mFilterRadius( MapperSettings["filter_radius"].GetDouble() ),
          mMaxNumberOfNeighbors( MapperSettings["max_nodes_in_filter_radius"].GetInt() ),
          mConsistentBackwardMapping (MapperSettings["consistent_mapping_to_geometry_space"].GetBool() )
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
        std::cout << "> Starting initialization of mapping..." << std::endl;

        if (mIsMappingInitialized == false)
        {
            CreateListOfNodesInOriginMdpa();
            CreateFilterFunction();
            InitializeMappingVariables();
            AssignMappingIds();
        }

        InitializeComputationOfMappingMatrix();
        CreateSearchTreeWithAllNodesInOriginMdpa();
        ComputeMappingMatrix();

        mIsMappingInitialized = true;

        std::cout << "> Finished initialization of mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rVariable, const Variable<array_3d> &rMappedVariable ) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer mapping_time;
        std::cout << "\n> Starting mapping of " << rVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mXValuesOrigin.clear();
        mYValuesOrigin.clear();
        mZValuesOrigin.clear();
        mXValuesDestination.clear();
        mYValuesDestination.clear();
        mZValuesDestination.clear();

        for(auto& node_i : mrOriginMdpa.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rVariable);
            mXValuesOrigin[i] = nodal_variable[0];
            mYValuesOrigin[i] = nodal_variable[1];
            mZValuesOrigin[i] = nodal_variable[2];
        }

        // Perform mapping
        noalias(mXValuesDestination) = prod(mMappingMatrix,mXValuesOrigin);
        noalias(mYValuesDestination) = prod(mMappingMatrix,mYValuesOrigin);
        noalias(mZValuesDestination) = prod(mMappingMatrix,mZValuesOrigin);

        // Assign results to nodal variable
        for(auto& node_i : mrDestinationMdpa.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);

            Vector node_vector = ZeroVector(3);
            node_vector(0) = mXValuesDestination[i];
            node_vector(1) = mYValuesDestination[i];
            node_vector(2) = mZValuesDestination[i];
            node_i.FastGetSolutionStepValue(rMappedVariable) = node_vector;
        }

        std::cout << "> Finished mapping in " << mapping_time.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rVariable, const Variable<array_3d> &rMappedVariable ) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer mapping_time;
        std::cout << "\n> Starting inverse mapping of " << rVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mXValuesOrigin.clear();
        mYValuesOrigin.clear();
        mZValuesOrigin.clear();
        mXValuesDestination.clear();
        mYValuesDestination.clear();
        mZValuesDestination.clear();

        for(auto& node_i : mrDestinationMdpa.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rVariable);
            mXValuesDestination[i] = nodal_variable[0];
            mYValuesDestination[i] = nodal_variable[1];
            mZValuesDestination[i] = nodal_variable[2];
        }

        // Perform mapping
        if (mConsistentBackwardMapping)
        {
            noalias(mXValuesOrigin) = prod(mMappingMatrix,mXValuesDestination);
            noalias(mYValuesOrigin) = prod(mMappingMatrix,mYValuesDestination);
            noalias(mZValuesOrigin) = prod(mMappingMatrix,mZValuesDestination);
        }
        else
        {
            SparseSpaceType::TransposeMult(mMappingMatrix,mXValuesDestination,mXValuesOrigin);
            SparseSpaceType::TransposeMult(mMappingMatrix,mYValuesDestination,mYValuesOrigin);
            SparseSpaceType::TransposeMult(mMappingMatrix,mZValuesDestination,mZValuesOrigin);
        }

        // Assign results to nodal variable
        for(auto& node_i : mrOriginMdpa.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);

            Vector node_vector = ZeroVector(3);
            node_vector(0) = mXValuesOrigin[i];
            node_vector(1) = mYValuesOrigin[i];
            node_vector(2) = mZValuesOrigin[i];
            node_i.FastGetSolutionStepValue(rMappedVariable) = node_vector;
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

    ModelPart& mrOriginMdpa;
    ModelPart& mrDestinationMdpa;
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

    // Initialized by class constructor
    const unsigned int mOriginNodeNumber;
    const unsigned int mDestinationNodeNumber;
    std::string mFilterType;
    double mFilterRadius;
    unsigned int mMaxNumberOfNeighbors;
    bool mConsistentBackwardMapping;

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesInOriginMdpa;
    KDTree::Pointer mpSearchTree;

    // Variables for mapping
    SparseMatrixType mMappingMatrix;
    Vector mXValuesOrigin, mYValuesOrigin, mZValuesOrigin;
    Vector mXValuesDestination, mYValuesDestination, mZValuesDestination;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CreateListOfNodesInOriginMdpa()
    {
        mListOfNodesInOriginMdpa.resize(mOriginNodeNumber);
        int counter = 0;
        for (ModelPart::NodesContainerType::iterator node_it = mrOriginMdpa.NodesBegin(); node_it != mrOriginMdpa.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesInOriginMdpa[counter++] = pnode;
        }
    }

    // --------------------------------------------------------------------------
    void CreateFilterFunction()
    {
        mpFilterFunction = Kratos::shared_ptr<FilterFunction>(new FilterFunction(mFilterType, mFilterRadius));
    }

    // --------------------------------------------------------------------------
    void InitializeMappingVariables()
    {
        mMappingMatrix.resize(mDestinationNodeNumber,mOriginNodeNumber,false);
        mMappingMatrix.clear();

        mXValuesOrigin.resize(mOriginNodeNumber,0.0);
        mYValuesOrigin.resize(mOriginNodeNumber,0.0);
        mZValuesOrigin.resize(mOriginNodeNumber,0.0);

        mXValuesDestination.resize(mDestinationNodeNumber,0.0);
        mYValuesDestination.resize(mDestinationNodeNumber,0.0);
        mZValuesDestination.resize(mDestinationNodeNumber,0.0);
    }

    // --------------------------------------------------------------------------
    void AssignMappingIds()
    {
        unsigned int i = 0;
        for(auto& node_i : mrOriginMdpa.Nodes())
            node_i.SetValue(MAPPING_ID,i++);
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesInOriginMdpa()
    {
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInOriginMdpa.begin(), mListOfNodesInOriginMdpa.end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix()
    {
        for(auto& node_i : mrDestinationMdpa.Nodes())
        {
            NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
            std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             mFilterRadius,
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(),
                                                                             mMaxNumberOfNeighbors );

            std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;

            ThrowWarningIfMaxNodeNeighborsReached( node_i, number_of_neighbors );
            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            FillMappingMatrixWithWeights( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
        }
    }

    // --------------------------------------------------------------------------
    void ThrowWarningIfMaxNodeNeighborsReached( ModelPart::NodeType& given_node, unsigned int number_of_neighbors )
    {
        if(number_of_neighbors >= mMaxNumberOfNeighbors)
            std::cout << "\n> WARNING!!!!! For node " << given_node.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << mMaxNumberOfNeighbors << " nodes) reached!" << std::endl;
    }

    // --------------------------------------------------------------------------
    virtual void ComputeWeightForAllNeighbors(  ModelPart::NodeType& node_of_interest,
                                                NodeVector& neighbor_nodes,
                                                unsigned int number_of_neighbors,
                                                std::vector<double>& list_of_weights,
                                                double& sum_of_weights )
    {
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double weight = mpFilterFunction->compute_weight( node_of_interest.Coordinates(), neighbor_node.Coordinates() );

            list_of_weights[neighbor_itr] = weight;
            sum_of_weights += weight;
        }
    }

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& node_of_interest,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        unsigned int row_id = node_of_interest.GetValue(MAPPING_ID);
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
