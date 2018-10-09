// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_GENERALIZED_VERTEX_MORPHING_H
#define MAPPER_GENERALIZED_VERTEX_MORPHING_H

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

class MapperGeneralizedVertexMorphing
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

    /// Pointer definition of MapperGeneralizedVertexMorphing
    KRATOS_CLASS_POINTER_DEFINITION(MapperGeneralizedVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperGeneralizedVertexMorphing( ModelPart& origin_mdpa, ModelPart& destination_mdpa, Parameters mapper_settings )
        : mrOriginMdpa(origin_mdpa),
          mrDestinationMdpa(destination_mdpa),
          mOriginSize(origin_mdpa.Nodes().size()),
          mDestinationSize(destination_mdpa.Nodes().size()),
          mFilterType( mapper_settings["filter_function_type"].GetString() ),
          mFilterRadius( mapper_settings["filter_radius"].GetDouble() ),
          mMaxNumberOfNeighbors( mapper_settings["max_nodes_in_filter_radius"].GetInt() )
    {
        CreateListOfNodesInOriginMdpa();
        CreateFilterFunction();
        InitializeMappingVariables();
        AssignMappingIds();
        InitializeComputationOfMappingMatrix();
        ComputeMappingMatrix();
    }

    /// Destructor.
    virtual ~MapperGeneralizedVertexMorphing()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rVariable, const Variable<array_3d> &rMappedVariable)
    {
        BuiltinTimer mapping_time;
        std::cout << "\n> Starting mapping..." << std::endl;

        PrepareVectorsForMapping( rVariable );
        PerformMultiplicationWithMappingMatrix();
        AssignMappingResultsToNodalVariable( rMappedVariable, false );

        std::cout << "> Time needed for mapping: " << mapping_time.ElapsedSeconds() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rVariable, const Variable<array_3d> &rMappedVariable)
    {
        BuiltinTimer mapping_time;
        std::cout << "\n> Starting inverse mapping..." << std::endl;

        PrepareVectorsForInverseMapping( rVariable );
        PerformMultiplicationWithTransposeMappingMatrix();
        AssignMappingResultsToNodalVariable( rMappedVariable, true );

        std::cout << "> Time needed for inverse mapping: " << mapping_time.ElapsedSeconds() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void UpdateMappingMatrix()
    {
        InitializeComputationOfMappingMatrix();
        ComputeMappingMatrix();
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
    virtual std::string Info() const
    {
        return "MapperGeneralizedVertexMorphing";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperGeneralizedVertexMorphing";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    ModelPart& mrOriginMdpa;
    ModelPart& mrDestinationMdpa;
    FilterFunction::Pointer mpFilterFunction;

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
    const unsigned int mOriginSize;
    const unsigned int mDestinationSize;
    std::string mFilterType;
    double mFilterRadius;
    unsigned int mMaxNumberOfNeighbors;

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesInOriginMdpa;
    KDTree::Pointer mpSearchTree;

    // Variables for mapping
    SparseMatrixType mMappingMatrix;
    Vector x_values_in_original_mdpa, y_values_in_original_mdpa, z_values_in_original_mdpa;
    Vector x_values_in_destination_mdpa, y_values_in_destination_mdpa, z_values_in_destination_mdpa;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CreateListOfNodesInOriginMdpa()
    {
        mListOfNodesInOriginMdpa.resize(mOriginSize);
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
        mMappingMatrix.resize(mDestinationSize,mOriginSize,false);
        mMappingMatrix.clear();

        x_values_in_original_mdpa.resize(mOriginSize,0.0);
        y_values_in_original_mdpa.resize(mOriginSize,0.0);
        z_values_in_original_mdpa.resize(mOriginSize,0.0);

        x_values_in_destination_mdpa.resize(mDestinationSize,0.0);
        y_values_in_destination_mdpa.resize(mDestinationSize,0.0);
        z_values_in_destination_mdpa.resize(mDestinationSize,0.0);
    }

    // --------------------------------------------------------------------------
    void AssignMappingIds()
    {
        unsigned int i = 0;
        for(auto& node_i : mrOriginMdpa.Nodes())
            node_i.SetValue(MAPPING_ID,i++);

        i = 0;
        for(auto& node_i : mrDestinationMdpa.Nodes())
            node_i.SetValue(MAPPING_ID,i++);
    }

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix()
    {
        BuiltinTimer timer;
        std::cout << "> Computing mapping matrix to perform mapping..." << std::endl;

        CreateSearchTreeWithAllNodesInOriginMdpa();
        ComputeEntriesOfMappingMatrix();

        std::cout << "> Mapping matrix computed in: " << timer.ElapsedSeconds() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesInOriginMdpa()
    {
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInOriginMdpa.begin(), mListOfNodesInOriginMdpa.end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void ComputeEntriesOfMappingMatrix()
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
    void PrepareVectorsForInverseMapping( const Variable<array_3d> &rNodalVariable )
    {
        x_values_in_original_mdpa.clear();
        y_values_in_original_mdpa.clear();
        z_values_in_original_mdpa.clear();
        x_values_in_destination_mdpa.clear();
        y_values_in_destination_mdpa.clear();
        z_values_in_destination_mdpa.clear();

        for(auto& node_i : mrDestinationMdpa.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
            x_values_in_destination_mdpa[i] = nodal_variable[0];
            y_values_in_destination_mdpa[i] = nodal_variable[1];
            z_values_in_destination_mdpa[i] = nodal_variable[2];
        }
    }

    // --------------------------------------------------------------------------
    void PrepareVectorsForMapping( const Variable<array_3d> &rNodalVariable )
    {
        x_values_in_original_mdpa.clear();
        y_values_in_original_mdpa.clear();
        z_values_in_original_mdpa.clear();
        x_values_in_destination_mdpa.clear();
        y_values_in_destination_mdpa.clear();
        z_values_in_destination_mdpa.clear();

        for(auto& node_i : mrOriginMdpa.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
            x_values_in_original_mdpa[i] = nodal_variable[0];
            y_values_in_original_mdpa[i] = nodal_variable[1];
            z_values_in_original_mdpa[i] = nodal_variable[2];
        }
    }

    // --------------------------------------------------------------------------
    void PerformMultiplicationWithMappingMatrix()
    {
        noalias(x_values_in_destination_mdpa) = prod(mMappingMatrix,x_values_in_original_mdpa);
        noalias(y_values_in_destination_mdpa) = prod(mMappingMatrix,y_values_in_original_mdpa);
        noalias(z_values_in_destination_mdpa) = prod(mMappingMatrix,z_values_in_original_mdpa);
    }

    // --------------------------------------------------------------------------
    void PerformMultiplicationWithTransposeMappingMatrix()
    {
        SparseSpaceType::TransposeMult(mMappingMatrix,x_values_in_destination_mdpa,x_values_in_original_mdpa);
        SparseSpaceType::TransposeMult(mMappingMatrix,y_values_in_destination_mdpa,y_values_in_original_mdpa);
        SparseSpaceType::TransposeMult(mMappingMatrix,z_values_in_destination_mdpa,z_values_in_original_mdpa);
    }

    // --------------------------------------------------------------------------
    void AssignMappingResultsToNodalVariable( const Variable<array_3d> &rNodalVariable, bool is_inverse_map )
    {
        if(is_inverse_map)
        {
            for(auto& node_i : mrOriginMdpa.Nodes())
            {
                int i = node_i.GetValue(MAPPING_ID);

                Vector node_vector = ZeroVector(3);
                node_vector(0) = x_values_in_original_mdpa[i];
                node_vector(1) = y_values_in_original_mdpa[i];
                node_vector(2) = z_values_in_original_mdpa[i];
                node_i.FastGetSolutionStepValue(rNodalVariable) = node_vector;
            }
        }
        else
        {
            for(auto& node_i : mrDestinationMdpa.Nodes())
            {
                int i = node_i.GetValue(MAPPING_ID);

                Vector node_vector = ZeroVector(3);
                node_vector(0) = x_values_in_destination_mdpa[i];
                node_vector(1) = y_values_in_destination_mdpa[i];
                node_vector(2) = z_values_in_destination_mdpa[i];
                node_i.FastGetSolutionStepValue(rNodalVariable) = node_vector;
            }
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
//      MapperGeneralizedVertexMorphing& operator=(MapperGeneralizedVertexMorphing const& rOther);

    /// Copy constructor.
//      MapperGeneralizedVertexMorphing(MapperGeneralizedVertexMorphing const& rOther);


    ///@}

}; // Class MapperGeneralizedVertexMorphing

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_GENERALIZED_VERTEX_MORPHING_H
