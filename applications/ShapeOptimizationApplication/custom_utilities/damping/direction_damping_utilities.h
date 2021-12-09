// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef DIRECTION_DAMPING_UTILITIES_H
#define DIRECTION_DAMPING_UTILITIES_H

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
#include "custom_utilities/filter_function.h"
#include "shape_optimization_application.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Classes
///@{

/// Class for direction damping of the shape update.
/** Class for damping a specidfied direction from the shape update (and
 * according to chain rule of differentiation also from sensitvities.)
 * Can be used to constrain certain nodes of the design surface to slide on a plane.

*/

class DirectionDampingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::NodeType::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef std::vector<NodeTypePointer>::iterator NodeVectorIterator;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<double>::iterator DoubleVectorIterator;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeVectorIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of DirectionDampingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DirectionDampingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DirectionDampingUtilities( ModelPart& modelPartToDamp, Parameters DampingSettings, int MaxNeighborNodes )
        : mrModelPartToDamp( modelPartToDamp ),
          mDampingSettings( DampingSettings ),
          mMaxNeighborNodes( MaxNeighborNodes )
    {
        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;

        Parameters default_parameters( R"(
            {
                "sub_model_part_name": "MODEL_PART_NAME",
                "damping_function_type": "cosine",
                "damping_radius": -1.0,
                "direction" : [0.0, 0.0, 0.0]
            }  )" );

        KRATOS_ERROR_IF(!mDampingSettings.Has("direction")) << "DirectionDampingUtilities: 'direction' vector is missing!" << std::endl;

        // Validate against defaults -- this ensures no type mismatch
        mDampingSettings.ValidateAndAssignDefaults(default_parameters);

        KRATOS_ERROR_IF(mDampingSettings["damping_radius"].GetDouble() < 0.0) << "DirectionDampingUtilities: 'damping_radius' is a mandatory setting and has to be > 0.0!" << std::endl;

        mDirection = mDampingSettings["direction"].GetVector();
        KRATOS_ERROR_IF(norm_2(mDirection) < std::numeric_limits<double>::epsilon()) << "DirectionDampingUtilities: 'direction' vector norm is 0!" << std::endl;
        mDirection /= norm_2(mDirection);

        KRATOS_INFO("ShapeOpt") << "Creating search tree to perform direction damping..." << std::endl;
        CreateListOfNodesOfModelPart();
        CreateSearchTreeWithAllNodesOfModelPart();
        KRATOS_INFO("ShapeOpt") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;
        InitalizeDampingFactorsToHaveNoInfluence();
        SetDampingFactors();
    }

    /// Destructor.
    virtual ~DirectionDampingUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CreateListOfNodesOfModelPart()
    {
        mListOfNodesOfModelPart.resize(mrModelPartToDamp.Nodes().size());
        int counter = 0;
        for (ModelPart::NodesContainerType::iterator node_it = mrModelPartToDamp.NodesBegin(); node_it != mrModelPartToDamp.NodesEnd(); ++node_it)
        {
            NodeTypePointer p_node = *(node_it.base());
            mListOfNodesOfModelPart[counter++] = p_node;
        }
    }

    void CreateSearchTreeWithAllNodesOfModelPart()
    {
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesOfModelPart.begin(), mListOfNodesOfModelPart.end(), mBucketSize));
    }

    void InitalizeDampingFactorsToHaveNoInfluence()
    {
        mDampingFactors = std::vector<double>(mrModelPartToDamp.Nodes().size(), 1.0);
    }

    void SetDampingFactors()
    {
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting to prepare direction damping..." << std::endl;

        const std::string sub_model_part_name = mDampingSettings["sub_model_part_name"].GetString();
        const ModelPart& damping_region = mrModelPartToDamp.GetRootModelPart().GetSubModelPart(sub_model_part_name);

        const std::string damping_function_type = mDampingSettings["damping_function_type"].GetString();
        const double damping_radius = mDampingSettings["damping_radius"].GetDouble();

        const auto p_damping_function = CreateDampingFunction( damping_function_type, damping_radius );

        // Loop over all nodes in specified damping sub-model part
        for(const auto& node_i : damping_region.Nodes())
        {
            NodeVector neighbor_nodes( mMaxNeighborNodes );
            DoubleVector resulting_squared_distances( mMaxNeighborNodes,0.0 );
            const unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                                damping_radius,
                                                                                neighbor_nodes.begin(),
                                                                                resulting_squared_distances.begin(),
                                                                                mMaxNeighborNodes );

            ThrowWarningIfNodeNeighborsExceedLimit( node_i, number_of_neighbors );

            // Loop over all nodes in radius (including node on damping region itself)
            for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
            {
                ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                const double damping_factor = 1.0 - p_damping_function->ComputeWeight( node_i.Coordinates(), neighbor_node.Coordinates());

                // We check if new damping factor is smaller than the assigned one for the current node.
                // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
                const int mapping_id = neighbor_node.GetValue(MAPPING_ID);
                if(damping_factor < mDampingFactors[mapping_id])
                    mDampingFactors[mapping_id] = damping_factor;
            }
        }

        KRATOS_INFO("ShapeOpt") << "Finished preparation of direction damping." << std::endl;
    }

    FilterFunction::Pointer CreateDampingFunction( std::string damping_type, double damping_radius ) const
    {
        return Kratos::make_unique<FilterFunction>(damping_type, damping_radius);
    }

    void ThrowWarningIfNodeNeighborsExceedLimit( const ModelPart::NodeType& given_node, const unsigned int number_of_neighbors ) const
    {
        if(number_of_neighbors >= mMaxNeighborNodes)
            KRATOS_WARNING("ShapeOpt::DirectionDampingUtilities") << "For node " << given_node.Id() << " and specified damping radius, maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) reached!" << std::endl;
    }

    void DampNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        for(auto& node_i : mrModelPartToDamp.Nodes())
        {
            const int mapping_id = node_i.GetValue(MAPPING_ID);
            const double damping_factor = mDampingFactors[mapping_id];

            if (damping_factor >= 1.0)
                continue;  // skipping for efficiency

            auto& nodal_value = node_i.FastGetSolutionStepValue(rNodalVariable);

            // mDirection is already normalized
            const double projected_length = inner_prod(nodal_value, mDirection);
            const array_3d correction = - projected_length * mDirection;
            const double damping_value = 1.0 - damping_factor;

            nodal_value += correction * damping_value;
        }
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "DirectionDampingUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DirectionDampingUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    // Initialized by class constructor
    ModelPart& mrModelPartToDamp;
    Parameters mDampingSettings;
    array_3d mDirection;
    std::vector<double> mDampingFactors;

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    unsigned int mMaxNeighborNodes = 10000;
    NodeVector mListOfNodesOfModelPart;
    KDTree::Pointer mpSearchTree;

    ///@}

}; // Class DirectionDampingUtilities

///@}

}  // namespace Kratos.

#endif // DIRECTION_DAMPING_UTILITIES_H
