// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner,
//                   Willer Matthias, https://github.com/matthiaswiller
//
// ==============================================================================

#ifndef DAMPING_UTILITIES_H
#define DAMPING_UTILITIES_H

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
#include "custom_utilities/filter_function.h"
#include "shape_optimization_application.h"

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

class DampingUtilities
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

    /// Pointer definition of DampingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DampingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DampingUtilities( ModelPart& modelPartToDamp, Parameters DampingSettings )
        : mrModelPartToDamp( modelPartToDamp ),
          mDampingSettings( DampingSettings ),
          mMaxNeighborNodes( DampingSettings["max_neighbor_nodes"].GetInt() )
    {
        Parameters default_parameters( R"(
            {
                "sub_model_part_name"   : "MODEL_PART_NAME",
                "damp_X"                : true,
                "damp_Y"                : true,
                "damp_Z"                : true,
                "damping_function_type" : "cosine",
                "damping_radius"        : -1.0
            }  )" );

        // Loop over all regions for which damping is to be applied
        for (auto& r_region_parameters : mDampingSettings["damping_regions"])
        {
            r_region_parameters.ValidateAndAssignDefaults(default_parameters);
            KRATOS_ERROR_IF(r_region_parameters["damping_radius"].GetDouble() < 0.0) << "DampingUtilities: 'damping_radius' is a mandatory setting and has to be > 0.0!" << std::endl;
        }

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Creating search tree to perform damping..." << std::endl;
        CreateListOfNodesOfModelPart();
        CreateSearchTreeWithAllNodesOfModelPart();
        KRATOS_INFO("ShapeOpt") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;
        InitalizeDampingFactorsToHaveNoInfluence();
        SetDampingFactorsForAllDampingRegions();
    }

    /// Destructor.
    virtual ~DampingUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void CreateListOfNodesOfModelPart()
    {
        mListOfNodesOfModelPart.resize(mrModelPartToDamp.Nodes().size());
        int counter = 0;
        for (ModelPart::NodesContainerType::iterator node_it = mrModelPartToDamp.NodesBegin(); node_it != mrModelPartToDamp.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesOfModelPart[counter++] = pnode;
        }
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesOfModelPart()
    {
        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesOfModelPart.begin(), mListOfNodesOfModelPart.end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void InitalizeDampingFactorsToHaveNoInfluence()
    {
        for(auto& node_i : mrModelPartToDamp.Nodes())
        {
            node_i.SetValue(DAMPING_FACTOR_X,1.0);
            node_i.SetValue(DAMPING_FACTOR_Y,1.0);
            node_i.SetValue(DAMPING_FACTOR_Z,1.0);
        }
    }

    // --------------------------------------------------------------------------
    void SetDampingFactorsForAllDampingRegions()
    {
        KRATOS_INFO("") << std::endl;

        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting to prepare damping..." << std::endl;

        // Loop over all regions for which damping is to be applied
        for (const auto& r_region_parameters : mDampingSettings["damping_regions"])
        {
            const std::string& sub_model_part_name = r_region_parameters["sub_model_part_name"].GetString();
            ModelPart& dampingRegion = mrModelPartToDamp.GetRootModelPart().GetSubModelPart(sub_model_part_name);

            const bool dampX = r_region_parameters["damp_X"].GetBool();
            const bool dampY = r_region_parameters["damp_Y"].GetBool();
            const bool dampZ = r_region_parameters["damp_Z"].GetBool();
            const auto& dampingFunctionType = r_region_parameters["damping_function_type"].GetString();
            const double dampingRadius = r_region_parameters["damping_radius"].GetDouble();

            const auto p_damping_function = CreateDampingFunction(dampingFunctionType, dampingRadius );

            // Loop over all nodes in specified damping sub-model part
            block_for_each(dampingRegion.Nodes(), [&](const ModelPart::NodeType& rNode) {
                NodeVector neighbor_nodes( mMaxNeighborNodes );
                const unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( rNode,
                                                                                 dampingRadius,
                                                                                 neighbor_nodes.begin(),
                                                                                 mMaxNeighborNodes );

                ThrowWarningIfNodeNeighborsExceedLimit( rNode, number_of_neighbors );

                // Loop over all nodes in radius (including node on damping region itself)
                for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
                {
                    ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                    const double dampingFactor = 1.0 - p_damping_function->ComputeWeight( rNode.Coordinates(), neighbor_node.Coordinates());

                    // For every specified damping direction we check if new damping factor is smaller than the assigned one for current node.
                    // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
                    array_3d& damping_factor_variable = neighbor_node.GetValue(DAMPING_FACTOR);
                    neighbor_node.SetLock();
                    if(dampX && dampingFactor < damping_factor_variable[0])
                        damping_factor_variable[0] = dampingFactor;
                    if(dampY && dampingFactor < damping_factor_variable[1])
                        damping_factor_variable[1] = dampingFactor;
                    if(dampZ && dampingFactor < damping_factor_variable[2])
                        damping_factor_variable[2] = dampingFactor;
                    neighbor_node.UnSetLock();
                }
            });
        }

        KRATOS_INFO("ShapeOpt") << "Finished preparation of damping in: " << timer.ElapsedSeconds() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    FilterFunction::Pointer CreateDampingFunction( std::string damping_type, double damping_radius ) const
    {
        return Kratos::make_unique<FilterFunction>(damping_type, damping_radius);
    }

    // --------------------------------------------------------------------------
    void ThrowWarningIfNodeNeighborsExceedLimit( const ModelPart::NodeType& given_node, const unsigned int number_of_neighbors ) const
    {
        if(number_of_neighbors >= mMaxNeighborNodes)
            KRATOS_WARNING("ShapeOpt::DampingUtilities") << "For node " << given_node.Id() << " and specified damping radius, maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) reached!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void DampNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        block_for_each(mrModelPartToDamp.Nodes(), [&](ModelPart::NodeType& rNode) {
            const array_3d& damping_factor = rNode.GetValue(DAMPING_FACTOR);
            array_3d& nodalVariable = rNode.FastGetSolutionStepValue(rNodalVariable);

            nodalVariable[0] *= damping_factor[0];
            nodalVariable[1] *= damping_factor[1];
            nodalVariable[2] *= damping_factor[2];
        });
    }

    // ==============================================================================

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
        return "DampingUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DampingUtilities";
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart& mrModelPartToDamp;
    Parameters mDampingSettings;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mBucketSize = 100;
    unsigned int mMaxNeighborNodes = 10000;
    NodeVector mListOfNodesOfModelPart;
    KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
//      DampingUtilities& operator=(DampingUtilities const& rOther);

    /// Copy constructor.
//      DampingUtilities(DampingUtilities const& rOther);


    ///@}

}; // Class DampingUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // DAMPIG_UTILITIES_H
