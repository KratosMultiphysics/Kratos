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
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/normal_calculation_utils.h"
#include "shape_optimization_application.h"
#include "damping_function.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    typedef boost::python::extract<ModelPart&> extractModelPart;

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
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeVectorIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;    

    /// Pointer definition of DampingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DampingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DampingUtilities( ModelPart& modelPartToDamp, boost::python::dict subModelPartsForDamping, Parameters optimizationSettings )
        : mrModelPartToDamp( modelPartToDamp ),
          mrDampingRegions( subModelPartsForDamping ),
          mDampingSettings( optimizationSettings["design_variables"]["damping"]["damping_regions"] )
    {
        boost::timer timer;        
        std::cout << "\n> Creating search tree to perform damping..." << std::endl;
        CreateListOfNodesOfModelPart();
        CreateSearchTreeWithAllNodesOfModelPart();
        std::cout << "> Search tree created in: " << timer.elapsed() << " s" << std::endl;

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
        mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(mListOfNodesOfModelPart.begin(), mListOfNodesOfModelPart.end(), mBucketSize));
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
        std::cout << "\n> Starting to prepare damping..." << std::endl;

        // Loop over all regions for which damping is to be applied
        for (unsigned int regionNumber = 0; regionNumber < len(mrDampingRegions); regionNumber++)
        {
            std::string dampingRegionSubModelPartName = mDampingSettings[regionNumber]["sub_model_part_name"].GetString();
            ModelPart& dampingRegion = extractModelPart( mrDampingRegions[dampingRegionSubModelPartName] );

            bool dampX = mDampingSettings[regionNumber]["damp_X"].GetBool();
            bool dampY = mDampingSettings[regionNumber]["damp_Y"].GetBool();
            bool dampZ = mDampingSettings[regionNumber]["damp_Z"].GetBool();
            std::string dampingFunctionType = mDampingSettings[regionNumber]["damping_function_type"].GetString();
            double dampingRadius = mDampingSettings[regionNumber]["damping_radius"].GetDouble();

            DampingFunction::Pointer mpDampingFunction = CreateDampingFunction( dampingFunctionType, dampingRadius );

            // Loop over all nodes in specified damping sub-model part 
            for(auto& node_i : dampingRegion.Nodes())
            {
                ModelPart::NodeType& currentDampingNode = node_i;
                NodeVector neighbor_nodes( mMaxNeighborNodes );
                DoubleVector resulting_squared_distances( mMaxNeighborNodes,0.0 );
                unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( currentDampingNode,
                                                                                 dampingRadius, 
                                                                                 neighbor_nodes.begin(),
                                                                                 resulting_squared_distances.begin(), 
                                                                                 mMaxNeighborNodes );

                ThrowWarningIfNodeNeighborsExceedLimit( currentDampingNode, number_of_neighbors );                                                                               

                // Loop over all nodes in radius (including node on damping region itself)
                for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
                {
                    ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                    double dampingFactor = mpDampingFunction->compute_damping_factor( currentDampingNode.Coordinates(), neighbor_node.Coordinates());

                    // For every specified damping direction we check if new damping factor is smaller than the assigned one for current node. 
                    // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
                    auto& damping_factor_variable = neighbor_node.GetValue(DAMPING_FACTOR);
                    if(dampX && dampingFactor < damping_factor_variable[0])
                        damping_factor_variable[0] = dampingFactor;
                    if(dampY && dampingFactor < damping_factor_variable[1])       
                        damping_factor_variable[1] = dampingFactor;   
                    if(dampZ && dampingFactor < damping_factor_variable[2])       
                        damping_factor_variable[2] = dampingFactor;                            
                }
            }
        }

        std::cout << "> Finished preparation of damping." << std::endl;
    }

    // --------------------------------------------------------------------------
    DampingFunction::Pointer CreateDampingFunction( std::string damping_type, double damping_radius )
    {
        return boost::shared_ptr<DampingFunction>(new DampingFunction(damping_type, damping_radius));
    }    

    // --------------------------------------------------------------------------
    void ThrowWarningIfNodeNeighborsExceedLimit( ModelPart::NodeType& given_node, unsigned int number_of_neighbors )
    {
        if(number_of_neighbors >= mMaxNeighborNodes)
            std::cout << "\n> WARNING!!!!! For node " << given_node.Id() << " and specified damping radius, maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) reached!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void DampNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        for(auto& node_i : mrModelPartToDamp.Nodes())
        {   
            auto& damping_factor = node_i.GetValue(DAMPING_FACTOR);
            auto& nodalVariable = node_i.FastGetSolutionStepValue(rNodalVariable);

            nodalVariable[0] *= damping_factor[0];
            nodalVariable[1] *= damping_factor[1];
            nodalVariable[2] *= damping_factor[2];  
        }  
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
    boost::python::dict mrDampingRegions;
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
