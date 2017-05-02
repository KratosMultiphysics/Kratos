// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_ITERATIVE_H
#define MAPPER_VERTEX_MORPHING_ITERATIVE_H

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
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "../../kratos/spatial_containers/spatial_containers.h"
#include "../../kratos/utilities/timer.h"
#include "../../kratos/processes/node_erase_process.h"
#include "../../kratos/utilities/binbased_fast_point_locator.h"
#include "../../kratos/utilities/normal_calculation_utils.h"
#include "../../kratos/spaces/ublas_space.h"
#include "shape_optimization_application.h"
#include "filter_function.h"
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

class MapperVertexMorphingIterative
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // Type definitions for linear algebra including sparse systems
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;    

    /// Pointer definition of MapperVertexMorphingIterative
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingIterative);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingIterative( ModelPart& designSurface, boost::python::dict dampingRegions, Parameters& optimizationSettings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables(designSurface.Nodes().size()),
          mFilterFunction( optimizationSettings["design_variables"]["filter"]["filter_function_type"].GetString() ),
          mFilterRadius( optimizationSettings["design_variables"]["filter"]["filter_radius"].GetDouble() ),
          mPerformDamping( optimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool() ),
          mpSearchTree()
    {
        createListOfNodesOfDesignSurface();
        assignMappingMatrixIds();
        initalizeDampingFactorsToHaveNoInfluence(); 
        if(mPerformDamping)   
            setDampingFactorsForAllDampingRegions( dampingRegions, optimizationSettings["design_variables"]["damping"]["damping_regions"] );         
    }

    /// Destructor.
    virtual ~MapperVertexMorphingIterative()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void createListOfNodesOfDesignSurface()
    {
        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesOfDesignSurface.push_back(pnode);
        }
    }
    
    // --------------------------------------------------------------------------
    void assignMappingMatrixIds()
    {
        unsigned int i = 0;
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            node_i->SetValue(MAPPING_MATRIX_ID,i++);
    }

    // --------------------------------------------------------------------------
    void initalizeDampingFactorsToHaveNoInfluence()
    {
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            node_i->SetValue(DAMPING_FACTOR_X,1.0);    
            node_i->SetValue(DAMPING_FACTOR_Y,1.0);  
            node_i->SetValue(DAMPING_FACTOR_Z,1.0);  
        } 
    }

    // --------------------------------------------------------------------------
    void setDampingFactorsForAllDampingRegions( boost::python::dict dampingRegions, Parameters dampingSettings )
    {
        std::cout << "\n> Starting to prepare damping..." << std::endl;

        resetSearchTreeIfAlreadyExisting( mpSearchTree );
        mpSearchTree = createSearchTreeWithAllNodesOnDesignSurface();

        // Loop over all regions for which damping is to be applied
        for (unsigned int regionNumber = 0; regionNumber < len(dampingRegions); regionNumber++)
        {
            // Extract sub-model part for damping
            std::string dampingRegionSubModelPartName = dampingSettings[regionNumber]["sub_model_part_name"].GetString();
            ModelPart& dampingRegion = extractModelPart( dampingRegions[dampingRegionSubModelPartName] );

            // Read settings for current edge
            bool dampX = dampingSettings[regionNumber]["damp_X"].GetBool();
            bool dampY = dampingSettings[regionNumber]["damp_Y"].GetBool();
            bool dampZ = dampingSettings[regionNumber]["damp_Z"].GetBool();
            std::string dampingFunctionType = dampingSettings[regionNumber]["damping_function_type"].GetString();
            double dampingRadius = dampingSettings[regionNumber]["damping_radius"].GetDouble();

            // Prepare damping function
            DampingFunction dampingFunction( dampingFunctionType, dampingRadius );

            // Loop over all nodes in specified damping sub-model part 
            for ( ModelPart::NodeIterator node_itr = dampingRegion.NodesBegin(); node_itr != dampingRegion.NodesEnd(); ++node_itr )
            {
                // Get node information
                ModelPart::NodeType& currentDampingNode = *node_itr;
                array_3d i_coord = currentDampingNode.Coordinates();

                // Arrays needed for spatial search
                unsigned int maxNodesAffected = 10000;
                NodeVector nodesAffected(maxNodesAffected);
                DoubleVector resulting_squared_distances(maxNodesAffected);

                // Perform spatial search for current node
                unsigned int numberOfNodesAffected;
                numberOfNodesAffected = mpSearchTree->SearchInRadius( currentDampingNode, 
                                                                      dampingRadius, 
                                                                      nodesAffected.begin(),
                                                                      resulting_squared_distances.begin(), 
                                                                      maxNodesAffected ); 
                
                // Throw a warning if specified (hard-coded) maximum number of neighbors is reached
                if(numberOfNodesAffected == maxNodesAffected)
                    std::cout << "\n> WARNING!!!!! For node " << currentDampingNode.Id() << " and specified damping radius, maximum number of neighbor nodes (=" << maxNodesAffected << " nodes) reached!" << std::endl;

                // Loop over all nodes in radius (including node on damping region itself)
                for(unsigned int j_itr = 0 ; j_itr<numberOfNodesAffected ; j_itr++)
                {
                    // Get node information
                    ModelPart::NodeType& node_j = *nodesAffected[j_itr];
                    array_3d j_coord = node_j.Coordinates();

                    // Computation of damping factor
                    double dampingFactor = dampingFunction.compute_damping_factor(i_coord,j_coord);

                    // For every specified damping direction we check if new damping factor is smaller than the assigned one for current node. 
                    // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
                    auto& damping_factor_variable = node_j.GetValue(DAMPING_FACTOR);
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
    void map_to_design_space( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    {
        if(mPerformDamping)
            damp_nodal_variable( rNodalVariable );
        perform_mapping_to_design_space( rNodalVariable, rNodalVariableInDesignSpace );
    }

    // --------------------------------------------------------------------------
    void map_to_geometry_space( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    {
        perform_mapping_to_geometry_space( rNodalVariable, rNodalVariableInGeometrySpace );
        if(mPerformDamping)
            damp_nodal_variable( rNodalVariable );
    }

    // --------------------------------------------------------------------------
    void perform_mapping_to_design_space( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    {
        // Measure time
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to design space..." << std::endl;

        // Perform mapping of objective sensitivities
        VectorType dJdsx, dJdsy, dJdsz;
        dJdsx.resize(mNumberOfDesignVariables);
        dJdsy.resize(mNumberOfDesignVariables);
        dJdsz.resize(mNumberOfDesignVariables);       
        
        resetSearchTreeIfAlreadyExisting( mpSearchTree );
        mpSearchTree = createSearchTreeWithAllNodesOnDesignSurface();

        // Prepare Weighting function to be used in the mapping
        FilterFunction filterFunction( mFilterFunction, mFilterRadius );

        // Loop over all design variables
        for (ModelPart::NodeIterator node_itr = mrDesignSurface.NodesBegin(); node_itr != mrDesignSurface.NodesEnd(); ++node_itr)
        {
            // Get node information
            ModelPart::NodeType& node_i = *node_itr;
            array_3d i_coord = node_i.Coordinates();

            // Arrays needed for spatial search
            unsigned int maxNodesAffected = 10000;            
            NodeVector nodesAffected(maxNodesAffected);
            DoubleVector resulting_squared_distances(maxNodesAffected);


            // Perform spatial search for current node
            unsigned int numberOfNodesAffected;
            numberOfNodesAffected = mpSearchTree->SearchInRadius( node_i, 
                                                                mFilterRadius, 
                                                                nodesAffected.begin(),
                                                                resulting_squared_distances.begin(), 
                                                                maxNodesAffected );

            // Throw a warning if specified (hard-coded) maximum number of neighbors is reached
            if(numberOfNodesAffected >= maxNodesAffected)
                std::cout << "\n> WARNING!!!!! For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << maxNodesAffected << " nodes) reached!" << std::endl;

            // Some lists to increase efficiency in the loop later
            DoubleVector list_of_weights(numberOfNodesAffected,0.0);
            std::vector<int> listOfNeighborMappingIds(numberOfNodesAffected,0);

            // Compute and assign weights in the mapping matrix
            double sum_weights = 0.0;
            for(unsigned int j_itr = 0 ; j_itr<numberOfNodesAffected ; j_itr++)
            {
                // Get node information
                ModelPart::NodeType& node_j = *nodesAffected[j_itr];
                array_3d j_coord = node_j.Coordinates();

                // Computation of weight according specified weighting function
                double weightForThisNeighbourNode = filterFunction.compute_weight(i_coord,j_coord);

                // Get the id of the node in the mapping matrix
                int j = node_j.GetValue(MAPPING_MATRIX_ID);

                // Add values to list
                listOfNeighborMappingIds[j_itr] = j;
                list_of_weights[j_itr] = weightForThisNeighbourNode;

                // Computed for integration of weighting function later using post-scaling
                sum_weights += weightForThisNeighbourNode;
            }


            // Get objective sensitivies for node_i 
            // (dJdX = node_i_obj_sens[0], dJdY = node_i_obj_sens[1], dJdZ = node_i_obj_sens[2])
            array_3d& node_i_obj_sens = node_i.FastGetSolutionStepValue(rNodalVariable);
            
            // Post scaling + sort in all matrix entries in mapping matrix
            // We sort in row by row using push_back. This is much more efficient than having only one loop and using a direct access
            for(unsigned int j_itr = 0 ; j_itr<numberOfNodesAffected ; j_itr++)
            {
                int j = listOfNeighborMappingIds[j_itr];
                double weightForThisNeighbourNode = list_of_weights[j_itr] / sum_weights;

                // Mapping of objective sensitivies
                dJdsx[j] += weightForThisNeighbourNode*node_i_obj_sens[0];
                dJdsy[j] += weightForThisNeighbourNode*node_i_obj_sens[1];
                dJdsz[j] += weightForThisNeighbourNode*node_i_obj_sens[2];
            }                          
        }

        // Assign results to nodal variables
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            VectorType dJds_i = ZeroVector(3);
            dJds_i(0) = dJdsx[i];
            dJds_i(1) = dJdsy[i];
            dJds_i(2) = dJdsz[i];
            node_i->FastGetSolutionStepValue(rNodalVariableInDesignSpace) = dJds_i;
        }

        // Console output for information
        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void perform_mapping_to_geometry_space( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    {
        // Measure time of mapping
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

        // Initialize vectors dx, dy, dz
        VectorType dx, dy, dz;
        dx.resize(mNumberOfDesignVariables);
        dy.resize(mNumberOfDesignVariables);
        dz.resize(mNumberOfDesignVariables);

        resetSearchTreeIfAlreadyExisting( mpSearchTree );
        mpSearchTree = createSearchTreeWithAllNodesOnDesignSurface();

        // Prepare Weighting function to be used in the mapping
        FilterFunction filterFunction( mFilterFunction, mFilterRadius );

        // Loop over all design variables
        for (ModelPart::NodeIterator node_itr = mrDesignSurface.NodesBegin(); node_itr != mrDesignSurface.NodesEnd(); ++node_itr)
        {
            // Get node information
            ModelPart::NodeType& node_i = *node_itr;
            array_3d i_coord = node_i.Coordinates();

            // Get tID of the node in the mapping matrix
            int i = node_i.GetValue(MAPPING_MATRIX_ID);

            // Arrays needed for spatial search
            unsigned int maxNodesAffected = 10000;            
            NodeVector nodesAffected(maxNodesAffected);
            DoubleVector resulting_squared_distances(maxNodesAffected);


            // Perform spatial search for current node
            unsigned int numberOfNodesAffected;
            numberOfNodesAffected = mpSearchTree->SearchInRadius( node_i, 
                                                                  mFilterRadius, 
                                                                  nodesAffected.begin(),
                                                                  resulting_squared_distances.begin(), 
                                                                  maxNodesAffected );

            // Throw a warning if specified (hard-coded) maximum number of neighbors is reached
            if(numberOfNodesAffected >= maxNodesAffected)
                std::cout << "\n> WARNING!!!!! For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << maxNodesAffected << " nodes) reached!" << std::endl;

            // Some list to increase efficiency in the loop later
            DoubleVector list_of_weights(numberOfNodesAffected,0.0);
            std::vector<int> listOfNeighborMappingIds(numberOfNodesAffected,0);

            // Compute and assign weights in the mapping matrix
            double sum_weights = 0.0;
            for(unsigned int j_itr = 0 ; j_itr<numberOfNodesAffected ; j_itr++)
            {
                // Get node information
                ModelPart::NodeType& node_j = *nodesAffected[j_itr];
                array_3d j_coord = node_j.Coordinates();

                // Get the id of the node in the mapping matrix
                int j = node_j.GetValue(MAPPING_MATRIX_ID);

                // Computation of weight according specified weighting function
                double weightForThisNeighbourNode = filterFunction.compute_weight(i_coord,j_coord);

                // Add values to list
                listOfNeighborMappingIds[j_itr] = j;
                list_of_weights[j_itr] = weightForThisNeighbourNode;

                // Computed for integration of weighting function later using post-scaling
                sum_weights += weightForThisNeighbourNode;
            }

            // Post scaling + sort in all matrix entries in mapping matrix
            // We sort in row by row using push_back. This is much more efficient than having only one loop and using a direct access
            for(unsigned int j_itr = 0 ; j_itr<numberOfNodesAffected ; j_itr++)
            {
                double weightForThisNeighbourNode = list_of_weights[j_itr] / sum_weights;

                // Get design update from node_j
                ModelPart::NodeType& node_j = *nodesAffected[j_itr];
                array_3d& node_design_update = node_j.FastGetSolutionStepValue(rNodalVariable);

                // add design update of node_j into vector of shape update
                dx[i] += weightForThisNeighbourNode*node_design_update[0];
                dy[i] += weightForThisNeighbourNode*node_design_update[1];
                dz[i] += weightForThisNeighbourNode*node_design_update[2];

            }                          
        }

        // Assign dx as nodal shape updates and update coordinates accordingly
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            // Read shape update from solution vector
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            array_3d shape_update;
            shape_update[0] = dx[i];
            shape_update[1] = dy[i];
            shape_update[2] = dz[i];

            // Assign shape update to nodal variable
            noalias(node_i->FastGetSolutionStepValue(rNodalVariableInGeometrySpace)) = shape_update;
        }

        // Console output for information
        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void damp_nodal_variable( const Variable<array_3d> &rNodalVariable )
    {
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {   
            auto& damping_factor = node_i->GetValue(DAMPING_FACTOR);
            auto& nodalVariable = node_i->FastGetSolutionStepValue(rNodalVariable);

            nodalVariable[0] *= damping_factor[0];
            nodalVariable[1] *= damping_factor[1];
            nodalVariable[2] *= damping_factor[2];  
        }  
    }

    // --------------------------------------------------------------------------
    void resetSearchTreeIfAlreadyExisting( KDTree::Pointer searchTree )
    {
        if(searchTree != 0)
            searchTree.reset();
    }

    // --------------------------------------------------------------------------
    KDTree::Pointer createSearchTreeWithAllNodesOnDesignSurface()
    {
        return boost::shared_ptr<KDTree>(new KDTree(mListOfNodesOfDesignSurface.begin(), mListOfNodesOfDesignSurface.end(), mBucketSize));
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
        return "MapperVertexMorphingIterative";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperVertexMorphingIterative";
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
    ModelPart& mrDesignSurface;
    const unsigned int mNumberOfDesignVariables;
    std::string mFilterFunction;
    const double mFilterRadius;
    bool mPerformDamping;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesOfDesignSurface;
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
//      MapperVertexMorphingIterative& operator=(MapperVertexMorphingIterative const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingIterative(MapperVertexMorphingIterative const& rOther);


    ///@}

}; // Class MapperVertexMorphingIterative

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_ITERATIVE_H
