// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef VERTEX_MORPHING_MAPPER_H
#define VERTEX_MORPHING_MAPPER_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

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

class VertexMorphingMapper
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

    /// Pointer definition of VertexMorphingMapper
    KRATOS_CLASS_POINTER_DEFINITION(VertexMorphingMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VertexMorphingMapper( ModelPart& designSurface, boost::python::dict dampingRegions, Parameters& optimizationSettings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables(designSurface.Nodes().size()),
          mFilterFunction( optimizationSettings["design_variables"]["filter"]["filter_function_type"].GetString() ),
          mFilterRadius( optimizationSettings["design_variables"]["filter"]["filter_radius"].GetDouble() ),
          mPerformDamping( optimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool() )
    {
        setPrecisionForOutput();
        assignMappingMatrixIds();
        initalizeDampingFactorsToHaveNoInfluence(); 
        if(mPerformDamping)   
            setDampingFactorsForAllDampingRegions( dampingRegions, optimizationSettings["design_variables"]["damping"]["damping_regions"] );
    }

    /// Destructor.
    virtual ~VertexMorphingMapper()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void setPrecisionForOutput()
    {
        std::cout.precision(12);
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

        // Creating an auxiliary list for the nodes to be searched on
        NodeVector listOfNodes;

        // Add nodes to list wich collects all nodes for neighbour-search
        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            listOfNodes.push_back(pnode);
        }

        // Compute tree with the node positions
        unsigned int bucket_size = 100;
        KDTree treeOfAllOptimizationNodes(listOfNodes.begin(), listOfNodes.end(), bucket_size);

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
                numberOfNodesAffected = treeOfAllOptimizationNodes.SearchInRadius( currentDampingNode, 
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
                    if(dampX == true && dampingFactor < node_j.GetValue(DAMPING_FACTOR_X))
                        node_j.SetValue(DAMPING_FACTOR_X, dampingFactor);   

                    if(dampY == true && dampingFactor < node_j.GetValue(DAMPING_FACTOR_Y))       
                        node_j.SetValue(DAMPING_FACTOR_Y, dampingFactor);   

                    if(dampZ == true && dampingFactor < node_j.GetValue(DAMPING_FACTOR_Z))       
                        node_j.SetValue(DAMPING_FACTOR_Z, dampingFactor);                            
                }
            }
        }

        std::cout << "> Finished preparation of damping." << std::endl;
    }

    // --------------------------------------------------------------------------
    void compute_mapping_matrix()
    {
        // Measure time
        boost::timer mapping_time;

        std::cout << "\n> Starting computation of complete mapping matrix..." << std::endl;

        // Initialize filter matrix
        mMappingMatrix.resize(mNumberOfDesignVariables,mNumberOfDesignVariables);

        // Initialize mapping matrix (zero possible entries)
        mMappingMatrix.clear();

        // Creating an auxiliary list for the nodes to be searched on
        NodeVector listOfNodes;

        // starting calculating time of construction of the kdtree
        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());

            // Putting the nodes of interest in an auxiliary list
            listOfNodes.push_back(pnode);
        }

        // Compute tree with the node positions
        unsigned int bucket_size = 100;
        KDTree treeOfAllOptimizationNodes(listOfNodes.begin(), listOfNodes.end(), bucket_size);

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
            numberOfNodesAffected = treeOfAllOptimizationNodes.SearchInRadius( node_i, 
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

                // Get the id of the node in the mapping matrix
                int j = node_j.GetValue(MAPPING_MATRIX_ID);

                // Computation of weight according specified weighting function
                double Aij = filterFunction.compute_weight(i_coord,j_coord);

                // Add values to list
                listOfNeighborMappingIds[j_itr] = j;
                list_of_weights[j_itr] = Aij;

                // Computed for integration of weighting function later using post-scaling
                sum_weights += Aij;
            }

            // Post scaling + sort in all matrix entries in mapping matrix
            // We sort in row by row using push_back. This is much more efficient than having only one loop and using a direct access
            for(unsigned int j_itr = 0 ; j_itr<numberOfNodesAffected ; j_itr++)
            {
                int j = listOfNeighborMappingIds[j_itr];
                double Aij = list_of_weights[j_itr] / sum_weights;
                mMappingMatrix.push_back(i,j,Aij);               
            }                          
        }

        // Console output for information
        std::cout << "> Time needed for computation of mapping matrix: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void map_sensitivities_to_design_space( bool constraint_given )
    {
        // Measure time
        boost::timer mapping_time;
        std::cout << "\n> Starting mapping of sensitivities to design space..." << std::endl;

        // First we apply edge damping to the sensitivities if specified
        if(mPerformDamping)
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY_X) *= node_i->GetValue(DAMPING_FACTOR_X);
                node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY_Y) *= node_i->GetValue(DAMPING_FACTOR_Y);
                node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY_Z) *= node_i->GetValue(DAMPING_FACTOR_Z);
                if(constraint_given)
                {
                    node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY_X) *= node_i->GetValue(DAMPING_FACTOR_X);
                    node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY_Y) *= node_i->GetValue(DAMPING_FACTOR_Y);
                    node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY_Z) *= node_i->GetValue(DAMPING_FACTOR_Z);                    
                }                
            }

        // Map objective sensitivities

        // Assign nodal sensitivities to vectors used for mapping
        VectorType dJdX, dJdY, dJdZ;
        dJdX.resize(mNumberOfDesignVariables);
        dJdY.resize(mNumberOfDesignVariables);
        dJdZ.resize(mNumberOfDesignVariables);
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            array_3d node_sens = node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            dJdX[i] = node_sens[0];
            dJdY[i] = node_sens[1];
            dJdZ[i] = node_sens[2];
        }
        VectorType dCdX, dCdY, dCdZ;
        if(constraint_given)
        {
            dCdX.resize(mNumberOfDesignVariables);
            dCdY.resize(mNumberOfDesignVariables);
            dCdZ.resize(mNumberOfDesignVariables);
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d node_sens = node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                dCdX[i] = node_sens[0];
                dCdY[i] = node_sens[1];
                dCdZ[i] = node_sens[2];
            }
        }

        // Perform mapping of objective sensitivities
        VectorType dJdsx, dJdsy, dJdsz;
        dJdsx.resize(mNumberOfDesignVariables);
        dJdsy.resize(mNumberOfDesignVariables);
        dJdsz.resize(mNumberOfDesignVariables);
        SparseSpaceType::TransposeMult(mMappingMatrix,dJdX,dJdsx);
        SparseSpaceType::TransposeMult(mMappingMatrix,dJdY,dJdsy);
        SparseSpaceType::TransposeMult(mMappingMatrix,dJdZ,dJdsz);

        // Assign results to nodal variables
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            VectorType dJds_i = ZeroVector(3);
            dJds_i(0) = dJdsx[i];
            dJds_i(1) = dJdsy[i];
            dJds_i(2) = dJdsz[i];
            node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY) = dJds_i;
        }

        // Repeat mapping to map constraint sensitivities
        if(constraint_given)
        {
            // Perform mapping
            VectorType dCdsx, dCdsy, dCdsz;
            dCdsx.resize(mNumberOfDesignVariables);
            dCdsy.resize(mNumberOfDesignVariables);
            dCdsz.resize(mNumberOfDesignVariables);
            SparseSpaceType::TransposeMult(mMappingMatrix,dCdX,dCdsx);
            SparseSpaceType::TransposeMult(mMappingMatrix,dCdY,dCdsy);
            SparseSpaceType::TransposeMult(mMappingMatrix,dCdZ,dCdsz);

            // Assign results to nodal variables
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                VectorType dCds_i = ZeroVector(3);
                dCds_i(0) = dCdsx[i];
                dCds_i(1) = dCdsy[i];
                dCds_i(2) = dCdsz[i];
                node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY) = dCds_i;
            }
        }

        // Console output for information
        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void map_design_update_to_geometry_space()
    {
        // Measure time of mapping
        boost::timer mapping_time;
        std::cout << "\n> Starting mapping of design update to geometry space..." << std::endl;

        // Assign design update to vector which shall be mapped (depending specified mapping matrix )
        VectorType dsx, dsy, dsz;
        dsx.resize(mNumberOfDesignVariables);
        dsy.resize(mNumberOfDesignVariables);
        dsz.resize(mNumberOfDesignVariables);
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            dsx[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_X);
            dsy[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Y);
            dsz[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Z);
        }

        // Perform mapping to compute shape update
        VectorType dx, dy, dz;
        dx.resize(mNumberOfDesignVariables);
        dy.resize(mNumberOfDesignVariables);
        dz.resize(mNumberOfDesignVariables);
        noalias(dx) = prod(mMappingMatrix,dsx);
        noalias(dy) = prod(mMappingMatrix,dsy);
        noalias(dz) = prod(mMappingMatrix,dsz);

        // Assign dx as nodal shape updates and update coordinates accordingly
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            // If shape update deactivated, set it to zero
            if(node_i->FastGetSolutionStepValue(SHAPE_UPDATES_DEACTIVATED))
            {
                array_3d zero_array(3,0.0);
                noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = zero_array;
            }
            // In case it is not deactivated, it is checked if it is on a specified boundary beyond which no update is wanted
            else
            {
                // Read shape update from solution vector
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d shape_update;
                shape_update[0] = dx[i];
                shape_update[1] = dy[i];
                shape_update[2] = dz[i];

                // If node is on specified boundary, project shape update on specified boundary plane
                // I.e. remove component that is normal to the boundary plane
                if(node_i->FastGetSolutionStepValue(IS_ON_BOUNDARY))
                {
                    array_3d boundary_plane = node_i->FastGetSolutionStepValue(BOUNDARY_PLANE);
                    shape_update = shape_update - (inner_prod(shape_update,boundary_plane))*boundary_plane/norm_2(boundary_plane);
                }

                // Assign shape update to nodal variable
                noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = shape_update;
            }
        }

        // Apply edge damping if specified
        if(mPerformDamping)
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {   
                node_i->FastGetSolutionStepValue(SHAPE_UPDATE_X) *= node_i->GetValue(DAMPING_FACTOR_X);
                node_i->FastGetSolutionStepValue(SHAPE_UPDATE_Y) *= node_i->GetValue(DAMPING_FACTOR_Y);
                node_i->FastGetSolutionStepValue(SHAPE_UPDATE_Z) *= node_i->GetValue(DAMPING_FACTOR_Z);
            }

        // Update model part and record absolute shape change 
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            array_3d shape_update = node_i->FastGetSolutionStepValue(SHAPE_UPDATE);
                
            // Update coordinates
            node_i->X() += shape_update[0];
            node_i->Y() += shape_update[1];
            node_i->Z() += shape_update[2];

            // Add shape update to previous updates
            noalias(node_i->FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE)) += shape_update;
        }

        // Console output for information
        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
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
        return "VertexMorphingMapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "VertexMorphingMapper";
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
    // General working arrays
    // ==============================================================================
    SparseMatrixType mMappingMatrix;

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
//      VertexMorphingMapper& operator=(VertexMorphingMapper const& rOther);

    /// Copy constructor.
//      VertexMorphingMapper(VertexMorphingMapper const& rOther);


    ///@}

}; // Class VertexMorphingMapper

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // VERTEX_MORPHING_MAPPER_H
