// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
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
    typedef boost::python::extract<double> extractDouble;
    typedef boost::python::extract<std::string> extractString;
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

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // ==========================================================================
    // Type definitions for linear algebra including sparse systems
    // ==========================================================================
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;
    typedef typename SparseSpaceType::VectorType VectorType;

    /// Pointer definition of VertexMorphingMapper
    KRATOS_CLASS_POINTER_DEFINITION(VertexMorphingMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VertexMorphingMapper( ModelPart& model_part,
                          std::string filter_function_type,
                          double filter_size,
                          bool perform_damping,
                          boost::python::list damping_regions )
        : mr_opt_model_part(model_part),
          m_filter_function_type(filter_function_type),
          m_filter_size(filter_size),
          m_perform_damping(perform_damping),
          m_damping_regions(damping_regions),
          m_number_of_design_variables(model_part.Nodes().size())
    {
        // Set precision for output
        std::cout.precision(12);

        assignMappingMatrixIds();

        initalizeDampingFactorsToHaveNoInfluence(); 
        if(perform_damping)   
            setDampingFactorsAccordingDampingFunction();
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
    void assignMappingMatrixIds()
    {
        unsigned int i = 0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            node_i->SetValue(MAPPING_MATRIX_ID,i++);
    }

    // --------------------------------------------------------------------------
    void initalizeDampingFactorsToHaveNoInfluence()
    {
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            node_i->SetValue(DAMPING_FACTOR_X,1.0);    
            node_i->SetValue(DAMPING_FACTOR_Y,1.0);  
            node_i->SetValue(DAMPING_FACTOR_Z,1.0);  
        } 
    }

    // --------------------------------------------------------------------------
    void setDampingFactorsAccordingDampingFunction()
    {
        std::cout << "\n> Starting to prepare damping..." << std::endl;

        // Creating an auxiliary list for the nodes to be searched on
        NodeVector list_of_nodes;

        // Start constructing and computing the kdtree
        typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;

        // Add nodes to list wich collects all nodes for neighbour-search
        for (ModelPart::NodesContainerType::iterator node_it = mr_opt_model_part.NodesBegin(); node_it != mr_opt_model_part.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            list_of_nodes.push_back(pnode);
        }

        // Arrays needed for spatial search
        unsigned int max_nodes_affected = 10000;
        NodeVector nodes_affected(max_nodes_affected);
        DoubleVector resulting_squared_distances(max_nodes_affected);

        // Compute tree with the node positions
        unsigned int bucket_size = 100;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

        // Loop over all regions for which damping is to be applied
        for (unsigned int region_itr = 0; region_itr < boost::python::len(m_damping_regions); region_itr++)
        {
            // Extract sub-model part for damping
            ModelPart& damping_sub_mdpa = extractModelPart(m_damping_regions[region_itr][0]);

            // Read settings for current edge
            bool damp_in_X = m_damping_regions[region_itr][1];
            bool damp_in_Y = m_damping_regions[region_itr][2];
            bool damp_in_Z = m_damping_regions[region_itr][3];
            std::string damping_function_type = extractString(m_damping_regions[region_itr][4]);
            double damping_radius = extractDouble(m_damping_regions[region_itr][5]);

            // Prepare damping function
            DampingFunction damping_function( damping_function_type, damping_radius );

            // Loop over all nodes in specified damping sub-model part 
            for ( ModelPart::NodeIterator node_itr = damping_sub_mdpa.NodesBegin(); node_itr != damping_sub_mdpa.NodesEnd(); ++node_itr )
            {
                // Get node information
                ModelPart::NodeType& damping_node_i = *node_itr;
                array_3d i_coord = damping_node_i.Coordinates();

                // Perform spatial search for current node
                unsigned int number_of_nodes_affected;
                number_of_nodes_affected = nodes_tree.SearchInRadius( damping_node_i, 
                                                                      damping_radius, 
                                                                      nodes_affected.begin(),
                                                                      resulting_squared_distances.begin(), 
                                                                      max_nodes_affected ); 
                
                // Throw a warning if specified (hard-coded) maximum number of neighbors is reached
                if(number_of_nodes_affected == max_nodes_affected)
                    std::cout << "\n> WARNING!!!!! For node " << damping_node_i.Id() << " and specified damping radius, maximum number of neighbor nodes (=" << max_nodes_affected << " nodes) reached!" << std::endl;

                // Loop over all nodes in radius (including node on damping region itself)
                for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
                {
                    // Get node information
                    ModelPart::NodeType& node_j = *nodes_affected[j_itr];
                    array_3d j_coord = node_j.Coordinates();

                    // Computation of damping factor
                    double damping_factor = damping_function.compute_damping_factor(i_coord,j_coord);

                    // For every specified damping direction we check if new damping factor is smaller than the assigned one for current node. 
                    // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
                    if(damp_in_X == true and damping_factor < node_j.GetValue(DAMPING_FACTOR_X))
                        node_j.SetValue(DAMPING_FACTOR_X, damping_factor);     
                    if(damp_in_Y == true and damping_factor < node_j.GetValue(DAMPING_FACTOR_Y))       
                        node_j.SetValue(DAMPING_FACTOR_Y, damping_factor);   
                    if(damp_in_Z == true and damping_factor < node_j.GetValue(DAMPING_FACTOR_Z))       
                        node_j.SetValue(DAMPING_FACTOR_Z, damping_factor);                            
                }
            }
        }

        std::cout << "\n> Finished preparation of damping." << std::endl;
    }

    // --------------------------------------------------------------------------
    void compute_mapping_matrix()
    {
        KRATOS_TRY;

        // Measure time
        boost::timer mapping_time;

        std::cout << "\n> Starting computation of complete mapping matrix..." << std::endl;

        // Initialize filter matrix
        m_mapping_matrix.resize(m_number_of_design_variables,m_number_of_design_variables);

        // Initialize mapping matrix (zero possible entries)
        m_mapping_matrix.clear();

        // Creating an auxiliary list for the nodes to be searched on
        NodeVector list_of_nodes;

        // Start constructing and computing the kdtree
        typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;

        // starting calculating time of construction of the kdtree
        for (ModelPart::NodesContainerType::iterator node_it = mr_opt_model_part.NodesBegin(); node_it != mr_opt_model_part.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());

            // Putting the nodes of interest in an auxiliary list
            list_of_nodes.push_back(pnode);
        }

        // Arrays needed for spatial search
        unsigned int max_nodes_affected = 10000;            
        NodeVector nodes_affected(max_nodes_affected);
        DoubleVector resulting_squared_distances(max_nodes_affected);

        // Compute tree with the node positions
        unsigned int bucket_size = 100;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

        // Prepare Weighting function to be used in the mapping
        FilterFunction filter_function( m_filter_function_type, m_filter_size );

        // Loop over all design variables
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            // Initialize list of affected nodes
            nodes_affected.clear();

            // Get node information
            ModelPart::NodeType& node_i = *node_itr;
            array_3d i_coord = node_i.Coordinates();

            // Get tID of the node in the mapping matrix
            int i = node_i.GetValue(MAPPING_MATRIX_ID);

            // Perform spatial search for current node
            unsigned int number_of_nodes_affected;
            number_of_nodes_affected = nodes_tree.SearchInRadius( node_i, 
                                                                  m_filter_size, 
                                                                  nodes_affected.begin(),
                                                                  resulting_squared_distances.begin(), 
                                                                  max_nodes_affected );

            // Throw a warning if specified (hard-coded) maximum number of neighbors is reached
            if(number_of_nodes_affected >= max_nodes_affected)
                std::cout << "\n> WARNING!!!!! For node " << node_i.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << max_nodes_affected << " nodes) reached!" << std::endl;

            // Some lists to increase efficiency in the loop later
            DoubleVector list_of_weights(number_of_nodes_affected,0.0);
            std::vector<int> list_of_neighbor_mapping_ids(number_of_nodes_affected,0);

            // Compute and assign weights in the mapping matrix
            double sum_weights = 0.0;
            for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
            {
                // Get node information
                ModelPart::NodeType& node_j = *nodes_affected[j_itr];
                array_3d j_coord = node_j.Coordinates();

                // Get the id of the node in the mapping matrix
                int j = node_j.GetValue(MAPPING_MATRIX_ID);

                // Computation of weight according specified weighting function
                double Aij = filter_function.compute_weight(i_coord,j_coord);

                // Add values to list
                list_of_neighbor_mapping_ids[j_itr] = j;
                list_of_weights[j_itr] = Aij;

                // Computed for integration of weighting function later using post-scaling
                sum_weights += Aij;
            }

            // Post scaling + sort in all matrix entries in mapping matrix
            // We sort in row by row using push_back. This is much more efficient than having only one loop and using a direct access
            for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
            {
                int j = list_of_neighbor_mapping_ids[j_itr];
                double Aij = list_of_weights[j_itr] / sum_weights;
                m_mapping_matrix.push_back(i,j,Aij);               
            }                          
        }

        // Console output for information
        std::cout << "> Time needed for computation of mapping matrix: " << mapping_time.elapsed() << " s" << std::endl;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_sensitivities_to_design_space( bool constraint_given )
    {
        KRATOS_TRY;

        // Measure time
        boost::timer mapping_time;
        std::cout << "\n> Starting mapping of sensitivities to design space..." << std::endl;

        // First we apply edge damping to the sensitivities if specified
        if(m_perform_damping)
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
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

        // Assign nodal sensitivities to vector used for mapping
        VectorType dJdX;
        VectorType dJdY;
        VectorType dJdZ;
        dJdX.resize(m_number_of_design_variables);
        dJdY.resize(m_number_of_design_variables);
        dJdZ.resize(m_number_of_design_variables);
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            array_3d node_sens = node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            dJdX[i] = node_sens[0];
            dJdY[i] = node_sens[1];
            dJdZ[i] = node_sens[2];
        }
        VectorType dCdX;
        VectorType dCdY;
        VectorType dCdZ;
        if(constraint_given)
        {
            dCdX.resize(m_number_of_design_variables);
            dCdY.resize(m_number_of_design_variables);
            dCdZ.resize(m_number_of_design_variables);
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d node_sens = node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                dCdX[i] = node_sens[0];
                dCdY[i] = node_sens[1];
                dCdZ[i] = node_sens[2];
            }
        }

        VectorType dJdsx;
        VectorType dJdsy;
        VectorType dJdsz;
        dJdsx.resize(m_number_of_design_variables);
        dJdsy.resize(m_number_of_design_variables);
        dJdsz.resize(m_number_of_design_variables);
        SparseSpaceType::TransposeMult(m_mapping_matrix,dJdX,dJdsx);
        SparseSpaceType::TransposeMult(m_mapping_matrix,dJdY,dJdsy);
        SparseSpaceType::TransposeMult(m_mapping_matrix,dJdZ,dJdsz);

        // Assign results to nodal variables
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
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
            VectorType dCdsx;
            VectorType dCdsy;
            VectorType dCdsz;
            dCdsx.resize(m_number_of_design_variables);
            dCdsy.resize(m_number_of_design_variables);
            dCdsz.resize(m_number_of_design_variables);
            SparseSpaceType::TransposeMult(m_mapping_matrix,dCdX,dCdsx);
            SparseSpaceType::TransposeMult(m_mapping_matrix,dCdY,dCdsy);
            SparseSpaceType::TransposeMult(m_mapping_matrix,dCdZ,dCdsz);

            // Assign results to nodal variables
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
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

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_design_update_to_geometry_space()
    {
        KRATOS_TRY;

        // Measure time of mapping
        boost::timer mapping_time;
        std::cout << "\n> Starting mapping of design update to geometry space..." << std::endl;

        // Assign design update to vector which shall be mapped (depending specified mapping matrix )
        VectorType dsx;
        VectorType dsy;
        VectorType dsz;
        dsx.resize(m_number_of_design_variables);
        dsy.resize(m_number_of_design_variables);
        dsz.resize(m_number_of_design_variables);
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            dsx[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_X);
            dsy[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Y);
            dsz[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Z);
        }

        // Perform mapping to compute shape update
        VectorType dx;
        VectorType dy;
        VectorType dz;
        dx.resize(m_number_of_design_variables);
        dy.resize(m_number_of_design_variables);
        dz.resize(m_number_of_design_variables);
        noalias(dx) = prod(m_mapping_matrix,dsx);
        noalias(dy) = prod(m_mapping_matrix,dsy);
        noalias(dz) = prod(m_mapping_matrix,dsz);

        // Assign dx as nodal shape updates and update coordinates accordingly
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
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
        if(m_perform_damping)
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {   
                node_i->FastGetSolutionStepValue(SHAPE_UPDATE_X) *= node_i->GetValue(DAMPING_FACTOR_X);
                node_i->FastGetSolutionStepValue(SHAPE_UPDATE_Y) *= node_i->GetValue(DAMPING_FACTOR_Y);
                node_i->FastGetSolutionStepValue(SHAPE_UPDATE_Z) *= node_i->GetValue(DAMPING_FACTOR_Z);
            }

        // Update model part and record absolute shape change 
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
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

        KRATOS_CATCH("");
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
    ModelPart& mr_opt_model_part;
    std::string m_filter_function_type;
    const double m_filter_size;
    bool m_perform_damping;
    boost::python::list m_damping_regions;
    const unsigned int m_number_of_design_variables;

    // ==============================================================================
    // General working arrays
    // ==============================================================================
    SparseMatrixType m_mapping_matrix;

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
