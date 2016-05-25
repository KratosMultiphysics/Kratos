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
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
//   Date:                $Date:                      March 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef VERTEX_MORPHING_UTILITIES_H
#define VERTEX_MORPHING_UTILITIES_H

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
#include "shape_optimization_application.h"
#include "../../kratos/spaces/ublas_space.h"

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

class VertexMorphingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > PointType;
    typedef Node < 3 > ::Pointer PointTypePointer;
    typedef std::vector<PointType::Pointer> PointVector;
    typedef std::vector<PointType::Pointer>::iterator PointIterator;
    typedef std::vector<double> DistanceVector;
    typedef std::vector<double>::iterator DistanceIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // ==========================================================================
    // Type definitions for linear algebra including sparse systems
    // ==========================================================================
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;

    /// Pointer definition of VertexMorphingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(VertexMorphingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VertexMorphingUtilities( ModelPart& model_part,
                             const int domain_size,
                             boost::python::dict py_objectives,
                             boost::python::dict py_constraints,
                             double filter_size,
                             const int max_nodes_affected )
        : mr_opt_model_part(model_part),
          m_domain_size(domain_size),
          m_filter_size(filter_size),
          m_number_of_design_variables(model_part.Nodes().size()),
          m_max_nodes_affected(max_nodes_affected)
    {
        // Set precision for output
        std::cout.precision(12);

        // Initialize filter matrix
        m_mapping_matrix.resize(m_number_of_design_variables*3,m_number_of_design_variables);

        // Create map to obtain local mapping matrix Id from global node Id
        unsigned int i = 0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            // Store local mapping matrix Id on the node
            node_i->SetValue(MAPPING_MATRIX_ID,i);

            // iterator design variable iterator i
            i++;
        }

        // Create a map between python dict string ids and C++ map integer ids

        // Create map for all objectives
        boost::python::list py_F_ids = py_objectives.keys();
        for (unsigned int i = 0; i < boost::python::len(py_F_ids); ++i)
        {
            boost::python::extract<const char*> py_F_id(py_F_ids[i]);
            py2cpp_F_id[py_F_id] = i;
        }

        // Create map for all constraints
        boost::python::list py_C_ids = py_constraints.keys();
        for (unsigned int i = 0; i < boost::python::len(py_C_ids); ++i)
        {
            boost::python::extract<const char*> py_C_id(py_C_ids[i]);
            py2cpp_C_id[py_C_id] = i;
        }
    }

    /// Destructor.
    virtual ~VertexMorphingUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    // General geometrical operations
    // ==============================================================================
    void compute_unit_surface_normals()
    {
        KRATOS_TRY;

        // Compute nodal are normal using given Kratos utilities (sets the variable "NORMAL")
        NormalCalculationUtils normal_util = NormalCalculationUtils();
        normal_util.CalculateOnSimplex(mr_opt_model_part,m_domain_size);

        // Take into account boundary conditions, normalize area normal and store in respective variable
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            // Normalize normal and assign to solution step value
            array_3d area_normal = node_i->FastGetSolutionStepValue(NORMAL);
            array_3d normalized_normal = area_normal / norm_2(area_normal);
            noalias(node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)) = normalized_normal;
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void project_grad_on_unit_surface_normal( bool constraint_given )
    {
        KRATOS_TRY;

        // We loop over all nodes and compute the part of the sensitivity which is in direction to the surface normal
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            // We compute dFdX_n = (dFdX \cdot n) * n
            array_3d node_sens = node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            array_3d node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
            double surface_sens = inner_prod(node_sens,node_normal);
            array_3d normal_node_sens = surface_sens * node_normal;

            // Assign resulting sensitivities back to node
            node_i->GetSolutionStepValue(OBJECTIVE_SURFACE_SENSITIVITY) = surface_sens;
            noalias(node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY)) = normal_node_sens;

            // Repeat for constraint
            if(constraint_given)
            {
                // We compute dFdX_n = (dFdX \cdot n) * n
                node_sens = node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
                surface_sens = inner_prod(node_sens,node_normal);
                normal_node_sens =  surface_sens * node_normal;

                // Assign resulting sensitivities back to node
                node_i->GetSolutionStepValue(CONSTRAINT_SURFACE_SENSITIVITY) = surface_sens;
                noalias(node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY)) = normal_node_sens;
            }
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For perfoming Vertex Morphing
    // ==============================================================================
    void compute_mapping_matrix()
    {
        KRATOS_TRY;

        // Measure time of mapping
        boost::timer mapping_time;
        std::cout << "> Start computing mapping matrix..." << std::endl;

        // Creating an auxiliary list for the nodes to be searched on
        PointVector list_of_nodes;

        // Start constructing and computing the kdtree
        typedef Bucket< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;

        // starting calculating time of construction of the kdtree
        for (ModelPart::NodesContainerType::iterator node_it = mr_opt_model_part.NodesBegin(); node_it != mr_opt_model_part.NodesEnd(); ++node_it)
        {
            PointTypePointer pnode = *(node_it.base());

            // Putting the nodes of interest in an auxiliary list
            list_of_nodes.push_back(pnode);
        }

        // Arrays needed for spatial search
        PointVector nodes_affected(m_max_nodes_affected);
        DistanceVector resulting_squared_distances(m_max_nodes_affected);

        // Compute tree with the node positions
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), m_max_nodes_affected);

        // Loop over all design variables
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
        	// Get node information
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];
            array_3d i_coord = node_i.Coordinates();

            // Get tID of the node in the mapping matrix
            int i = node_i.GetValue(MAPPING_MATRIX_ID);

            // Perform spatial search for current node
            unsigned int number_of_nodes_affected;
            number_of_nodes_affected = nodes_tree.SearchInRadius(node_i, m_filter_size, nodes_affected.begin(),resulting_squared_distances.begin(), m_max_nodes_affected);

            // Store results to reuse later in the forward mapping
            m_listOf_nodesAffected[i_ID] = nodes_affected;
            m_number_of_nodes_affected[i_ID] = number_of_nodes_affected;

            // Compute and assign weights in the mapping matrix
            double sum_weights = 0.0;
            for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
            {
            	// Get node information
                int j_ID = nodes_affected[j_itr]->Id();
                ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];
                array_3d j_coord(3,0.0);
                j_coord[0] = node_j.X();
                j_coord[1] = node_j.Y();
                j_coord[2] = node_j.Z();

                // Get the id of the node in the mapping matrix
                int j = node_j.GetValue(MAPPING_MATRIX_ID);

                array_3d dist_vector = i_coord - j_coord;
                double squared_scalar_distance = dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2];

                // Computation of weight according specified weighting function
                // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                double Aij = exp(-squared_scalar_distance/(2*m_filter_size*m_filter_size/9.0));

                // Multiplication by the node normal (nodal director)
                // In this way we implicitly preserve the in-plane mesh quality (pure Heuristic)
                m_mapping_matrix(3*i+0,j) = Aij*node_j.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)[0];
                m_mapping_matrix(3*i+1,j) = Aij*node_j.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)[1];
                m_mapping_matrix(3*i+2,j) = Aij*node_j.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)[2];

                // Computed for integration of weighting function later using post-scaling
                sum_weights += Aij;
            }

            // Here we perform the post scaling by the sum of all weights
            for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
            {
            	// Get node information
                int j_ID = nodes_affected[j_itr]->Id();
                ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];

                // Get the id of the node in the mapping matrix
                int j = node_j.GetValue(MAPPING_MATRIX_ID);

                // Scaling of the weights
                m_mapping_matrix(3*i+0,j) /= sum_weights;
                m_mapping_matrix(3*i+1,j) /= sum_weights;
                m_mapping_matrix(3*i+2,j) /= sum_weights;
            }
        }

        // Console output for information
        std::cout << "> Finished computing mapping matrix!" << std::endl;
        std::cout << "> Time needed for computation of mapping matrix: " << mapping_time.elapsed() << std::endl;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_sensitivities_to_design_space( bool constraint_given )
    {
    	KRATOS_TRY;

    	// Measure time of mapping
    	boost::timer mapping_time;
    	std::cout << "\n> Start mapping sensitivities to design space..." << std::endl;

    	// Loop over all design variables
    	for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
    	{
    		// Get node information
    		int i_ID = node_i->Id();

    		// Get the id of the node in the mapping matrix
    		int i = node_i->GetValue(MAPPING_MATRIX_ID);

    		//Instead of performing spatial search again, we read the results obtained from the computation of the mapping matrix
    		unsigned int number_of_nodes_affected = m_number_of_nodes_affected[i_ID];
    		PointVector nodes_affected(number_of_nodes_affected);
    		nodes_affected = m_listOf_nodesAffected[i_ID];

    		// Mapping gradients to design space - looping only over neighbor nodes
    		double dFds_i = 0.0;
    		for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
    		{
    			// Get node information
    			int j_ID = nodes_affected[j_itr]->Id();
    			ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];

    			// Get the id of the node in the mapping matrix
    			unsigned int j = node_j.GetValue(MAPPING_MATRIX_ID);

    			// Ask for the sensitiviteis
    			array_3d node_sens = node_j.FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);

    			// Compute dot-product
    			dFds_i += m_mapping_matrix(j*3+0,i) * node_sens[0]
						+ m_mapping_matrix(j*3+1,i) * node_sens[1]
						+ m_mapping_matrix(j*3+2,i) * node_sens[2];
    		}
    		m_filtered_dFdX[i_ID] = dFds_i;

    		// If constrainted optimization, then also compute filtered gradients for the constraints (dCds)
    		if(constraint_given)
    		{
    			// Mapping sensitivities to design space - looping only over neighbor nodes
    			double dCds_i = 0.0;
    			for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
    			{
    				// Get node information
    				int j_ID = nodes_affected[j_itr]->Id();
    				ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];

    				// Get the id of the node in the mapping matrix
    				unsigned int j = node_j.GetValue(MAPPING_MATRIX_ID);

    				// Ask for the sensitiviteis
    				array_3d node_sens = node_j.FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);

    				// Compute dot-product - Here we use the transpose of the mapping matrix (exchanged indices)
    				dCds_i += m_mapping_matrix(j*3+0,i) * node_sens[0]
							+ m_mapping_matrix(j*3+1,i) * node_sens[1]
							+ m_mapping_matrix(j*3+2,i) * node_sens[2];
    			}
    			m_filtered_dCdX[i_ID] = dCds_i;
    		}
    	}

    	// Console output for information
    	std::cout << "> Finished mapping sensitivities to design space!" << std::endl;
    	std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << std::endl;

    	KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_design_update_to_geometry_space()
    {
    	KRATOS_TRY;

    	// Measure time of mapping
    	boost::timer mapping_time;
    	std::cout << "\n> Start mapping design update to geometry space..." << std::endl;

    	// Loop over all design variables
    	for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
    	{
    		// Get node information
    		int i_ID = node_i->Id();

    		// Get the id of the node in the mapping matrix
    		int i = node_i->GetValue(MAPPING_MATRIX_ID);

    		// Instead of performing spatial search again, we read the results obtained from the computation of the mapping matrix
    		unsigned int number_of_nodes_affected = m_number_of_nodes_affected[i_ID];
    		PointVector nodes_affected(number_of_nodes_affected);
    		nodes_affected = m_listOf_nodesAffected[i_ID];

    		// Mapping of design variable update to geometry space
    		// - looping only over neighbor nodes and summ all shape update contributions form neighbor nodes
    		array_3d shape_update(3,0.0);
    		for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
    		{
    			// Get node information
    			int j_ID = nodes_affected[j_itr]->Id();
    			ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];

    			// Get the id of the node in the mapping matrix
    			int j = node_j.GetValue(MAPPING_MATRIX_ID);

    			// Store shape update contribution in global node list to update after all updates have been computed
    			shape_update[0] += m_mapping_matrix(3*i+0,j) * m_design_variable_update[j_ID];
    			shape_update[1] += m_mapping_matrix(3*i+1,j) * m_design_variable_update[j_ID];
    			shape_update[2] += m_mapping_matrix(3*i+2,j) * m_design_variable_update[j_ID];
    		}

    		// Store shape update
    		noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = shape_update;
    	}

    	// Update of coordinates and absolute values AFTER shape update is modified according to special conditions
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
    			// Project shape update at boundary on specified boundary plane (Remove component that is normal to the boundary plane)
    			if(node_i->FastGetSolutionStepValue(IS_ON_BOUNDARY))
    			{
    				array_3d boundary_plane = node_i->FastGetSolutionStepValue(BOUNDARY_PLANE);
    				array_3d original_update = node_i->FastGetSolutionStepValue(SHAPE_UPDATE);
    				array_3d projected_update = original_update - (inner_prod(original_update,boundary_plane))*boundary_plane/norm_2(boundary_plane);
    				noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = projected_update;
    			}
    		}

    		// Update coordinates
    		array_3d shape_update = node_i->FastGetSolutionStepValue(SHAPE_UPDATE);
    		node_i->X() += shape_update[0];
    		node_i->Y() += shape_update[1];
    		node_i->Z() += shape_update[2];

    		// Add final shape update to previous updates
    		noalias(node_i->FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE)) += shape_update;
    	}

    	// Console output for information
    	std::cout << "> Finished mapping design update to geometry space!" << std::endl;
    	std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << std::endl;

    	KRATOS_CATCH("");
    }

    // ==============================================================================
    // General optimization operations
    // ==============================================================================
    void compute_design_update( double step_size )
    {
        KRATOS_TRY;

        // Compute max norm of search direction
        double max_norm_search_dir = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();
            if(fabs(m_search_direction[i_ID])>max_norm_search_dir)
                max_norm_search_dir = fabs(m_search_direction[i_ID]);
        }

        // computation of update of design variable & update search direction to its maxnorm to eliminate dependency of
        // value of constant step size on lenght of search direction
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();
            m_design_variable_update[i_ID] = step_size * ( m_search_direction[i_ID] / max_norm_search_dir );
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For running unconstrained descent methods
    // ==============================================================================
    void compute_search_direction_steepest_descent()
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "\n> No constraints given or active. The negative objective gradient is chosen as search direction..." << std::endl;

        // Clear search direction
        m_search_direction.clear();

        // search direction is negative of filtered gradient
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();
            m_search_direction[i_ID] = -1.0 * m_filtered_dFdX[i_ID];
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For running augmented Lagrange method
    // ==============================================================================
    void initialize_augmented_lagrange( boost::python::dict py_constraints,
                                        double penalty_fac,
                                        double gamma,
                                        double penalty_fac_max,
                                        double lambda_0 )
    {
        KRATOS_TRY;

        // Store factors
        m_penalty_fac = penalty_fac;
        m_gamma = gamma;
        m_penalty_fac_max = penalty_fac_max;

        // Initialize a lagrange multiplier for each constraint
        boost::python::list py_C_ids = py_constraints.keys();
        for (unsigned int i = 0; i < boost::python::len(py_C_ids); ++i)
        {
            // Read constraint id from python dict and transform to cpp id
            boost::python::extract<const char*> py_C_id(py_C_ids[i]);
            unsigned int cpp_C_id = py2cpp_C_id[py_C_id];

            // Initialize multiplier using the cpp constraint id
            m_lambda[cpp_C_id] = lambda_0;
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void compute_search_direction_augmented_lagrange(  boost::python::dict py_constraints,
                                                       boost::python::dict py_response )
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "> Search direction is computed as steepest descent direction of augmented Lagrange function..." << std::endl;

        // For the computation of the eucledian norm of the search direction (for normalization)
        double norm_search_direction = 0.0;

        // Compute search direction for each node
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();

            // First term in the augmented lagrange function
            m_search_direction[i_ID] = m_filtered_dFdX[i_ID];

            // Loop over all constraints and add corresponding terms to Lagrange lagrange function
            boost::python::list py_C_ids = py_constraints.keys();
            for (unsigned int i = 0; i < boost::python::len(py_C_ids); ++i)
            {
                // Extract values to right type
                boost::python::extract<double> c(py_response[py_C_ids[i]]["func"]);
                boost::python::extract<const char*> C_type(py_constraints[py_C_ids[i]]["type"]);

                // Map python constraint id to cpp id
                boost::python::extract<const char*> py_C_id(py_C_ids[i]);
                unsigned int cpp_C_id = py2cpp_C_id[py_C_id];

                // Different operations corresponding to equality and inequality constraints
                if( std::strcmp(C_type,"eq")==0 )
                {
                    // Contributions from constraints
                    m_search_direction[i_ID] += m_lambda[cpp_C_id] * m_filtered_dCdX[i_ID] + 2 * m_penalty_fac * c * m_filtered_dCdX[i_ID];
                }
                else if( std::strcmp(C_type,"eq")==0 )
                {
                    KRATOS_THROW_ERROR(std::runtime_error, "Type of constraint not yet implemented!",C_type);
                }
                else
                    KRATOS_THROW_ERROR(std::runtime_error, "Wrong type of constraint. Choose either EQ or INEQ!",C_type);

                // Add contribution to the eucledian norm
                norm_search_direction += m_search_direction[i_ID] * m_search_direction[i_ID];
            }
        }

        // Compute eucledian norm
        norm_search_direction = sqrt(norm_search_direction);

        // Search direction is usually chosen as negative of gradient of augmented lagrange function
        // Furthermore we normalize the search direction
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
             int i_ID = node_i->Id();
            m_search_direction[i_ID] = -1.0 * m_search_direction[i_ID] / norm_search_direction;
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void udpate_augmented_lagrange_parameters( boost::python::dict py_constraints,
                                               boost::python::dict py_response )
    {
        KRATOS_TRY;

        // Loop over all constraints
        boost::python::list py_C_ids = py_constraints.keys();
        for (unsigned int i = 0; i < boost::python::len(py_C_ids); ++i)
        {
            // Extract values to right type
            boost::python::extract<double> c(py_response[py_C_ids[i]]["func"]);
            boost::python::extract<const char*> C_type(py_constraints[py_C_ids[i]]["type"]);

            // Map python constraint id to cpp id
            boost::python::extract<const char*> py_C_id(py_C_ids[i]);
            unsigned int cpp_C_id = py2cpp_C_id[py_C_id];

            // Different operations corresponding to equality and inequality constraints
            if( std::strcmp(C_type,"eq")==0 )
            {
                // Update lambda corresponding to current constraint
                m_lambda[cpp_C_id] += 2 * m_penalty_fac * c;
            }
            else if( std::strcmp(C_type,"ineq")==0 )
            {
                KRATOS_THROW_ERROR(std::runtime_error, "Type of constraint not yet implemented!",C_type);
            }
            else
                KRATOS_THROW_ERROR(std::runtime_error, "Wrong type of constraint. Choose either EQ or INEQ!",C_type);
        }

        // Update penalty factor
        m_penalty_fac = m_gamma * m_penalty_fac;
        if(m_penalty_fac>m_penalty_fac_max)
            m_penalty_fac = m_penalty_fac_max;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    double get_penalty_fac()
    {
        KRATOS_TRY;

        return m_penalty_fac;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    double get_lambda(const char* py_C_id)
    {
        KRATOS_TRY;

        // Translate python dictionary id to cpp map id
        unsigned int cpp_C_id = py2cpp_C_id[py_C_id];

        // Read value
        return m_lambda[cpp_C_id];

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    double get_value_of_augmented_lagrangian( std::string only_F_id,
                                              boost::python::dict py_constraints,
                                              boost::python::dict py_response )
    {
        KRATOS_TRY;

        // Compute value of augmented lagrange function and store on the following variable
        // First contribution is value of objective function
        boost::python::extract<double> f(py_response[only_F_id]["func"]);
        double return_value = f;

        // Loop over all constraints
        boost::python::list py_C_ids = py_constraints.keys();
        for (unsigned int i = 0; i < boost::python::len(py_C_ids); ++i)
        {
            // Extract values to right type
            boost::python::extract<double> c(py_response[py_C_ids[i]]["func"]);
            boost::python::extract<const char*> C_type(py_constraints[py_C_ids[i]]["type"]);

            // Map python constraint id to cpp id
            boost::python::extract<const char*> py_C_id(py_C_ids[i]);
            unsigned int cpp_C_id = py2cpp_C_id[py_C_id];

            // Different operations corresponding to equality and inequality constraints
            if( std::strcmp(C_type,"eq")==0 )
            {
                // Contributions from constraints
                return_value += m_lambda[cpp_C_id] * c + m_penalty_fac * c * c;
            }
            else if( std::strcmp(C_type,"ineq")==0 )
            {
                KRATOS_THROW_ERROR(std::runtime_error, "Type of constraint not yet implemented!",C_type);
            }
            else
                KRATOS_THROW_ERROR(std::runtime_error, "Wrong type of constraint. Choose either EQ or INEQ!",C_type);
        }

        return return_value;

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For running penalized projection method
    // ==============================================================================
    void compute_search_direction_penalized_projection( double c )
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "> Constraint is active. Modified search direction on the constraint hyperplane is computed..." << std::endl;

        // Compute norm of constraint and objective gradient
        double norm_2_dFdX = 0.0;
        double norm_2_dCdX = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();
            norm_2_dFdX += m_filtered_dFdX[i_ID] * m_filtered_dFdX[i_ID];
            norm_2_dCdX += m_filtered_dCdX[i_ID] * m_filtered_dCdX[i_ID];
        }
        norm_2_dFdX = sqrt(norm_2_dFdX);
        norm_2_dCdX = sqrt(norm_2_dCdX);

        // Compute dot product of objective gradient and normalized constraint gradient
        double dot_dFds_dCds = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();
            dot_dFds_dCds += m_filtered_dFdX[i_ID] * (m_filtered_dCdX[i_ID] / norm_2_dCdX);
        }

        // Compute modified search direction (negative of modified objective derivative)
        double norm_2_search_dir = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();
            m_search_direction[i_ID] = -1 * (m_filtered_dFdX[i_ID] - dot_dFds_dCds * (m_filtered_dCdX[i_ID] / norm_2_dCdX) + c * norm_2_dFdX * (m_filtered_dCdX[i_ID] / norm_2_dCdX));

            // Compute norm
            norm_2_search_dir += m_search_direction[i_ID] * m_search_direction[i_ID];
        }
        norm_2_search_dir = sqrt(norm_2_search_dir);

        // Normalize search direction
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i_ID = node_i->Id();
            m_search_direction[i_ID] /= norm_2_search_dir;
        }

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
        return "VertexMorphingUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "VertexMorphingUtilities";
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
    const int m_domain_size;
    std::map<const char*,unsigned int> py2cpp_F_id;
    std::map<const char*,unsigned int> py2cpp_C_id;
    const double m_filter_size;
    const unsigned int m_number_of_design_variables;
    const int m_max_nodes_affected;

    // ==============================================================================
    // General working arrays
    // ==============================================================================
    std::map<int,PointVector> m_listOf_nodesAffected;
    std::map<int,unsigned int> m_number_of_nodes_affected;
    std::map<int,double> m_design_variable_update;
    std::map<int,double> m_filtered_dFdX;
    std::map<int,double> m_filtered_dCdX;
    std::map<int,double> m_search_direction;

    SparseMatrixType m_mapping_matrix;


    // ==============================================================================
    // For running augmented Lagrange method
    // ==============================================================================
    double m_penalty_fac;
    double m_gamma;
    double m_penalty_fac_max;
    std::map<int,double> m_lambda;

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
//      VertexMorphingUtilities& operator=(VertexMorphingUtilities const& rOther);

    /// Copy constructor.
//      VertexMorphingUtilities(VertexMorphingUtilities const& rOther);


    ///@}

}; // Class VertexMorphingUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // VERTEX_MORPHING_UTILITIES_H
