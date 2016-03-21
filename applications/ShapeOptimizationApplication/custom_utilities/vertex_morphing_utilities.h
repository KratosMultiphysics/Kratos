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
// Note that in the following "../kratos/" was inserted to allow for proper
// code-highlighting and function following in qtCreator
#include "../kratos/includes/define.h"
#include "../kratos/includes/define.h"
#include "../kratos/processes/process.h"
#include "../kratos/includes/node.h"
#include "../kratos/includes/element.h"
#include "../kratos/includes/model_part.h"
#include "../kratos/includes/kratos_flags.h"
#include "../kratos/spatial_containers/spatial_containers.h"
#include "../kratos/utilities/timer.h"
#include "../kratos/processes/node_erase_process.h"
#include "../kratos/utilities/binbased_fast_point_locator.h"
#include "../kratos/utilities/normal_calculation_utils.h"
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

        // Normalize area normal and store in respective variable
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];

            // Normalize normal and assign to solution step value
            array_3d area_normal = node_i.FastGetSolutionStepValue(NORMAL);
            array_3d normalized_normal = area_normal / norm_2(area_normal);
            noalias(node_i.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)) = normalized_normal;

        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void project_grad_on_unit_surface_normal( bool constraint_given )
    {
        KRATOS_TRY;

        // Initialize variables for normalization
        double J_norm = 0.0;
        double C_norm = 0.0;

        // We loop over all nodes and compute the part of the sensitivity which is in direction
        // to the surface normal
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];

            // We compute dFdX_n = (dFdX \cdot n) * n
            array_3d node_sens = node_i.FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            array_3d node_normal = node_i.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
            array_3d normal_node_sens =  inner_prod(node_sens,node_normal) * node_normal;

            // Assign resulting sensitivity back to node
            noalias(node_i.FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY)) = normal_node_sens;

            // Compute contribution to objective norm denominator
            J_norm += inner_prod(normal_node_sens,normal_node_sens);

            // Repeat for constraint
            if(constraint_given)
            {
                // We compute dFdX_n = (dFdX \cdot n) * n
                node_sens = node_i.FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                node_normal = node_i.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
                normal_node_sens =  inner_prod(node_sens,node_normal) * node_normal;

                // Assign resulting sensitivity back to node
                noalias(node_i.FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY)) = normal_node_sens;

                // Compute contribution to objective norm denominator
                C_norm += inner_prod(normal_node_sens,normal_node_sens);
            }
        }

        // Compute contribution to objective norm denominator
        J_norm = sqrt(J_norm);
        C_norm = sqrt(C_norm);

        // Normalize sensitivities
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];

            array_3d node_sens = node_i.FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            noalias(node_i.FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY)) = node_sens/J_norm;

          // Repeat for constraint
            if(constraint_given)
            {
                node_sens = node_i.FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                noalias(node_i.FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY)) = node_sens/C_norm;
            }
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For perfoming Vertex Morphing
    // ==============================================================================
    void filter_gradients( bool constraint_given )
    {
        KRATOS_TRY;

        // Measure time of filtering
        boost::timer filtering_time;

        // Some output for information
        std::cout << "> Start backward filtering..." << std::endl;

        // Initialization of needed variables
        std::map<int,double> dFds; // Objective gradients
        std::map<int,double> dCds; // Constraint gradients


        // Creating an auxiliary list for the nodes to be searched on
        PointVector list_of_nodes;

        // Start constructing and computing the kdtree
        typedef Bucket< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;

        // starting calculating time of construction of the kdtree
        boost::timer kdtree_construction;
        for (ModelPart::NodesContainerType::iterator node_it = mr_opt_model_part.NodesBegin(); node_it != mr_opt_model_part.NodesEnd(); ++node_it)
        {
            PointTypePointer pnode = *(node_it.base());

            // Putting the nodes of interest in an auxiliary list
            list_of_nodes.push_back(pnode);
        }
        std::cout << "> Construction time of KDTree: " << kdtree_construction.elapsed() << std::endl;

        // Arrays needed for spatial search
        PointVector nodes_affected(m_max_nodes_affected);
        DistanceVector resulting_squared_distances(m_max_nodes_affected);

        // compute the tree with the position of the nodes
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), m_max_nodes_affected);

        // Loop over all design variables
        double max_dFds = 0.0;
        double max_dCds = 0.0;
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];

            // perform spatial search for current node
            unsigned int number_of_nodes_affected;
            number_of_nodes_affected = nodes_tree.SearchInRadius(node_i, m_filter_size, nodes_affected.begin(),resulting_squared_distances.begin(), m_max_nodes_affected);

            // Store results to reuse later in the forward mapping
            m_listOf_nodesAffected[i_ID] = nodes_affected;
            m_number_of_nodes_affected[i_ID] = number_of_nodes_affected;

            // Compute weights
            std::map<int,array_3d> Ai;
            Ai = compute_weights(node_i,nodes_affected,number_of_nodes_affected);

            // Compute filtered gradients for objective function (Do backward mapping, filtering)
            double dFds_i = 0.0;
            for(unsigned int j = 0 ; j<number_of_nodes_affected ; j++)
            {
                int j_ID = nodes_affected[j]->Id();
                ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];

                // Ask for the sensitiviteis
                array_3d node_sens = node_j.FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);

                // Compute dot-product
                dFds_i += Ai[j_ID][0] * node_sens[0] + Ai[j_ID][1] * node_sens[1] + Ai[j_ID][2] * node_sens[2];
            }
            dFds[i_ID] = dFds_i;

            // Store max value for later normalization
            if(fabs(dFds_i)>max_dFds)
                max_dFds = fabs(dFds_i);

            // If constrainted optimization, then also compute filtered gradients for the constraints (dCds)
            if(constraint_given)
            {
                // Do backward mapping (filtering)
                double dCds_i = 0.0;
                for(unsigned int j = 0 ; j<number_of_nodes_affected ; j++)
                {
                    int j_ID = nodes_affected[j]->Id();
                    ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];

                    // Ask for the sensitiviteis
                    array_3d node_sens = node_j.FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);

                    // Compute dot-product
                    dCds_i += Ai[j_ID][0] * node_sens[0] + Ai[j_ID][1] * node_sens[1] + Ai[j_ID][2] * node_sens[2];
                }
                dCds[i_ID] = dCds_i;

                // Store max value for later normalization
                if(fabs(dCds_i)>max_dCds)
                    max_dCds = fabs(dCds_i);
            }
        }

        // Filtered gradients are computed as dFds (dCds) normalized by its max-norm
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            m_filtered_dFdX[i_ID] = dFds[i_ID] / max_dFds;

            // If constraints are active, also normalize the filtered constraint gradients
            if(constraint_given)
                m_filtered_dCdX[i_ID] = dCds[i_ID] / max_dCds;
        }

        std::cout << "> Finished backward filtering!" << std::endl;
        std::cout << "> Time needed for backward filtering: " << filtering_time.elapsed() << std::endl;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    std::map<int,array_3d> compute_weights( ModelPart::NodeType& given_node,
                                            PointVector nodes_affected,
                                            unsigned int number_of_nodes_affected )
    {
        KRATOS_TRY;

        // Initializize variables
        std::map<int,array_3d> Ai;
        double summ_of_all_weights = 0.0;
        array_3d i_coord = given_node.Coordinates();

        // compute weights
        array_3d dist_vector(3,0.0);
        array_3d j_coord(3,0.0);;
        for(unsigned int j = 0 ; j<number_of_nodes_affected ; j++)
        {
            unsigned int j_ID = nodes_affected[j]->Id();
            j_coord[0] = nodes_affected[j]->X();
            j_coord[1] = nodes_affected[j]->Y();
            j_coord[2] = nodes_affected[j]->Z();

            dist_vector = i_coord - j_coord;
            double squared_scalar_distance = dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2];

            // Computation of weight according specified weighting function (here it is the implementation of Carat)
            // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
            double Aij = exp(-squared_scalar_distance/(2*m_filter_size*m_filter_size/9.0));
            Ai[j_ID][0] = Aij;
            Ai[j_ID][1] = Aij;
            Ai[j_ID][2] = Aij;

            // Computed for necessary integration of weighting function (post-scaling)
            summ_of_all_weights += Aij;
        }

        // Here we perform the normalization by the sum of all weights & multiply by the node normal (nodal director)
        // In this way we implicitly preserve the in-plane mesh quality implicitely (pure Heuristic of Majid Hojjat)
        for(unsigned int j = 0 ; j<number_of_nodes_affected ; j++)
        {
            unsigned int j_ID = nodes_affected[j]->Id();
            Ai[j_ID] /= summ_of_all_weights;

            Ai[j_ID][0] *= given_node.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)[0];
            Ai[j_ID][1] *= given_node.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)[1];
            Ai[j_ID][2] *= given_node.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)[2];
        }
        return Ai;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void update_design_variable( double step_size )
    {
        KRATOS_TRY;

        // computation of update of design variable
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            m_design_variable_update[i_ID] = step_size * m_search_direction[i_ID];
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // General optimization operations
    // ==============================================================================
    void update_shape()
    {
        KRATOS_TRY;

        // Initialize variables
        std::map<int,array_3d> old_surface_nodes;
        std::map<int,array_3d> shape_update;

        // Store old surface nodes to reproduce later. Furthermore initialize the shape update
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];

            // Store old coordinates
            array_3d old_node(3,0.0);
            old_node[0] = node_i.X();
            old_node[1] = node_i.Y();
            old_node[2] = node_i.Z();
            old_surface_nodes[i_ID] = old_node;

            // Initialize array
            array_3d zero_array(3,0.0);
            shape_update[i_ID] = zero_array;
        }
        // Perform the forward mapping

        // Loop over all design variables
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];

            // Instead of performing spatial search again, we read the results obtained from "computeSearchDirection()"
            unsigned int number_of_nodes_affected = m_number_of_nodes_affected[i_ID];
            PointVector nodes_affected(number_of_nodes_affected);
            nodes_affected = m_listOf_nodesAffected[i_ID];

            // Compute weights
            std::map<int,array_3d> Ai;
            Ai = compute_weights(node_i,nodes_affected,number_of_nodes_affected);

            // Do forward mapping (filtering)
            for(unsigned int j = 0 ; j<number_of_nodes_affected ; j++)
            {
                unsigned int affectedNode_j_ID = nodes_affected[j]->Id();
                array_3d shape_update_affectedNode_j_ID(3,0.0);

                shape_update_affectedNode_j_ID[0] = Ai[affectedNode_j_ID][0] * m_design_variable_update[i_ID];
                shape_update_affectedNode_j_ID[1] = Ai[affectedNode_j_ID][1] * m_design_variable_update[i_ID];
                shape_update_affectedNode_j_ID[2] = Ai[affectedNode_j_ID][2] * m_design_variable_update[i_ID];

//                // Update coordinates in the mesh (note that this is only in accordance to carat, it is not yet known why
//                // this is done since by that the iteration sequence through the design variables plays an role)
//                mr_opt_model_part.Nodes()[affectedNode_j_ID].X() += shape_update_affectedNode_j_ID[0];
//                mr_opt_model_part.Nodes()[affectedNode_j_ID].Y() += shape_update_affectedNode_j_ID[1];
//                mr_opt_model_part.Nodes()[affectedNode_j_ID].Z() += shape_update_affectedNode_j_ID[2];

                // Store shape update contribution in global node list to update after all updates have been computed
                // Note that this is different as in carat. It is assumed that carat does a mistake here
                shape_update[affectedNode_j_ID][0] += shape_update_affectedNode_j_ID[0];
                shape_update[affectedNode_j_ID][1] += shape_update_affectedNode_j_ID[1];
                shape_update[affectedNode_j_ID][2] += shape_update_affectedNode_j_ID[2];
            }
        }

        // Store shape updates and update coordinates in the mesh AFTER all shape updates have been computed (as mentioned above, this is different as in carat)
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];

            // Update coordinates
            node_i.X() += shape_update[i_ID][0];
            node_i.Y() += shape_update[i_ID][1];
            node_i.Z() += shape_update[i_ID][2];

            // Save shape update
            noalias(node_i.FastGetSolutionStepValue(SHAPE_UPDATE)) = shape_update[i_ID];
            noalias(node_i.FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE)) += shape_update[i_ID];
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
        std::cout << "> No constraints given or active. The negative objective gradient is chosen as search direction..." << std::endl;

        // Clear search direction
        m_search_direction.clear();

        // search direction is negative of filtered gradient
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
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
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();

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
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
             int i_ID = node_itr->Id();
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
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            norm_2_dFdX += m_filtered_dFdX[i_ID] * m_filtered_dFdX[i_ID];
            norm_2_dCdX += m_filtered_dCdX[i_ID] * m_filtered_dCdX[i_ID];
        }
        norm_2_dFdX = sqrt(norm_2_dFdX);
        norm_2_dCdX = sqrt(norm_2_dCdX);

        // Compute dot product of objective gradient and normalized constraint gradient
        double dot_dFds_dCds = 0.0;
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            dot_dFds_dCds += m_filtered_dFdX[i_ID] * (m_filtered_dCdX[i_ID] / norm_2_dCdX);
        }

        // Compute modified search direction (negative of modified objective derivative)
        double norm_2_search_dir = 0.0;
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
            m_search_direction[i_ID] = -1 * (m_filtered_dFdX[i_ID] - dot_dFds_dCds * (m_filtered_dCdX[i_ID] / norm_2_dCdX) + c * norm_2_dFdX * (m_filtered_dCdX[i_ID] / norm_2_dCdX));

            // Compute norm
            norm_2_search_dir += m_search_direction[i_ID] * m_search_direction[i_ID];
        }
        norm_2_search_dir = sqrt(norm_2_search_dir);

        // Normalize search direction
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            int i_ID = node_itr->Id();
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
