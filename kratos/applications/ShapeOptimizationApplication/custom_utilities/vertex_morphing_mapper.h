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

class VertexMorphingMapper
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
                             const int max_nodes_affected )
        : mr_opt_model_part(model_part),
          m_filter_function_type(filter_function_type),
          m_filter_size(filter_size),
          m_max_nodes_affected(max_nodes_affected),
          m_number_of_design_variables(model_part.Nodes().size())
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
    void compute_mapping_matrix()
    {
        KRATOS_TRY;

        // Measure time of mapping
        boost::timer mapping_time;
        std::cout << "> Start computing mapping matrix..." << std::endl;

        // Initialize mapping matrix (zero possible entries)
        m_mapping_matrix.clear();

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

        // Prepare Weighting function to be used in the mapping
        FilterFunction filter_function( m_filter_function_type, m_filter_size );

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

            // Compute and assign weights in the mapping matrix
            double sum_weights = 0.0;
            for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
            {
            	// Get node information
                int j_ID = nodes_affected[j_itr]->Id();
                ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];
                array_3d j_coord = node_j.Coordinates();

                // Get the id of the node in the mapping matrix
                int j = node_j.GetValue(MAPPING_MATRIX_ID);

                // Computation of weight according specified weighting function
                // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                double Aij = filter_function.compute_weight(i_coord,j_coord);

                // Multiplication by the node normal (nodal director)
                array_3d n = node_j.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

                // In this way we implicitly preserve the in-plane mesh quality (pure Heuristic)
                m_mapping_matrix(3*i+0,j) = Aij*n[0];
                m_mapping_matrix(3*i+1,j) = Aij*n[1];
                m_mapping_matrix(3*i+2,j) = Aij*n[2];

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
        std::cout << "> Time needed for computation of mapping matrix: " << mapping_time.elapsed() << " s" << std::endl;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_sensitivities_to_design_space( bool constraint_given )
    {
        KRATOS_TRY;

    	// Measure time of mapping
    	boost::timer mapping_time;
    	std::cout << "\n> Start mapping sensitivities to design space..." << std::endl;

        // Map objective sensitivities

        // Assign nodal sensitivities to vector used vor mapping
        VectorType dJdX;
        dJdX.resize(m_number_of_design_variables*3);
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            array_3d node_sens = node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            dJdX[i*3+0] = node_sens[0];
            dJdX[i*3+1] = node_sens[1];
            dJdX[i*3+2] = node_sens[2];
        }

        // Perform mapping
        VectorType dJds;
        dJds.resize(m_number_of_design_variables);
        SparseSpaceType::TransposeMult(m_mapping_matrix,dJdX,dJds);

        // Assign results to nodal variables
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
               int i = node_i->GetValue(MAPPING_MATRIX_ID);
               node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY) = dJds[i];
        }

        // Repeat mapping to map constraint sensitivities
        if(constraint_given)
        {
            // Assign nodal sensitivities to vector used vor mapping
            VectorType dCdX;
            dCdX.resize(m_number_of_design_variables*3);
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d node_sens = node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                dCdX[i*3+0] = node_sens[0];
                dCdX[i*3+1] = node_sens[1];
                dCdX[i*3+2] = node_sens[2];
            }

            // Perform mapping
            VectorType dJds;
            dJds.resize(m_number_of_design_variables);
            SparseSpaceType::TransposeMult(m_mapping_matrix,dCdX,dJds);

            // Assign results to nodal variables
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                   int i = node_i->GetValue(MAPPING_MATRIX_ID);
                   node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY) = dJds[i];
            }
        }

    	// Console output for information
    	std::cout << "> Finished mapping sensitivities to design space!" << std::endl;
        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;

    	KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_design_update_to_geometry_space()
    {
    	KRATOS_TRY;

    	// Measure time of mapping
    	boost::timer mapping_time;
    	std::cout << "\n> Start mapping design update to geometry space..." << std::endl;

        // Assign design update to vector which shall be mapped
        VectorType ds;
        ds.resize(m_number_of_design_variables);
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            ds[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE);
        }

        // Perform mapping to compute shape update
        VectorType dx;
        dx.resize(m_number_of_design_variables*3);
        SparseSpaceType::Mult(m_mapping_matrix,ds,dx);

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
                // Read shape update from solution vector
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d shape_update;
                shape_update[0] = dx[3*i+0];
                shape_update[1] = dx[3*i+1];
                shape_update[2] = dx[3*i+2];

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
    const int m_max_nodes_affected;
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
