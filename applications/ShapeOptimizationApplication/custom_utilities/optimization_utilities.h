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

#ifndef OPTIMIZATION_UTILITIES_H
#define OPTIMIZATION_UTILITIES_H

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
#include "../../kratos/utilities/timer.h"
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

class OptimizationUtilities
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


    /// Pointer definition of OptimizationUtilities
    KRATOS_CLASS_POINTER_DEFINITION(OptimizationUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OptimizationUtilities( ModelPart& model_part,
                           boost::python::dict py_objectives,
                           boost::python::dict py_constraints,
                           double step_size,
                           bool normalize_search_direction )
        : mr_opt_model_part(model_part),
          m_step_size(step_size),
          m_normalize_search_direction(normalize_search_direction)
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

        // Initialize variables
        m_correction_scaling = 0.1;
    }

    /// Destructor.
    virtual ~OptimizationUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    // General optimization operations
    // ==============================================================================
    void compute_design_update()
    {
        KRATOS_TRY;

        // Computation of update of design variable. Normalization is applied if specified.
        if(m_normalize_search_direction)
        {
            // Compute max norm of search direction
            double max_norm_search_dir = 0.0;
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                array_3d search_dir = node_i->FastGetSolutionStepValue(SEARCH_DIRECTION);
                double squared_length = inner_prod(search_dir,search_dir);
                
                if(squared_length>max_norm_search_dir)
                    max_norm_search_dir = squared_length;
            }
            max_norm_search_dir = sqrt(max_norm_search_dir);

            // Compute update
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                array_3d design_update = m_step_size * ( node_i->FastGetSolutionStepValue(SEARCH_DIRECTION)/max_norm_search_dir );
                noalias(node_i->FastGetSolutionStepValue(DESIGN_UPDATE)) = design_update;

                // Sum design updates to obtain control point position
                noalias(node_i->FastGetSolutionStepValue(DESIGN_CHANGE_ABSOLUTE)) += design_update;
            }
        }
        else
        {
            // Compute update
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                array_3d design_update = m_step_size * ( node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) );
                noalias(node_i->FastGetSolutionStepValue(DESIGN_UPDATE)) = design_update;

                // Sum design updates to obtain control point position
                noalias(node_i->FastGetSolutionStepValue(DESIGN_CHANGE_ABSOLUTE)) += design_update;
            }
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

        // search direction is negative of filtered gradient
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) = -1.0 * node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For running penalized projection method
    // ==============================================================================
    void compute_projected_search_direction( double c )
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "\n> Constraint is active. Modified search direction on the constraint hyperplane is computed..." << std::endl;

        // Compute norm of constraint gradient
        double norm_2_dCds_i = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
        	array_3d dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            norm_2_dCds_i += inner_prod(dCds_i,dCds_i);
        }
       norm_2_dCds_i = sqrt(norm_2_dCds_i);

        // Compute dot product of objective gradient and normalized constraint gradient
        double dot_dFds_dCds = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
        	array_3d dFds_i = node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
        	array_3d dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            dot_dFds_dCds += inner_prod(dFds_i,(dCds_i / norm_2_dCds_i));
        }

        // Compute and assign projected search direction
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
        	array_3d dFds_i = node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
        	array_3d dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);

        	array_3d projection_term = dot_dFds_dCds * (dCds_i / norm_2_dCds_i);

            node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) = -1 * (dFds_i - projection_term);
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void correct_projected_search_direction( double c, double previous_c, boost::python::list& py_correction_scaling )
    {
        // Check correction necessary
        if(c==0)
         return;

        // Calucation of vector norms 
    	double norm_correction_term = 0.0;
    	double norm_search_direction = 0.0;
    	for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
    	{
    		array_3d correction_term = c * node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
    		norm_correction_term += inner_prod(correction_term,correction_term);

    		array_3d ds = node_i->FastGetSolutionStepValue(SEARCH_DIRECTION);
    		norm_search_direction += inner_prod(ds,ds);
    	}
    	norm_correction_term = sqrt(norm_correction_term);
    	norm_search_direction = sqrt(norm_search_direction);

    	// Three cases need to be covered

		// 1) In case we have two subsequently decreasing constraint values --> correction is fine --> leave current correction scaling

    	// 2) In case the correction jumps over the constraint (change of sign) --> correction was too big --> reduce
    	if(c*previous_c<0)
    	{
    		m_correction_scaling *= 0.5;
    		std::cout << "Correction scaling needs to decrease...." << std::endl;
    	}

    	// 3) In case we have subsequently increasing constraint value --> correction was too low --> increase
    	if(std::abs(c)>std::abs(previous_c) and c*previous_c>0)
    	{
    		std::cout << "Correction scaling needs to increase...." << std::endl;
    		m_correction_scaling = std::min(m_correction_scaling*2,1.0);
    	}
    	double correction_factor = m_correction_scaling * norm_search_direction / norm_correction_term;

        // Rewrite value in container to keep track of the value in python
    	py_correction_scaling[0] = m_correction_scaling;

        // Perform correction
    	for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
    	{
    		array_3d correction_term = correction_factor * c * node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            // KRATOS_WATCH(correction_term);
            // KRATOS_WATCH(node_i->FastGetSolutionStepValue(SEARCH_DIRECTION));
    		node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) -= correction_term;
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
        return "OptimizationUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "OptimizationUtilities";
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
    std::map<const char*,unsigned int> py2cpp_F_id;
    std::map<const char*,unsigned int> py2cpp_C_id;
    double m_step_size;
    bool m_normalize_search_direction;

    // ==============================================================================
    // For running penalized projection method
    // ==============================================================================
    double m_correction_scaling;

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
//      OptimizationUtilities& operator=(OptimizationUtilities const& rOther);

    /// Copy constructor.
//      OptimizationUtilities(OptimizationUtilities const& rOther);


    ///@}

}; // Class OptimizationUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // OPTIMIZATION_UTILITIES_H
