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
                           boost::python::dict py_constraints )
        : mr_opt_model_part(model_part)
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
    void compute_design_update( double step_size )
    {
        KRATOS_TRY;

        // Compute max norm of search direction
        double max_norm_search_dir = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            if(fabs(node_i->FastGetSolutionStepValue(SEARCH_DIRECTION))>max_norm_search_dir)
                max_norm_search_dir = fabs(node_i->FastGetSolutionStepValue(SEARCH_DIRECTION));
        }

        // Computation of update of design variable. Before search direction is scaled to its maxnorm to improve dependency of
        // step size on lenght of search direction
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            node_i->FastGetSolutionStepValue(DESIGN_UPDATE) = step_size * ( node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) / max_norm_search_dir );
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

        // Working variable
        double search_direction_i;

        // For the computation of the eucledian norm of the search direction (for normalization)
        double norm_search_direction = 0.0;

        // Compute search direction for each node
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            // Get information about sensitivities in design space
            double dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);

            // First term in the augmented lagrange function
            search_direction_i = node_i->FastGetSolutionStepValue(DESIGN_UPDATE);

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
                    search_direction_i += m_lambda[cpp_C_id] * dCds_i + 2 * m_penalty_fac * c * dCds_i;
                }
                else if( std::strcmp(C_type,"eq")==0 )
                {
                    KRATOS_THROW_ERROR(std::runtime_error, "Type of constraint not yet implemented!",C_type);
                }
                else
                    KRATOS_THROW_ERROR(std::runtime_error, "Wrong type of constraint. Choose either EQ or INEQ!",C_type);

                // Add contribution to the eucledian norm
                norm_search_direction += search_direction_i * search_direction_i;
            }
        }

        // Compute eucledian norm
        norm_search_direction = sqrt(norm_search_direction);

        // Search direction is usually chosen as negative of gradient of augmented lagrange function
        // Furthermore we normalize the search direction
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) = -1.0 * search_direction_i / norm_search_direction;
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
        std::cout << "\n> Constraint is active. Modified search direction on the constraint hyperplane is computed..." << std::endl;

        // Compute norm of constraint and objective gradient
        double norm_2_dFds_i = 0.0;
        double norm_2_dCds_i = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            double dFds_i = node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
            double dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            norm_2_dFds_i += dFds_i * dFds_i;
            norm_2_dCds_i += dCds_i * dCds_i;
        }
        norm_2_dFds_i = sqrt(norm_2_dFds_i);
        norm_2_dCds_i = sqrt(norm_2_dCds_i);

        // Compute dot product of objective gradient and normalized constraint gradient
        double dot_dFds_dCds = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            double dFds_i = node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
            double dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            dot_dFds_dCds += dFds_i * (dCds_i / norm_2_dCds_i);
        }

        // Compute modified search direction (negative of modified objective derivative)
        double norm_2_search_dir = 0.0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            double dFds_i = node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
            double dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            double search_dir_i = -1 * (dFds_i - dot_dFds_dCds * (dCds_i / norm_2_dCds_i) + c * norm_2_dFds_i * (dCds_i / norm_2_dCds_i));

            // Assign search direction
            node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) = search_dir_i;

            // Compute norm
            norm_2_search_dir += search_dir_i * search_dir_i;
        }
        norm_2_search_dir = sqrt(norm_2_search_dir);

        // Normalize search direction
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) = node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) / norm_2_search_dir;
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
