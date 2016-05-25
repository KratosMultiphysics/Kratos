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

#ifndef GEOMETRY_UTILITIES_H
#define GEOMETRY_UTILITIES_H

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
#include "../../kratos/utilities/normal_calculation_utils.h"
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

class GeometryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double,3> array_3d;

    /// Pointer definition of GeometryUtilities
    KRATOS_CLASS_POINTER_DEFINITION(GeometryUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometryUtilities( ModelPart& model_part,
                       const int domain_size )
        : mr_opt_model_part(model_part),
          m_domain_size(domain_size)
    {
        // Set precision for output
        std::cout.precision(12);
    }

    /// Destructor.
    virtual ~GeometryUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

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
        return "GeometryUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometryUtilities";
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
//      GeometryUtilities& operator=(GeometryUtilities const& rOther);

    /// Copy constructor.
//      GeometryUtilities(GeometryUtilities const& rOther);


    ///@}

}; // Class GeometryUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // GEOMETRY_UTILITIES_H
