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

#ifndef DAMPING_FUNCTION_H
#define DAMPING_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip> // for std::setprecision

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

class DampingFunction
{
  public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double, 3> array_3d;

    /// Pointer definition of DampingFunction
    KRATOS_CLASS_POINTER_DEFINITION(DampingFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DampingFunction(std::string damping_function_type, double damping_radius)
        : m_damping_radius(damping_radius)
    {
        // Set precision for output
        std::cout.precision(12);

        // Create strings to compare to
        std::string cosine("cosine");
        std::string linear("linear");

        // Set type of dymping function

        // Type 1: Cosine function
        if (damping_function_type.compare(cosine) == 0)
            m_damping_function_type = 1;

        // Type 2: Linear function
        else if (damping_function_type.compare(linear) == 0)
            m_damping_function_type = 2;

        // Throw error message in case of wrong specification
        else
            KRATOS_THROW_ERROR(std::invalid_argument, "Specified damping function type not recognized. Options are: linear , cosine. Specified: ",damping_function_type);
    }

    /// Destructor.
    virtual ~DampingFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================

    double compute_damping_factor(array_3d i_coord, array_3d j_coord)
    {
        KRATOS_TRY;

        // Compute distance vector
        array_3d dist_vector = i_coord - j_coord;

        // Depending on which damping function is chosen, compute damping factor
        double damping_factor = 0.0;
        switch (m_damping_function_type)
        {
        // Cosine damping function
        case 1:
        {
            // Compute distance
            double scalar_distance = sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
            // Compute damping factor
            damping_factor = std::min(1.0, 0.5*(1-cos(PI/m_damping_radius*scalar_distance)));
            break;
        }
        // Linear damping function
        case 2:
        {
            // Compute distance
            double distance = sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
            // Compute damping factor
            damping_factor = std::min(1.0, distance / m_damping_radius);
            break;
        }
        }

        return damping_factor;

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
        return "DampingFunction";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "DampingFunction";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
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

    double m_damping_radius;
    unsigned int m_damping_function_type;

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
    //      DampingFunction& operator=(DampingFunction const& rOther);

    /// Copy constructor.
    //      DampingFunction(DampingFunction const& rOther);

    ///@}

}; // Class DampingFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // DAMPING_FUNCTION_H
