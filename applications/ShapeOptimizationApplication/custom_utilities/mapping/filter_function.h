// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef FILTER_FUNCTION_H
#define FILTER_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
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

class FilterFunction
{
  public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double, 3> array_3d;

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(FilterFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FilterFunction(std::string filter_function_type, double filter_size)
        : m_filter_size(filter_size)
    {
        // Create strings to compare to
        std::string gaussian("gaussian");
        std::string linear("linear");
        std::string constant("constant");

        // Set type of weighting function

        // Type 1: Gaussian function
        if (filter_function_type.compare(gaussian) == 0)
            m_filter_function_type = 1;

        // Type 2: Linear function
        else if (filter_function_type.compare(linear) == 0)
            m_filter_function_type = 2;

        // Type 3: Constant function
        else if (filter_function_type.compare(constant) == 0)
            m_filter_function_type = 3;

        // Throw error message in case of wrong specification
        else
            KRATOS_ERROR << "Specified filter function type not recognized. Options are: constant, linear , gaussian. Specified: " << std::endl;
    }

    /// Destructor.
    virtual ~FilterFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================

    double compute_weight(array_3d i_coord, array_3d j_coord)
    {
        KRATOS_TRY;

        // Compute distance vector
        array_3d dist_vector = i_coord - j_coord;

        // Depending on which weighting function is chosen, compute weight
        double weight_ij = 0.0;
        switch (m_filter_function_type)
        {
        // Gaussian filter
        case 1:
        {
            // Compute squared distance
            double squared_scalar_distance = dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2];
            // Compute weight
            // Note that we do not compute the square root of the distances to save this expensive computation (it is not needed here)
            weight_ij = std::max(0.0, exp(-squared_scalar_distance / (2 * m_filter_size * m_filter_size / 9.0)));
            break;
        }
        // Linear filter
        case 2:
        {
            // Compute distance
            double distance = sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
            // Compute weight
            weight_ij = std::max(0.0, (m_filter_size - distance) / m_filter_size);
            break;
        }
        // Constant filter
        case 3:
        {
            // Compute weight
            weight_ij = 1.0;
            break;
        }
        }

        return weight_ij;

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
        return "FilterFunction";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "FilterFunction";
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

    double m_filter_size;
    unsigned int m_filter_function_type;

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
    //      FilterFunction& operator=(FilterFunction const& rOther);

    /// Copy constructor.
    //      FilterFunction(FilterFunction const& rOther);

    ///@}

}; // Class FilterFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // FILTER_FUNCTION_H
