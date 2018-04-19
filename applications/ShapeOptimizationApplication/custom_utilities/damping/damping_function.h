// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef DAMPING_FUNCTION_H
#define DAMPING_FUNCTION_H

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

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
        // Create strings to compare to
        std::string cosine("cosine");
        std::string linear("linear");
        std::string quartic("quartic");

        // Set type of dymping function

        // Type 1: Cosine function
        if (damping_function_type.compare(cosine) == 0)
            m_damping_function_type = 1;

        // Type 2: Linear function
        else if (damping_function_type.compare(linear) == 0)
            m_damping_function_type = 2;

        // Type 3: Quartic function
        else if (damping_function_type.compare(quartic) == 0)
            m_damping_function_type = 3;

        // Throw error message in case of wrong specification
        else
            KRATOS_ERROR << "Specified damping function type not recognized. Options are: linear , cosine. Specified: " << damping_function_type << std::endl;
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
            double distance = sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
            // Compute damping factor
            damping_factor = std::min(1.0, 0.5*(1-cos(PI/m_damping_radius*distance)));
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
        // Quartic damping function
        case 3:
        {
            // Compute distance
            double distance = sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
            double numerator = distance-m_damping_radius;
            // Compute damping factor
            damping_factor = std::min(1.0, (1-pow(numerator,4.0)/pow(m_damping_radius,4.0)));
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
