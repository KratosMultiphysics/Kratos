// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Aditya Ghantasala, https://github.com/adityaghantasala
//
// ==============================================================================

#ifndef FILTER_FUNCTION_H
#define FILTER_FUNCTION_H

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
/**
 * FilterFunction implementations.
*/

class FilterFunction
{
  public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double, 3> Array3DType;
    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(FilterFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FilterFunction(const std::string FilterFunctionType, const double Radius)
        : mRadius(Radius)
    {
        // Set type of weighting function

        // Type 1: Gaussian function
        if (FilterFunctionType == "gaussian")
            mFilterFunctional =  [](double radius, double distance) {return std::max(0.0, exp(-(distance*distance) / (2 * radius * radius / 9.0)));};

        // Type 2: Linear function
        else if (FilterFunctionType == "linear")
            mFilterFunctional =  [](double radius, double distance) {return std::max(0.0, (radius - distance) / radius);};

        // Type 3: Constant function
        else if (FilterFunctionType == "constant")
            mFilterFunctional = [](double radius, double distance) {return 1.0;};

        // Type 4: Cosine function
        else if (FilterFunctionType == "cosine")
            mFilterFunctional = [](double radius, double distance) {return std::max(0.0, 1-0.5*(1-cos(PI/radius*distance)));};

        // Type 5: Quartic function
        else if (FilterFunctionType == "quartic")
            mFilterFunctional = [](double radius, double distance) {return std::max(0.0, (pow(distance-radius,4.0)/pow(radius,4.0)));};

        // Throw error message in case of wrong specification
        else
            KRATOS_ERROR << "Specified kernel function of type : "<< FilterFunctionType << " is not recognized. \n \t Options are: constant, linear , gaussian, cosine, quartic." << std::endl;
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

    double ComputeWeight(const Array3DType& ICoord, const Array3DType& JCoord) const
    {
        KRATOS_TRY;

        // Compute distance vector
        const double distance = GetDistance(ICoord, JCoord);

        // Depending on which weighting function is chosen, compute weight
        return mFilterFunctional(mRadius, distance);

        KRATOS_CATCH("");
    }

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
        rOStream << "FilterFunction ";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mRadius;
    std::function<double (double, double)> mFilterFunctional;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    double inline GetDistance(const Array3DType ICoord, const Array3DType JCoord) const
    {
        const Array3DType dist_vector = ICoord - JCoord;
        return sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

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
