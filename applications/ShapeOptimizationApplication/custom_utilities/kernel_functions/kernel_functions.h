// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala https://github.com/adityaghantasala
//
//   Based on the previous implementations of filter functions and damping functions.
//   This just unifies them so they can be used everywhere in the code from the same source.
//
// ==============================================================================

#ifndef KERNEL_FUNCTIONS_H
#define KERNEL_FUNCTIONS_H

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

/// Short class definition.
/**
 * Base class for all the kernel functions
*/
class KernelFunction
{
  public:
    ///@name Type Definitions
    ///@{
    typedef array_1d<double, 3> array_3d;

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(KernelFunction);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    KernelFunction(const double FilterSize)
        : mFilterSize(FilterSize)
    {}
    /// Destructor.
    virtual ~KernelFunction()
    {}
    ///@}
    ///@name Operations
    ///@{
        virtual double ComputeWeight(const array_3d iCoord, const array_3d jCoord) = 0;
    ///@}

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KernelFunction";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "KernelFunction";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }
    ///@}
  private:
    ///@name Member Variables
    ///@{
    double mFilterSize;
    ///@}
}; // Class KernelFunction


} // namespace Kratos.

#endif // KERNEL_FUNCTIONS_H
