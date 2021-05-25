
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

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_utilities/kernel_functions.h"

namespace Kratos {

KernelFunction::UniquePointer KernelFunction::New(const std::string FunctionType, const double Radius)
    {
        if(FunctionType.compare("linear"))
            return Kratos::make_unique<LinearKernelFunction>(Radius);
        else if (FunctionType.compare("cosine"))
            return Kratos::make_unique<CosineKernelFunction>(Radius);
        else if (FunctionType.compare("quartic"))
            return Kratos::make_unique<QuarticKernelFunction>(Radius);
        else if (FunctionType.compare("gaussian"))
            return Kratos::make_unique<GaussianKernelFunction>(Radius);
        else if (FunctionType.compare("constant"))
            return Kratos::make_unique<ConstantKernelFunction>(Radius);
        else
            KRATOS_ERROR<<"Unknown kernel function of type : "<<FunctionType<<std::endl;
    }
}