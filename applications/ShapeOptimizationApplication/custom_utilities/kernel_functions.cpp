
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

KernelFunction::Pointer KernelFunction::New(const std::string FunctionType, const double Radius)
    {
        if(FunctionType.compare("linear"))
            return Kratos::make_shared<LinearKernelFunction>(Radius);
        else if (FunctionType.compare("cosine"))
            return Kratos::make_shared<CosineKernelFunction>(Radius);
        else if (FunctionType.compare("quartic"))
            return Kratos::make_shared<QuarticKernelFunction>(Radius);
        else if (FunctionType.compare("gaussian"))
            return Kratos::make_shared<GaussianKernelFunction>(Radius);
        else if (FunctionType.compare("constant"))
            return Kratos::make_shared<ConstantKernelFunction>(Radius);
        else
            KRATOS_ERROR<<"Unknown kernel function of type : "<<FunctionType<<std::endl;
    }
}