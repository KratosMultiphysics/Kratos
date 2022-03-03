// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

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
#include "containers/model.h"
#include "includes/model_part.h"

// ==============================================================================

namespace Kratos
{

class OptimizationAlgorithm
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(OptimizationAlgorithm);

    OptimizationAlgorithm(std::string OptName, std::string OptType, Model& rModel, Parameters OptSettings)
    : mOptName(OptName), mOptType(OptType), mrModel(rModel), mSettings(OptSettings)
    {
    }

    virtual ~OptimizationAlgorithm() {};



    std::string mOptName;
    std::string mOptType;
    Model& mrModel;
    Parameters mSettings;

}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // ALGORITHM_BASE_H
