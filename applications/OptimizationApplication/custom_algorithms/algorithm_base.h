//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

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

    OptimizationAlgorithm(std::string OptName, std::string OptType, Model& rModel, Parameters& OptSettings)
    : mOptName(OptName), mOptType(OptType), mrModel(rModel), mrSettings(OptSettings)
    {
    }

    virtual ~OptimizationAlgorithm() {};

    // --------------------------------------------------------------------------
    virtual void Initialize(){};
    // --------------------------------------------------------------------------
    virtual void CalculateSolutionStep(){};    

    std::string mOptName;
    std::string mOptType;
    Model& mrModel;
    Parameters& mrSettings;

}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // ALGORITHM_BASE_H
