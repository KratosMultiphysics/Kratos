// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef STEEPEST_DESCENT_H
#define STEEPEST_DESCENT_H

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
#include "algorithm_base.h"

// ==============================================================================

namespace Kratos
{

class KRATOS_API(OPTIMIZATION_APPLICATION) AlgorithmSteepestDescent : public OptimizationAlgorithm
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(AlgorithmSteepestDescent);

    AlgorithmSteepestDescent(std::string OptName, Model& rModel, Parameters& rOptSettings)
    : OptimizationAlgorithm(OptName,"steepest_desceent",rModel,rOptSettings)
    {
    }

    virtual ~AlgorithmSteepestDescent() {};

    // --------------------------------------------------------------------------
    void Initialize() override {

        std::cout<<"Hi rEza you called me "<<std::endl;
        std::cout<<mrSettings<<std::endl;
    };


}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // STEEPEST_DESCENT_H
