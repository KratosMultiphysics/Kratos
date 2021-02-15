//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes

// Application includes

// Include base h
#include "adjoint_test_utilities.h"

namespace Kratos
{
void AdjointTestUtilities::PerturbNodalVariable(
    NodeType& rNode,
    const Variable<double>& rVariable,
    const double Perturbation)
{
    if (rVariable.IsComponent() && rVariable.GetSourceVariable() == SHAPE_SENSITIVITY) {
        rNode.Coordinates()[rVariable.GetComponentIndex()] += Perturbation;
    } else {
        rNode.FastGetSolutionStepValue(rVariable) += Perturbation;
    }
}

} // namespace Kratos