//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_SHALLOW_WATER_TESTS_UTILITIES_H_INCLUDED)
#define KRATOS_SHALLOW_WATER_TESTS_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos {

namespace Testing {

namespace ShallowWaterTestsUtilities {

void AssembleRHS(
    Vector& rRHS_element,
    const Vector& rRHS_condition,
    const std::vector<IndexType>& rIds);

void AddVariables(ModelPart& rModelPart);

void CreateGeometry(
    ModelPart& rModelPart,
    const std::string& rElementName,
    const std::string& rConditionName);

void CalculateAndAssembleRHS(
    ModelPart& rModelPart,
    Vector& rRHS);

} // namespace ShallowWaterTestsUtilities

} // namespace Testing

} // namespace Kratos

#endif // KRATOS_SHALLOW_WATER_TESTS_UTILITIES_H_INCLUDED  defined
