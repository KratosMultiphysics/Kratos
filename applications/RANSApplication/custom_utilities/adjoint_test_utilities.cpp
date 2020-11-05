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
#include <functional>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/model_part.h"

// Application includes

// Include base h
#include "adjoint_test_utilities.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{

template<>
std::function<double&(ModelPart::NodeType&, const IndexType)> GetPerturbationMethod(
    const Variable<double>& rPerturbationVariable)
{
    auto perturbation_method = [&rPerturbationVariable](ModelPart::NodeType& rNode,
                                                        const IndexType Dimension) -> double& {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(Dimension != 0)
            << "Dimension should be always 0 for scalar perturbation "
               "variables. [ Perturbation variable = "
            << rPerturbationVariable.Name() << ", Dimension = " << Dimension << " ].\n";

        return rNode.FastGetSolutionStepValue(rPerturbationVariable);

        KRATOS_CATCH("");
    };

    return perturbation_method;
}

template<>
std::function<double&(ModelPart::NodeType&, const IndexType)> GetPerturbationMethod(
    const Variable<array_1d<double, 3>>& rPerturbationVariable)
{
    if (rPerturbationVariable == SHAPE_SENSITIVITY)
    {
        auto perturbation_method = [](ModelPart::NodeType& rNode, const IndexType iDim) -> double& {
            array_1d<double, 3>& r_coordinates = rNode.Coordinates();
            return r_coordinates[iDim];
        };
        return perturbation_method;
    }
    else
    {
        auto perturbation_method = [&rPerturbationVariable](ModelPart::NodeType& rNode,
                                                            const IndexType Dimension) -> double& {
            KRATOS_TRY

            KRATOS_DEBUG_ERROR_IF(Dimension > 2)
                << "Dimension should be always less than 3 for scalar "
                   "perturbation "
                   "variables. [ Perturbation variable = "
                << rPerturbationVariable.Name() << ", Dimension = " << Dimension << " ].\n";

            array_1d<double, 3>& r_vector =
                rNode.FastGetSolutionStepValue(rPerturbationVariable);
            return r_vector[Dimension];

            KRATOS_CATCH("");
        };
        return perturbation_method;
    }
}

template<>
void CompareValues(const double& rA, const double& rB, const double Tolerance)
{
    KRATOS_CHECK_RELATIVE_NEAR(rA, rB, Tolerance);
}

template<>
void CompareValues(const array_1d<double, 3>& rA, const array_1d<double, 3>& rB, const double Tolerance)
{
    KRATOS_CHECK_VECTOR_NEAR(rA, rB, Tolerance);
}

template<>
IndexType GetVariableDimension(const Variable<double>& rVariable, const ProcessInfo& rProcessInfo)
{
    return 1;
}

template<>
IndexType GetVariableDimension(const Variable<array_1d<double, 3>>& rVariable, const ProcessInfo& rProcessInfo)
{
    return rProcessInfo[DOMAIN_SIZE];
}

} // namespace RansApplicationTestUtilities
} // namespace Kratos